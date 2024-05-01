from scipy.sparse import coo_array
import dask.bag as db
import re
import os
import pandas as pd
import io
import yaml
from .format.meta import FormatMeta

# ASSUMPTIONS:
# GRID is never used
# LP is type dependent; have to ask
# N is constant ; good assumption by and large
# Molecule ID is useless ; have to ask


# Exceptions


class InvalidItem(Exception):
    pass


class InvalidFormat(Exception):
    pass


class InvalidChunks(Exception):
    pass


# TODO
#
# Write intermediate bag step formas
# Overridable blocksize
class Simulation:
    """
    Representation of an entire simulation run, with its metadata.

    The Simulation class is initialized with a metadata file describing
    the parameters of a given simulation's chosen chunks. Trajectories
    and other outputs are stored as attributes (default: None) which
    can be initialized through the appropriate methods.

    Args:
        meta_file (os.PathLike): File specifying simulation metadata
        chunks (list[int]): Chosen simulation chunks to handle
        block (str): Block size for dask
        eager (bool): Whether to compute attributes immediately

    Attributes:
        trajectory: Atomic trajectories
        bonds: ReaxFF bond data
        species: ReaxFF species data
        thermo: Thermodyanmic data

    """

    def __init__(
        self,
        meta_file: os.PathLike,
        chunks: list[int] = None,
        block_size: str = "250MiB",
        eager: bool = False,
    ) -> None:

        with open(meta_file, "r") as f:
            self.meta = FormatMeta(**yaml.safe_load(f)["Metadata"]).model_dump()

        valid_chunks = [i for i in range(self.meta["partition"]["n_chunks"])]

        self.eager = eager
        self.block = block_size

        if chunks is not None:
            if set(chunks) <= set(valid_chunks):
                self.chunks = chunks
            else:
                raise InvalidChunks(
                    f"Valid chunks are in range [0,...,{valid_chunks[-1]}]. Was supplied {chunks}"
                )
        else:
            self.chunks = valid_chunks

        self.trajectory = None
        self.bonds = None
        self.species = None
        self.thermo = None
        self.other_data = (
            {k: None for k in self.meta["output"]["other"].keys()}
            if self.meta["output"]["other"]
            else None
        )

    # End user methods

    def read_bonds(self, data_path: os.PathLike = None):
        """
        Read bond files, parse and store in class attribute

        Args:
        data_path (os.PathLike): alternate base path containing chosen chunks
        """
        base_path = self.meta["data_path"] if data_path is None else data_path
        file_paths = [
            os.path.join(
                base_path,
                f"{chunk}/dat_bonds_{self.meta['sim_id']}_{chunk}.reaxff",
            )
            for chunk in self.chunks
        ]

        corpus = (
            db.read_text(file_paths, linedelimiter="# Timestep", blocksize=self.block)
            .remove(lambda x: x == "# Timestep")
            .map(
                lambda x: [
                    line
                    for line in x.split("\n")
                    if not (line.startswith("#") or line == "")
                ]
            )
            .remove(lambda x: x == [])
            .map(self.__process_bond_step)
            .distinct(key=lambda x: x["timestep"])
        )

        self.bonds = corpus.compute() if self.eager else corpus

    def read_trajectory(
        self, data_path: os.PathLike = None, atomic_format: str = "frame"
    ):
        """
        Read trajectory files, parse and store in class attribute

        Args:
        data_path (os.PathLike): alternate base path containing chosen chunks
        atomic_format: format to project trajectories into
        """
        base_path = self.meta["data_path"] if data_path is None else data_path
        file_paths = [
            os.path.join(
                base_path,
                f"{chunk}/dat_trajectory_{self.meta['sim_id']}_{chunk}.dump",
            )
            for chunk in self.chunks
        ]

        corpus = (
            db.read_text(file_paths, linedelimiter="TIMESTEP", blocksize=self.block)
            .remove(lambda x: x == "ITEM: TIMESTEP")
            .map(lambda x: x.split("ITEM: "))
            .map(lambda x: x[:-1] if (x[-1] == "TIMESTEP") else x)
            .map(self.__process_traj_step, atomic_format=atomic_format)
            # .distinct(key=lambda x: x["timestep"]) ; causes memory leak on nanosecond scale data
        )

        self.trajectory = corpus.compute() if self.eager else corpus

    def read_species(self):
        pass

    def read_ave(self):
        pass

    def read_thermo(self):
        pass

    # Intermediate processing steps

    def __process_traj_step(self, step_text: str, atomic_format: str):
        """
        Parse raw trajectory data text of one frame into chosen format
        """
        frame = {"timestep": "", "n_atoms": "", "atomic": ""}
        item_regex = "([A-Z ]*)([A-z ]*)\n((.*[\n]?)*)"
        valid_items = ["NUMBER OF ATOMS", "BOX BOUNDS", "ATOMS", "DIMENSIONS"]
        # valid_bounds = []

        timestep = int(step_text.pop(0).strip())
        frame["timestep"] = timestep

        for item in step_text:
            label, header, data = re.search(item_regex, item).group(1, 2, 3)
            label = label.strip()

            if label not in valid_items:
                raise InvalidItem("Not a valid LAMMPS data item")

            elif label == "NUMBER OF ATOMS":
                frame["n_atoms"] = int(data.strip())

            elif label == "BOX BOUNDS":
                frame["box"] = {
                    # "dim": int(len(header.split())),
                    "style": header.split(),
                    "bounds": [
                        tuple(map(lambda x: float(x), i.split()))
                        for i in data.split("\n")[:-1]
                    ],
                }

            elif label == "ATOMS":

                # Almost intentionally bad implementation, never use this
                if atomic_format == "frame":
                    frame["atomic"] = [
                        [
                            (
                                float(num)
                                if any(mark in num for mark in [".", "e"])
                                else int(num)
                            )
                            for num in x.split()
                        ]
                        for x in data.split("\n")
                    ]

                elif atomic_format == "pandas":
                    dataf = pd.read_csv(
                        io.StringIO(data), sep=" ", names=header.split()
                    ).set_index("id")
                    frame["atomic"] = dataf

                else:
                    raise InvalidFormat("Select a valid atomic output format")

            elif label == "DIMENSIONS":
                # ??Grid??
                pass

        return frame

    def __process_bond_step(self, step_text: str):
        """
        Parse raw bond data text of one frame into chosen format
        """
        # TODO Leverage symmetry and halve time
        # atomids start from 0 (actualid -1)
        timestep = int(step_text.pop(0))

        i, j, v = [], [], []

        for line in step_text:
            line = line.split()[:-3]
            i += [int(line[0]) - 1] * int(line[2])
            j += [int(x) - 1 for x in line[3 : 3 + int(line[2])]]
            v += [float(val) for val in line[int(line[2]) + 4 :]]

        return {
            "timestep": timestep,
            "bonds": coo_array(
                (v, (i, j)),
                shape=(self.meta["box"]["n_atoms"], self.meta["box"]["n_atoms"]),
            ),
        }

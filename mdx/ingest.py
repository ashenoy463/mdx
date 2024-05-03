from scipy.sparse import coo_array
import dask.bag as db
import re
import os
import pandas as pd
import io
import yaml
from mdx.models.core import check_path
from mdx.models.meta import FormatMeta
from pydantic import PositiveInt, ValidationError, validate_call
from typing import Union


# Exceptions


class InvalidItem(Exception):
    pass


class InvalidFormat(Exception):
    pass


class InvalidChunks(Exception):
    pass


class Simulation:
    """
    Representation of an entire simulation run, with its metadata.

    The Simulation class is initialized with a metadata file describing
    the parameters of a given simulation's chosen chunks. Trajectories
    and other outputs are stored as attributes (default: None) which
    can be initialized through the appropriate methods.

    Args:
        meta_file (os.PathLike): File specifying simulation metadata
        chunks (list[int] / int): Chosen simulation chunks to handle
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

        self.eager = eager
        self.block = block_size
        self.chunks = self.__decide_chunks(chunks, self.meta["partition"]["n_chunks"])

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

    def read_trajectory(
        self,
        data_path: os.PathLike = None,
        atomic_format: str = "frame",
        blocksize: str = None,
    ) -> None:
        """
        Read trajectory files, parse and store in class attribute

        Args:
        data_path (os.PathLike): alternate base path containing chosen chunks
        atomic_format: format to project trajectories into
        """
        corpus = (
            db.read_text(
                self.__get_data_files(data_path, "trajectory"),
                linedelimiter="TIMESTEP",
                blocksize=f"{self.block if blocksize is None else blocksize}",
            )
            .remove(lambda x: x == "ITEM: TIMESTEP")
            .map(lambda x: x.split("ITEM: "))
            .map(lambda x: x[:-1] if (x[-1] == "TIMESTEP") else x)
            .map(self.__process_traj_step, atomic_format=atomic_format)
            # .distinct(key=lambda x: x["timestep"]) ; causes memory leak on nanosecond scale data
        )

        self.trajectory = corpus.compute() if self.eager else corpus

    def read_bonds(self, data_path: os.PathLike = None, blocksize: str = None) -> None:
        """
        Read bond files, parse and store in class attribute

        Args:
        data_path (os.PathLike): alternate base path containing chosen chunks
        """
        corpus = (
            db.read_text(
                self.__get_data_files(data_path, "bonds"),
                linedelimiter="# Timestep",
                blocksize=f"{self.block if blocksize is None else blocksize}",
            )
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

    def read_species(
        self, data_path: os.PathLike = None, blocksize: str = None
    ) -> None:
        """
        Read species files, parse and store species data in class attribute

        Args:
        data_path (os.PathLike): alternate base path containing chosen chunks
        """
        corpus = (
            db.read_text(
                self.__get_data_files(data_path, "species"),
                linedelimiter="# Timestep",
                blocksize=f"{self.block if blocksize is None else blocksize}",
            )
            .map(lambda x: x[1:].split("\n")[:-1])
            .remove(lambda x: x == [])
            .map(self.__process_species_step)
        )

        self.species = corpus.compute() if self.eager else corpus

    def read_thermo(self, data_path: os.PathLike = None) -> None:
        """
        Read log files, parse and store thermodynamic data in class attribute

        Args:
        data_path (os.PathLike): alternate base path containing chosen chunks
        """
        thermo_data = pd.DataFrame()

        for log_path in self.__get_data_files(data_path, "log"):
            with open(log_path, "r", encoding="UTF-8") as logfile:
                # Get relevant part of log
                corpus = logfile.read().split("Loop")[0].split("Step")[1]
                # Extract header
                headers = corpus.split("\n")[0].split()
                headers.insert(0, "Step")
                # Generate and write list of row-lists
                rows = [
                    [float(itm) for itm in line.split()]
                    for line in corpus.split("\n")[1:]
                ][:-1]
                thermo_data = pd.concat(
                    [thermo_data, pd.DataFrame(rows, columns=headers)]
                )

        thermo_data["Boxtime"] = thermo_data["Step"].map(
            lambda x: x * self.meta["partition"]["step_size"]
        )
        thermo_data.drop_duplicates(["Step"], inplace=True)
        thermo_data.reset_index(drop=True, inplace=True)

        self.thermo = thermo_data

    def read_ave(self):
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

    def __process_species_step(self, step_text: str):
        frame = {"timestep": "", "no_moles": "", "no_species": "", "species": {}}
        header, data = list(map(lambda x: x.split(), step_text))

        frame["timestep"] = int(data.pop(0))
        frame["no_moles"] = int(data[0])
        frame["no_species"] = int(data[1])

        for specie, amount in zip(header[2:], data[2:]):
            frame["species"][specie] = int(amount)

        return frame

    # Helper methods

    @validate_call
    def __decide_chunks(
        self, given_chunks: Union[list[int], int, None], meta_chunks: PositiveInt
    ) -> list[int]:
        """
        Validate given chunks against metadata
        """
        valid_chunks = [i for i in range(meta_chunks)]
        if given_chunks is not None:
            if type(given_chunks) is int:
                given_chunks = [given_chunks]
            if set(given_chunks) <= set(valid_chunks):
                return given_chunks
            else:
                raise InvalidChunks(
                    f"Valid chunks are in range [0,...,{valid_chunks[-1]}]. Was supplied {given_chunks}"
                )
                # raise ValidationError()
        else:
            return valid_chunks

    @validate_call
    def __get_data_files(
        self, data_path: Union[None, os.PathLike], type: str, exts=None
    ):
        """
        Get files across simulation chunks
        """
        base_path = self.meta["data_path"] if data_path is None else data_path

        # TEMPFIX: file prefixes and extensions
        filetypes = (
            {
                "trajectory": ("dat_trajectory" ".dump"),
                "bonds": ("dat_bonds", ".reaxff"),
                "species": ("dat_species", ".out"),
                "log": ("log_out", ""),
            }
            if exts is None
            else exts
        )

        file_paths = [
            check_path(
                os.path.join(
                    base_path,
                    f"{chunk}/{filetypes[type][0]}_{self.meta['sim_id']}_{chunk}{filetypes[type][1]}",
                )
            )
            for chunk in self.chunks
        ]
        return file_paths

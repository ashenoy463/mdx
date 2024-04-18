import dask
from scipy.sparse import coo_array
import dask.bag as db
import dask.dataframe as dd
import dask.array as da
from toolz.functoolz import compose
import re
import os
import numpy as np
import pandas as pd
import io
import yaml

# ASSUMPTIONS:
#
# GRID is never used
# LP is type dependent
# N is constant
# Molecule ID is useless


class InvalidItem(Exception):
    pass


class Ingest:
    # filereading methods
    # specialchunks
    # better metafiles
    # handle invalid glob error
    # client management
    # classdocstring
    # record traj resolution
    # write comments at each step
    # i/o interfaces
    # bond rehydration should respect groups
    # UTEST
    def __init__(self, meta_file, chunks=[], block="250MiB", eager=False) -> None:

        self.chunks = chunks
        self.eager = eager
        self.block = block

        with open(meta_file, "r") as f:
            self.meta = yaml.safe_load(f)["Metadata"]

        if not chunks:
            self.chunks = [i for i in range(self.meta["partition"]["n_chunks"])]

    # End user methods

    def read_bonds(self, bond_path="."):

        bond_paths = [
            os.path.join(
                bond_path,
                f"{self.meta['sim_id']}/{chunk}/dat_bonds_{self.meta['sim_id']}_{chunk}.reaxff",
            )
            for chunk in self.chunks
        ]
        corpus = db.read_text(
            bond_paths, linedelimiter="# Timestep", blocksize=self.block
        )
        corpus = corpus.remove(lambda x: x == "# Timestep")
        corpus = corpus.map(
            lambda x: [
                line
                for line in x.split("\n")
                if not (line.startswith("#") or line == "")
            ]
        )
        corpus = corpus.remove(lambda x: x == [])
        corpus = corpus.map(self.__process_bond_step)
        result = corpus.distinct(key=lambda x: x["timestep"])

        return result

    def read_species(self, spec_path="."):
        # SERIAL
        spec_paths = [
            os.path.join(
                spec_path,
                f"{self.meta['sim_id']}/{chunk}/dat_species_{self.meta['sim_id']}_{chunk}.out",
            )
            for chunk in self.chunks
        ]
        corpus = list(map(lambda f: open(f).read(), spec_paths))
        return corpus.take(5)

    def read_ave(self):
        pass

    def read_log(self, log_path="."):
        log_glob = os.path.join(".", f"log_out_prelim_*")

    def read_trajectory(self, traj_path=".", atomic_format="frame"):
        traj_paths = [
            os.path.join(
                traj_path,
                f"{self.meta['sim_id']}/{chunk}/dat_trajectory_{self.meta['sim_id']}_{chunk}.dump",
            )
            for chunk in self.chunks
        ]
        corpus = db.read_text(
            traj_paths, linedelimiter="TIMESTEP", blocksize=self.block
        )
        steps = corpus.remove(lambda x: x == "ITEM: TIMESTEP")
        items = steps.str.split("ITEM: ")
        items = items.map(lambda x: x[:-1] if (x[-1] == "TIMESTEP") else x)
        result = items.map(self.__process_traj_step, atomic_format=atomic_format)
        # Drop duplicate frames across chunks
        result = result.distinct(key=lambda x: x["timestep"])

        if self.eager:
            return result.compute()
        else:
            return result

    # Intermediate processing steps

    def __process_traj_step(self, step_text, atomic_format):
        """
        Form a frame from list of items in a trajectory timestep
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
                # Atomwise data
                if atomic_format == "frame":
                    frame["atomic"] = [x for x in data.split("\n")]
                elif atomic_format == "pandas":
                    dataf = pd.read_csv(
                        io.StringIO(data), sep=" ", names=header.split()
                    ).set_index("id")
                    frame["atomic"] = dd.from_pandas(dataf, npartitions=1)
                    frame["atomic"] = dataf

            elif label == "DIMENSIONS":
                # ??Grid??
                pass

        return frame

    def __process_bond_step(self, step_text):

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


def extract(f, bag):
    """
    Eagerly map function onto bag and return an array
    """
    array = np.array((bag.map(f).compute()))
    return array

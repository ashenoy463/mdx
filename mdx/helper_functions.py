import os
import h5py
from collections.abc import MutableMapping
import numpy as np


# thints
class DataFileExists(Exception):
    pass


def dict_merge(dct, merge_dct):
    """
    Recursive dict merge. Inspired by :meth:``dict.update()``, instead of
    updating only top-level keys, dict_merge recurses down into dicts nested
    to an arbitrary depth, updating keys. The ``merge_dct`` is merged into
    ``dct``.
    :param dct: dict onto which the merge is executed
    :param merge_dct: dct merged into dct
    :return: None
    """
    # From https://gist.github.com/angstwad/bf22d1822c38a92ec0a9
    for k, v in merge_dct.items():
        if (
            k in dct and isinstance(dct[k], dict) and isinstance(merge_dct[k], dict)
        ):  # noqa
            dict_merge(dct[k], merge_dct[k])
        else:
            dct[k] = merge_dct[k]


def dict_flatten(dictionary, parent_key="", separator="_"):
    # From https://stackoverflow.com/a/6027615
    items = []
    for key, value in dictionary.items():
        new_key = parent_key + separator + key if parent_key else key
        if isinstance(value, MutableMapping):
            items.extend(dict_flatten(value, new_key, separator=separator).items())
        else:
            items.append((new_key, value))
    return dict(items)


def write_rdf(rdf_obj, sim_id: str, chunk_id: str, c: str, s: str, suffix=None):
    """
    Store MDAnalysis InterRDF results to HDF5 file

    rdf_obj:    MDAnalaysis InterRDF object with RDF data
    sim_id:     Unique simulation ID
    chunk_id:   Chunk ID
    c:          Central atom name
    s:          Surrounding atom name

    suffix:     Optional suffix to standard filename
    """
    # TODO
    # combine chunk data and remove param
    # batched creation
    # path handling
    # default chunkid?
    # typehints
    # write some metadata
    #
    filename = f"./rdf_{sim_id}_{chunk_id}_{c}_{s}.hdf5"
    print(filename)
    if suffix:
        filename = filename[:-5] + f"_{suffix}" + ".hdf5"
    if not os.path.isfile(filename):
        with h5py.File(filename, "w") as f:
            f.create_dataset("rdf", data=rdf_obj.results.rdf, compression="gzip")
            f.create_dataset("bins", data=rdf_obj.results.bins, compression="gzip")
            f.create_dataset("edges", data=rdf_obj.results.edges, compression="gzip")
            f.create_dataset("count", data=rdf_obj.results.count, compression="gzip")
    else:
        raise DataFileExists(
            f"RDF dataset for {sim_id}_{chunk_id} {c}-{s} already exists"
        )


def read_rdf(sim_id: str, chunk_id: str, c: str, s: str):
    """
    sim_id:     Unique simulation ID
    chunk_id:   Chunk ID
    c:          Central atom name
    s:          Surrounding atom name
    """
    filename = f"./rdf_{sim_id}_{chunk_id}_{c}_{s}.hdf5"
    try:
        with h5py.File(filename, "r") as f:
            rdf = f.get("rdf")[:]
            bins = f.get("bins")[:]
            edges = f.get("edges")[:]
            count = f.get("count")[:]
    except FileNotFoundError:
        raise Exception(f"RDF dataset for {sim_id}_{chunk_id} {c}-{s} not found")
    return {"rdf": rdf, "bins": bins, "edges": edges, "count": count}


# def bond_unpack()
# sym
# filltocorners


def extract(f, bag):
    """
    Eagerly map function onto bag and return an array
    """
    array = np.array((bag.map(f).compute()))
    return array

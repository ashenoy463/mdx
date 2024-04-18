import os
import h5py
import dask


# thints
class DataFileExists(Exception):
    pass


def bundle_data(label):
    pass


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
    #sym
    #filltocorners
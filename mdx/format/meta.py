from pydantic import BaseModel, PositiveFloat, PositiveInt, ValidationError
from typing import Literal, Dict, Any, Optional, Union
from datetime import datetime

import os
from pydantic.functional_validators import AfterValidator
from typing_extensions import Annotated

# Allowed values

# fmt: off
valid_elements = [
"H","He",
"Li","Be","B","C","N","O","F","Ne",
"Na","Mg","Al","Si","P","S","Cl","Ar",
"K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr",
"Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe",
"Cs","Ba",
"La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb",
"Lu","Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn",
"Fr","Ra","Ac","Th","Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No",
"Lr","Rf","Db","Sg","Bh","Hs","Mt","Ds","Rg","Cn","Nh","Fl","Mc","Lv","Ts","Og"
]
# fmt: on

# Validators


def check_path(path: os.PathLike) -> os.PathLike:
    assert os.path.exists(path), f"{path} is not a valid path"
    return path


ValidPath = Annotated[os.PathLike, AfterValidator(check_path)]

# Format for simulation metadata


class MetaPartition(BaseModel):
    step_size: PositiveFloat  # Femtoseconds in a step
    chunk_size: PositiveInt  # Steps in a chunk
    n_chunks: PositiveInt  # Chunks in a simulation


# TEMPFIX: Everything that would be here is left to the user currently
class MetaExperimental(BaseModel):
    pass


class MetaOutput(BaseModel):
    # Interval of steps at which output files are written
    trajectory: PositiveInt
    bonds: PositiveInt
    species: PositiveInt
    thermo: PositiveInt
    other: Dict[str, PositiveInt] = None


class MetaBox(BaseModel):
    n_atoms: PositiveInt
    elements: list[tuple[int, Literal[tuple(valid_elements)]]]  # type: ignore


class FormatMeta(BaseModel):
    sim_id: str
    sim_desc: str
    exec_times: list[tuple[datetime, datetime]]
    data_path: ValidPath
    partition: MetaPartition
    box: MetaBox
    experimental: Optional[Dict[str, Any]]
    output: MetaOutput

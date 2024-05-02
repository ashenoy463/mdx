from pydantic.functional_validators import AfterValidator
from typing_extensions import Annotated
import os

# ValidPath


def check_path(path: os.PathLike) -> os.PathLike:
    assert os.path.exists(path), f"{path} is not a valid path"
    return path


ValidPath = Annotated[os.PathLike, AfterValidator(check_path)]

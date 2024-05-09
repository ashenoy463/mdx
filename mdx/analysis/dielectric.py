from scipy.constants import epsilon_0, Boltzmann, e, angstrom
import numpy as np
import xarray as xr
from . import norm


def dipole_moment(ds: xr.Dataset) -> xr.Dataset:
    """
    Compute total dipole moment about box origin

    !!SI Units!!

    Required Dimensions:
        r: Atomic position (...,atom,3)
        q: Atomic charge (...,atom,1)
        atom: Atom index

    Args:
        ds (xr.Dataset): Trajectory dataset

    Returns:
        xr.Dataset: Dipole moment vectors (...,3)
    """
    return (ds["q"] * e * ds["r"] * angstrom).sum(dim="atom")


def box_volume(ds: xr.Dataset):
    """
    Compute box volume from bounds

    !!SI Units!!

    Required Dimensions:
        box: Box bounds (...,pos,case)
        case: Upper/lower bound (2)
        pos: Coordinate direction (3)

    Args:
        ds (xr.Dataset): Trajectory dataset

    Returns:
        xr.Dataset: Box volumes (...,1)
    """
    return ds["box"].diff("case").prod(dim="pos") * (angstrom**3)


def dielectric_const(ds: xr.Dataset, T: float) -> xr.Dataset:
    """
    Compute dielectric constant from dipole moment fluctuations
    as per Gereben 2011

    !!SI Units!!

    Formula:
        epsilon = [ <M^2> - <M>^2 ] / [ 3*epsilon_0*V*k_B*T ]

    Required Dimensions:
        m: Box bounds (step,3)
        volume: Upper/lower bound (step,1)
        step: Timestep

    Args:
        ds (xr.Dataset): Trajectory dataset
        T (float): Temperature

    Returns:
        xr.Dataset: Dielectric constants (1)
    """
    return 1 + (
        norm(ds["m"]).var(dim="step")
        / (3 * epsilon_0 * ds["volume"].mean(dim="step") * Boltzmann * T)
    )

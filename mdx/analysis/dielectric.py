from scipy.constants import epsilon_0, Boltzmann, e, angstrom
import numpy as np


def compute_eps(M: np.ndarray, V: float, T: float) -> float:
    """
    Compute dielectric constant from dipole moment fluctuations
    as per Gereben 2011

    M : Dipole moment vectors  (n,3)
    V : Box Volume             (1,1)
    T : Temperature            (1,1)

    Formula:
    epsilon = [ <M^2> - <M>^2 ] / [ 3*epsilon_0*V*k_B*T ]

    Returns: Dielectric constant [scalar]
    """
    # <M^2>
    M1 = np.average(np.square(np.linalg.norm(M, axis=1)))
    # <M>^2
    M2 = np.average(np.linalg.norm(M, axis=1)) ** 2
    # Gereben 2011 , eqn. 1
    epsilon = 1 + ((M1 - M2) / (3 * epsilon_0 * V * Boltzmann * T))

    return epsilon


def compute_volume(frame):
    """
    Compute the volume of the box in a frame
    !! SI Units !!

    Returns: Volume [scalar]
    """
    box_dims = np.array(frame["box"]["bounds"])

    return np.prod((box_dims[:, 1] - box_dims[:, 0]), axis=0) * (angstrom**3)


def compute_dipole(frame) -> np.ndarray:
    """
    Compute the net dipole moment of the box in a frame
    !! SI Units !!

    Returns: Dipole moment vector [ndarray of shape (1,3)]
    """
    traj = frame["atomic"]
    r = np.array(traj[["xu", "yu", "zu"]].values) * angstrom
    qr = np.multiply(traj["q"].values * e, r.T).T
    M = qr.sum(axis=0)

    return M

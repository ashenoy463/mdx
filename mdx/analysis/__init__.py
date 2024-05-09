__all__ = [
    "dielectric",
    # 'kinetics',
    # 'rdf'
]

norm = lambda x: np.sqrt(x[..., 0] ** 2 + x[..., 1] ** 2 + x[..., 2] ** 2)

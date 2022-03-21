
# TODO(Alex)
# This is more work than I would like:
# Need to ensure the cell origin is at the corner
# Would prefer fractional coordinates
# Need to test that the expression works on a simple 2D example I can work out the result to
# dr - np.rint(dr / cell_lengths) * cell_lengths

# def minimum_image_distance_matrix(a: np.ndarray, lattice: np.ndarray):
#
#     assert a.shape[1] == 3, "Expect vectors to be Euclidean"
#     n_vectors = a.shape[0]
#     cell_lengths = np.asarray([np.linalg.norm(lattice[i, :]) for i in range(0, 3)])
#
#     d = np.empty(shape=(n_vectors, n_vectors))
#     np.fill_diagonal(d, 0.)
#
#     # Upper triangle
#     for i in range(0, n_vectors):
#         for j in range(i + 1, n_vectors):
#             dr = a[j, :] - a[i, :]
#             print(dr - np.rint(dr / cell_lengths) * cell_lengths)
#             d[i, j] = np.linalg.norm(dr - np.rint(dr / cell_lengths) * cell_lengths)
#             d[j, i] = d[i, j]
#
#     return d
#
#
# def distance_matrix(a: np.ndarray):
#
#     assert a.shape[1] == 3, "Expect vectors to be Euclidean"
#     n_vectors = a.shape[0]
#     d = np.empty(shape=(n_vectors, n_vectors))
#     np.fill_diagonal(d, 0.)
#
#     # Upper triangle
#     for i in range(0, n_vectors):
#         for j in range(i + 1, n_vectors):
#             d[i, j] = np.linalg.norm(a[j, :] - a[i, :])
#             d[j, i] = d[i, j]
#
#     return d
#
# def optimal_muffin_tin_radii(positions: np.ndarray, lattice: np.ndarray, atomic_numbers: List[int]):
#
#     # Use rgkmax as proxy for MT radius (seem to be correlated)
#     species_min_mt = np.amin([fixed_precision_rgkmax(x) for x in atomic_numbers])
#
#     dm = scipy_distance_matrix(positions, positions)
#     dm2 = distance_matrix(positions)
#     dm3 = minimum_image_distance_matrix(positions - positions[0, :], lattice)
#     print(dm - dm2)
#     print(dm2)
#     print(dm3)
#

# def mufftin_radius_to_conserve_precision():
#     """ Given the minimum MT radius in the system, and the input RGKMAX,
#     return the MT radius to assure a consistent RGKMAX for element x.
#
#     :return:
#     """

# def consistent_muffin_tin_radius(rgkmax_input: int, atomic_number_x: int, minimum_mt_radius: float):
#     """ Given the smallest MT radius of species in a system, return the MT radius of element X
#
#     Different `rgkmax` values are required for different species to obtain the same accuracy (in same quantity - I
#     assume total energy). Because `rgkmax` is set once according to the smallest muffin tin radius, this determines
#     the maximum $G$ vector for all species. As such, to obtain the desired `rgkmax` values for all other species in the
#     calculation, their muffin tin radii must be set proportionally to min(MT), to obtain the target `rgkmax` values
#     which ensure a consistent accuracy.
#
#      \textrm{rgkmax}_{input} = \min(\textrm{MT}) * \textrm{G}_{max}
#
#      therefore:
#
#      \begin{equation}
#      \textrm{G}_{max} = \textrm{rgkmax}_{input} / \min(MT).
#      \end{equation}
#
#      So for a given species, $X$, in the system:
#
#      \begin{align}
#      \textrm{rgkmax}_X &= MT_X * \textrm{G}_{max},  \\
#      \textrm{rgkmax}_X &= MT_X * \textrm{rgkmax}_{input} / \min(MT).
#      \end{align}
#
#      Therefore, to obtain the desired $\textrm{rgkmax}$, $MT_X$ is the free parameter:
#
#      \begin{equation}
#      MT_X = \frac{\textrm{rgkmax}_X}{\textrm{rgkmax}_{input}} * min(MT)
#      \end{equation}
#
#     :return:
#     """
#     # return (fixed_precision_rgkmax(atomic_number_x) / rgkmax_input) * minimum_mt_radius


# def minimum_muffin_tin_radius(positions: np.ndarray, lattice: np.ndarray, atomic_numbers: List[int]):
#     """
#
#
#
#     :param positions:
#     :param lattice:
#     :param atomic_numbers:
#     :return:
#     """
#     bonds = find_minimum_bond_lengths(positions, lattice, atomic_numbers)
#     an_min = atomic_number_species_with_mt_min(atomic_numbers)
#     rgkmax_min = fixed_precision_rgkmax(an_min)
#
#     for an, bond_length in bonds.items():
#         denominator = 1. + (fixed_precision_rgkmax(an) / rgkmax_min)
#         # Percentage of the bond that should be taken by the radius of an_min and an, respectively
#         percentage_min = 1. / denominator
#         percentage_x = (fixed_precision_rgkmax(an) / rgkmax_min) / denominator
#         # print(percentage_min, percentage_x)
#         print((an_min, an), percentage_min * bond_length)
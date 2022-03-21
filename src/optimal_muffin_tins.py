import numpy as np
from typing import List

from ase.neighborlist import mic
# from ase.geometry import complete_cell, find_mic, wrap_positions

from src.materials import MoS2WS2Bilayer, ZnOWurzite, ZrO2Primitive, SiliconPrimitive, TiO2Rutile


def fixed_precision_rgkmax(atomic_number: int) -> float:
    """ Get an rgkmax value, given an atomic number

    rgkmax that give a consistent total_energy precision, per elemental crystal.
    Found empirically by Sven working on the Delta-DFT project.
    TODO(Alex) Document what precision in total energy, these rgkmax resulted in

    :param atomic_number: Atomic number
    :return: rgkmax value
    """
    fixed_rgkmax = np.array([
        5.835430, 8.226318, 8.450962, 8.307929,
        8.965808, 9.376204, 9.553568, 10.239864, 10.790975, 10.444355, 10.636286, 10.579793,
        10.214125, 10.605334, 10.356352, 9.932381, 10.218153, 10.466519, 10.877475, 10.774763,
        11.580691, 11.800971, 11.919804, 12.261896, 12.424606, 12.571031, 12.693836, 12.781331,
        12.619806, 12.749802, 12.681350, 12.802838, 12.785680, 12.898916, 12.400000, 10.596757,
        11.346060, 10.857573, 11.324413, 11.664200, 11.859519, 11.892673, 12.308470, 12.551024,
        12.740728, 12.879424, 13.027090, 13.080576, 13.230621, 13.450665, 13.495632, 13.261039,
        13.432654, 11.329591, 13.343047, 13.011835, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan,
        np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, 15.134859, 14.955721,
        14.607311, 13.930505, 13.645267, 13.629439, 13.450805, 13.069046, 13.226699, 13.261342,
        13.365992, 13.557571, 13.565048, 13.579543, np.nan, 12.273924])
    return fixed_rgkmax[atomic_number + 1]


def atomic_number_species_with_mt_min(atomic_numbers: List[int]) -> int:
    """ Get atomic number of the species expected to have the smallest MT radius in the system.

    Use rgkmax as proxy for MT radius (seem to be correlated).

    :param atomic_numbers: List of atomic numbers
    :return: Atomic number of element expected to have the smallest MT radus.
    """
    index_rgkmax_min = np.argmin([fixed_precision_rgkmax(x) for x in atomic_numbers])
    an_min_mt = atomic_numbers[index_rgkmax_min]
    return an_min_mt


def find_minimum_bond_lengths(positions: np.ndarray, lattice: np.ndarray, atomic_numbers: List[int]) -> dict:
    """Find the minimum bond length between the species with (the anticipated) smallest MT radius,
    and each other species in a periodic cell.

    :return minimum_bond_lengths: Minimum bond length between species X and Y, for all Y, where X
    is the species with the smallest MT radius.
    """
    np_atomic_numbers = np.asarray(atomic_numbers)
    unique_atomic_numbers = set(atomic_numbers)
    an_min_mt = atomic_number_species_with_mt_min(atomic_numbers)

    # Get the atomic indices for each species
    indices = {}
    for an in unique_atomic_numbers:
        indices[an] = np.where(np_atomic_numbers == an)[0]

    # Find all pair vectors in the unit cell, between two species (X_min_MT, Y)
    # and store the shortest corresponding bond_length.
    minimum_bond_lengths = {}
    r_x = positions[indices[an_min_mt]]
    for an in unique_atomic_numbers:
        r_y = positions[indices[an]]
        d_vectors = wrapped_displacement_vectors(r_x, r_y, lattice, remove_self_interaction=True)
        minimum_bond_lengths[an] = minimum_vector_length(d_vectors)

    return minimum_bond_lengths


def wrapped_displacement_vectors(r_x: np.ndarray,
                                 r_y: np.ndarray,
                                 lattice: np.ndarray,
                                 remove_self_interaction=True) -> np.ndarray:
    """ Find all displacement vectors between atom/s at position/s r_x
    and atom/s at position/s r_y, in the minimum image convention.

    wrapped_vectors = [[r0 - r0], [r1 - r0], ... [r_Ny - r_Nx]]

    :return wrapped_vectors: Displacement vectors between positions r_x and r_y, in the
     minimum image convention.
    """
    displacement_vectors = []
    for i in range(r_x.shape[0]):
        for j in range(r_y.shape[0]):
            displacement_vectors.append(r_y[j, :] - r_x[i, :])

    if remove_self_interaction:
        zeros = (displacement_vectors == np.array([0., 0., 0.])).all(-1)
        non_zeros = [not x for x in zeros]
        displacement_vectors = np.asarray(displacement_vectors)
        displacement_vectors = displacement_vectors[non_zeros, :]

    # Apply minimum image convention to these vectors
    # TODO(Alex) Check this works as expected
    wrapped_vectors = mic(displacement_vectors, lattice, pbc=True)
    return wrapped_vectors


def minimum_vector_length(vectors: np.ndarray) -> float:
    """ Find the minimum magnitude from a set of displacement vectors.

    :param vectors: Array of vectors.
    :return: minimum vector length.
    """
    if vectors.shape[1] != 3:
        raise ValueError('Vectors should be stored (n_vectors, 3)')

    n_vectors = vectors.shape[0]
    norms = np.empty(shape=n_vectors)

    for i in range(n_vectors):
        norms[i] = np.linalg.norm(vectors[i, :])

    return np.min(norms)


def optimal_smallest_muffin_tin_radius(an_x, an_y, bond_length):
    """For each bond length, compute the smaller MT radius associated with the two
    elements comprising the bond, according to the ratio of rgkmaxs.

    MT_y = (rgkmax_y / rgkmax_x) * MT_x
    where MT_x corresponds to the smallest MT in the system, and rgkmax_i are a consistent set of tabulated
    rgkmax values (leading to the same precision in total energy), the smallest MT radius can be used
    to determine consistent MT radii for all other elements in the system.

    The sum of two muffin tin radii cannot exceed the bond length:
    MT_y + MY_x = bond_length

    therefore substituting for MT_y and rearranging for MT_x:
    MT_x = bond_length / (1 + (rgkmax_y / rgkmax_x))

    :param an_x:
    :param an_y:
    :param bond_length:
    :return:
    """
    mt_x = bond_length / (1 + fixed_precision_rgkmax(an_y) / fixed_precision_rgkmax(an_x))
    return mt_x


def main(system, scaling_factor=1.0):
    """Get Muffin Tin radii with optimal ratios w.r.t. the smallest MT.

    The smallest MT radius determines all other MT radii in the system.

    This routine proceeds by:
        * Using rgkmax values to establish the species that will have the smallest MT radius in the system.
        * Compute the smallest bond lengths between X (element with smallest RMT) and all Y (other elements).
        * Given these minimum bond lengths, one can use a ratio of rgkmaxs to determine the opitmal ratio of MT radii
          for each bond.
        * MT_X are tabulated for each bond - one selects min(MT_X), and then recomputes MT_Y with this value,
          such that:
         a) MTs do not overlap in any instance
         b) Given min(MT_X), a value of MT_Y is chosen in every case such that the rgkmax ratio is preserved.

    TODO Run this for a few scaling factors, put into the code, check for when the density overlap becomes small

    Basically want to get the touching sphere recommendations, select the smallest MT radius for the element
    expected to have the smallest MT radius, else I cannot maintain valid ratios with the other species without overlap
    Then I just reduce all MT radii by 10-30% (printing out the radii) and do the charge plotting to view an optimal

    :return:
    """
    an_mt_min = atomic_number_species_with_mt_min(system.atomic_numbers)
    min_bond_lengths = find_minimum_bond_lengths(system.positions, system.lattice, system.atomic_numbers)

    for an, bond in min_bond_lengths.items():
        print(f'Minimum bond length between atomic numbers ({an_mt_min}, {an}): {bond}')

    print("For each bond, compute the MT radius associated with the element with the smallest MT radius. \n"
          "Do this for each bond, and choose the minimum.")
    mt_mins = []
    for an_y, bond_length in min_bond_lengths.items():
        radius = optimal_smallest_muffin_tin_radius(an_mt_min, an_y, bond_length)
        mt_mins.append(round(radius, ndigits=4))
    print(mt_mins)

    mt_min = scaling_factor * min(mt_mins)
    print(f"Smallest muffin tin radius for element {an_mt_min} is {mt_min}, scaled by {scaling_factor}")

    print("Use this to determine the MT radii for all other elements in the system, with optimal ratio \n"
          "determined by tabulated rkgmax ratios")
    for an_y in set(system.atomic_numbers):
        mt_y = fixed_precision_rgkmax(an_y) / fixed_precision_rgkmax(an_mt_min) * mt_min
        print(f"Atomic Number {an_y}, MT radius {mt_y}, sum of MT radii {mt_min + mt_y}, min_bond_length {min_bond_lengths[an_y]}")


if __name__ == "__main__":
    #main(MoS2WS2Bilayer(), scaling_factor=1.0)
    main(ZrO2Primitive(), scaling_factor=1.0)

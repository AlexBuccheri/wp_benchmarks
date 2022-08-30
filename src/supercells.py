"""
Silicon supercells for exciting input.

Form a set of 3 benchmarks. 2, 8 and 12 atoms (could consider 16)
"""
import numpy as np
from typing import List

from ase.atoms import Atoms
from ase.build import make_supercell
from ase.io import write

from materials import SiliconPrimitive, angstrom_to_bohr


def primitive_to_supercell(translations: List[int]) -> Atoms:
    """
    Only performs translations along a1, a2 and a3 to create supercells.
    See https://en.wikipedia.org/wiki/Supercell_(crystal)

    :param translations: List of length 3, with values > 0
    :return:
    """
    assert len(translations) == 3, "Must supply 3 integer translations"

    assert SiliconPrimitive.unit == 'angstrom', "silicon primitive cell expected to be tabulated in angstrom"
    primitive_cell = Atoms(cell=SiliconPrimitive.lattice,
                           scaled_positions=SiliconPrimitive.fractional_positions,
                           pbc=[1, 1, 1])
    # write('si_prim.xyz', primitive_cell)

    P_matrix = np.zeros(shape=(3, 3))
    np.fill_diagonal(P_matrix, translations)

    return make_supercell(primitive_cell, P_matrix)


if __name__ == "__main__":

    # Generate supercells by extending translation vectors
    for n1 in [1, 4, 5, 6, 8]:
        super_cell = primitive_to_supercell([n1, 1, 1])
        # write('si_cell.xyz', super_cell)

        fractional_pos = super_cell.get_scaled_positions()
        n_atoms = fractional_pos.shape[0]

        print(f'Supercell containing {n_atoms} atoms')
        print('Lattice vectors in bohr')
        print(super_cell.get_cell() * angstrom_to_bohr)
        print()

        print('Atomic positions in fractional coordinates')
        print(fractional_pos)
        print()


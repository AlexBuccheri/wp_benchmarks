"""
WP2 Benchmarks, written in Python data structures
"""
from dataclasses import dataclass
import numpy as np

from pymatgen.io.cif import CifParser

# Constants
bohr_to_angstrom = 0.529177210903
angstrom_to_bohr = 1. / bohr_to_angstrom


def cif_to_pymatgen(file: str, get_primitive: bool):
    """Extract CIF Data.

    atoms = ase.atoms.Atoms(numbers=structure.atomic_numbers,
                            cell=structure.lattice.matrix,
                            scaled_positions=structure.frac_coords,
                            pbc=True)

    :param file: file name
    """
    structure = CifParser(file).get_structures()[0]

    # TODO(Alex) Does not appear to give the primitive lattice vectors
    if get_primitive:
        structure = structure.get_primitive_structure()

    print(structure)
    print("Atomic numbers: ", structure.atomic_numbers)
    print("Lattice: ", structure.lattice)
    print("Fractional coordinates:", structure.frac_coords)


@dataclass(frozen=True)
class ZrO2Primitive:
    """A1 System

    Ref: https://github.com/nomad-coe/greenX-wp2/tree/main/Benchmarks/Y-ZrO2/PBEsol/A_Cubic/B_Primitive
    Lattice in angstrom, as extracted from FHI-aims input file
    """
    lattice = np.array([[0.00000000,  2.53574055,   2.53574055],
                        [2.53574055,  0.00000000,   2.53574055],
                        [2.53574055,  2.53574055,   0.00000000]])
    fractional_positions = np.array([[0.00000000,  0.00000000,  0.00000000],
                                     [0.25000000,  0.25000000,  0.2500000],
                                     [0.75000000,  0.75000000,  0.7500000]])

    unit = 'angstrom'
    positions = np.transpose(np.transpose(lattice) @ np.transpose(fractional_positions))

    elements = ['Zr', 'O', 'O']
    atomic_numbers = [40, 8, 8]


@dataclass(frozen=True)
class MoS2WS2Bilayer:
    """ Tabulated E1 System.

    Ref: https://github.com/nomad-coe/greenX-wp2/blob/main/Benchmarks/Heterobilayers/MoS2-WS2.xyz
    Values in angstrom, as extracted from .xyz file.
    """
    # lattice = np.array([[3.168394160510246, -3.3524146226426544e-10,  0.0],
    #                     [-1.5841970805453853, 2.7439098312114987, 0.0],
    #                     [6.365167633892244e-18 -3.748609667104484e-30, 39.58711265]])

    # Lattice with ~ zero terms rounded. Values in angstrom
    lattice = np.array([[3.168394160510246,   0.0,                0.0],
                        [-1.5841970805453853, 2.7439098312114987, 0.0],
                        [0.0,                 0.0,                39.58711265]])

    # Positions in angstrom
    positions = np.array([[0.00000000, 0.00000000, 16.68421565],
                          [1.58419708, 0.91463661, 18.25982194],
                          [1.58419708, 0.91463661, 15.10652203],
                          [1.58419708, 0.91463661, 22.90251866],
                          [0.00000000, 0.00000000, 24.46831689],
                          [0.00000000, 0.00000000, 21.33906353]])

    fractional_positions = np.transpose(np.linalg.inv(np.transpose(lattice)) @ np.transpose(positions))

    unit = 'angstrom'
    elements = ['W', 'S', 'S', 'Mo', 'S', 'S']
    atomic_numbers = [74, 16, 16, 42, 16, 16]


@dataclass(frozen=True)
class TiO2Rutile:
    """
    Copied from:
    https://github.com/nomad-coe/greenX-wp2/blob/main/Benchmarks/TiO2-rutile/TiO2.abi

    A materials project reference:
    https://materialsproject.org/materials/mp-2657/
    """
    fractional_positions = np.array([[0.000000000,  0.000000000,  0.000000000],
                                     [0.500000000,  0.500000000,  0.500000000],
                                     [0.303779258,  0.303779258,  0.000000000],
                                     [0.696220742,  0.696220742,  0.000000000],
                                     [0.803779258,  0.196220742,  0.500000000],
                                     [0.196220742,  0.803779258,  0.500000000]])

    lattice = np.array([[8.680645000000,  0.000000000000,  0.000000000000],
                        [0.000000000000,  8.680645000000,  0.000000000000],
                        [0.000000000000,  0.000000000000,  5.591116638050]])
    unit = 'Bohr'
    positions = np.transpose(np.transpose(lattice) @ np.transpose(fractional_positions))

    elements = ['Ti', 'Ti', 'O', 'O', 'O', 'O']
    atomic_numbers = [22, 22, 8, 8, 8, 8]


@dataclass(frozen=True)
class ZnOWurzite:
    """ZnO Wurzite.

    Full Formula (Zn2 O2)
    Reduced Formula: ZnO
    abc   :   3.249400   3.249400   5.203800
    angles:  90.000000  90.000000 120.000000
    Sites (4)
      #  SP           a         b       c
    ---  ----  --------  --------  ------
      0  Zn    0.666667  0.333333  0
      1  Zn    0.333333  0.666667  0.5
      2  O     0.666667  0.333333  0.6179
      3  O     0.333333  0.666667  0.1179

    Vesta is consistent with what's in the CIF - I should just take the result from the CIF

    CIF data extracted using:
      cif_to_pymatgen('/Users/alexanderbuccheri/Python/pycharm_projects/wp_benchmarks/ZnO_WZ/ZnO-Wz.cif')
    """
    # Angstrom (as from cif)
    lattice = np.array([[-1.624700,  -2.814063, -0.000000],
                        [-1.624700,   2.814063,  0.000000],
                        [ 0.000000,   0.000000, -5.203800]])

    fractional_positions = np.array([[6.66666667e-01, 3.33333333e-01, 1.23259516e-32],
                                     [3.33333333e-01, 6.66666667e-01, 5.00000000e-01],
                                     [6.66666667e-01, 3.33333333e-01, 6.17900000e-01],
                                     [3.33333333e-01, 6.66666667e-01, 1.17900000e-01]])

    # In angstrom
    positions = np.transpose(np.transpose(lattice) @ np.transpose(fractional_positions))

    unit = 'angstrom'
    elements = ['Zn', 'Zn', 'O', 'O']
    atomic_numbers = [30, 30, 8, 8]


@dataclass(frozen=True)
class SiliconPrimitive:
    """Silicon primitive  cell

      Full Formula (Si2)
      Reduced Formula: Si
      abc   :   5.430000   5.430000   5.430000
      angles:  90.000000  90.000000  90.000000
      Sites (2)
        #  SP       a     b     c
      ---  ----  ----  ----  ----
        0  Si    0     0     0
        1  Si    0.25  0.25  0.25

    Note, the corresponding lattice vectors are the convention, and do not correspond to
    the atomic basis (for is for the primitive cell)

    Data extracted from cif using:
        cif_to_pymatgen('/Users/alexanderbuccheri/Python/pycharm_projects/wp_benchmarks/silicon/si.cif')
    """
    # Angstrom (from cif)
    lattice = 5.430000 * np.array([[0.000000,  0.500000,  0.500000],
                                   [0.500000,  0.000000,  0.500000],
                                   [0.500000,  0.500000,  0.000000]])

    fractional_positions = np.array([[0.00,   0.00,   0.00],
                                     [0.25,   0.25,   0.25]])

    # In angstrom
    positions = np.transpose(np.transpose(lattice) @ np.transpose(fractional_positions))
    unit='angstrom'
    elements = ['Si', 'Si']
    atomic_numbers = [14, 14]


#cif_to_pymatgen('/Users/alexanderbuccheri/Python/pycharm_projects/wp_benchmarks/silicon/si.cif', True)

import xml.etree.ElementTree as ET
import numpy as np
from typing import Tuple, Optional, Union, List

import ase


# TODO(Alex0 Move to excitingtools
def extract_charge_density(file_name: str) -> Tuple[np.ndarray, np.ndarray]:
    """ Extract charge density from RHO1D.xml file.

    :param file_name: File name
    :return: distance and charge_density
    """
    tree = ET.parse(file_name)
    root = tree.getroot()

    distance = []
    charge_density = []
    for child in root[1][2]:
        distance.append(child.attrib['distance'])
        charge_density.append(eval(child.attrib['value']))

    return np.asarray(distance), np.asarray(charge_density)


def get_standardised_band_path(lattice_vectors) -> Tuple[np.ndarray, dict]:
    """ ASE standardised band path and a fixed k-grid sampling the path.

    Notes:
      Check out band_path.plot()

    :param lattice_vectors: Lattice vectors stored row wise np array or as [a, b, c]
    and (most likely) in angstrom.
    :return: Tuple of the band path and high symmetry points {symbol: k_point/s}
    """
    cell = ase.atoms.Cell(lattice_vectors)
    band_path: ase.dft.kpoints.BandPath = cell.bandpath()
    return band_path.path, band_path.special_points


def exciting_band_path_xml(symbolic_path, high_symmetry_points, steps:Optional[int]=100):
    """ XML-formatted band structure path for exciting.

      <bandstructure>
         <plot1d>
            <path steps="100">
               <point coord="1.0     0.0     0.0" label="Gamma"/>
               <point coord="0.625   0.375   0.0" label="K"/>
               <point coord="0.5     0.5     0.0" label="X"/>
               <point coord="0.0     0.0     0.0" label="Gamma"/>
               <point coord="0.5     0.0     0.0" label="L"/>
            </path>
         </plot1d>
      </bandstructure>

    :param symbolic_path:
    :param high_symmetry_points:
    :param steps:
    :return:
    """
    string = f"""<bandstructure>
    <plot1d>
      <path steps="{int(steps)}">\n"""

    indent = ' ' * 8

    # Identify which high symmetry points are followed by a discontinuity in the band path
    indices = [i-1 for i, value in enumerate(symbolic_path) if value == ","]
    break_point_str = [''] * len(symbolic_path)
    for i in indices:
        break_point_str[i] = 'breakafter="true"'

    # Iterate over high-symmetry points
    for i, symbol in enumerate(symbolic_path):
        if symbol == ',': continue
        point_str = " ".join(str(x) for x in high_symmetry_points[symbol]).strip()
        line = f'<point coord="{point_str}" label="{symbol}" {break_point_str[i]}/>'
        string += indent + line + '\n'

    string += """     </path>
   </plot1d>
</bandstructure>
    """

    return string


# TODO(Alex0 Move to excitingtools
def parse_bandstructure_xml(file_name: str) -> Tuple[np.ndarray, np.ndarray]:
    """ Parse bandstructure.xml

    TODO(Alex) Write a test
              Could extend to parse the vertices, if desired

    :param file_name: bandstructure.xml prepended by path.
    :return: Tuple of discrete set of points sampling the k-path, and band energies
    of .shape = (n_kpts, n_bands).
    """
    tree = ET.parse(file_name)
    root = tree.getroot()

    # Split bands and vertices
    elements = list(root)
    # title = elements[0].text
    bs_xml = {'band': [], 'vertex': []}
    for item in elements[1:]:
        bs_xml[item.tag].append(item)

    n_bands = len(bs_xml['band'])
    first_band = bs_xml['band'][0]
    n_kpts = len(list(first_band))

    # Same set of flattened k-points, per band - so parse once
    k_vector = np.empty(shape=n_kpts)
    for ik, point in enumerate(list(first_band)):
        k_vector[ik] = point.get('distance')

    # Useful to note - then remove
    # print(point.tag, point.items(), point.keys(), point.get('distance'))

    # Read E(k), per band
    band_energies = np.empty(shape=(n_kpts, n_bands))
    for ib, band in enumerate(bs_xml['band']):
        points = list(band)
        for ik, point in enumerate(points):
            band_energies[ik, ib] = point.get('eval')

    return k_vector, band_energies


class BandStructure:
    """ Wrap free functions for processing exciting band structures.
    """
    flat_k_path: Union[list, np.ndarray]
    band_energies:  np.ndarray

    def __init__(self, file_name: Optional[str]):
        self.file_name = file_name

    def parse(self):
        self.flat_k_path, band_energies = parse_bandstructure_xml(self.file_name)

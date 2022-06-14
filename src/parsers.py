import xml.etree.ElementTree as ET
import numpy as np
from typing import Tuple, Optional, List

import ase


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


def find_discontinuities(path: list) -> List[bool]:
    """ Find and return high-symmetry points that end in a discontinuous band path

    For example:
       WLK,UX
    would return the mask [False, False, True, False, False, False],
    indicating that K is the end of a continuous path.

    :return: List of
    """
    mask = [False] * len(path)
    indices = [i-1 for i, value in enumerate(path) if value == ","]
    for i in indices:
        mask[i] = True
    return mask

import numpy as np
import os
import re
from typing import Tuple, Optional, List
import ase

from excitingtools.dataclasses.band_structure import BandData
from excitingtools.dataclasses.data_structs import BandIndices, PointIndex
from excitingtools.dataclasses.eigenvalues import EigenValues
from excitingtools.exciting_obj_parsers.gw_eigenvalues import NitrogenEvalQPColumns, gw_eigenvalue_parser

from conversions import ha_to_ev


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


def exciting_band_path_xml(symbolic_path, high_symmetry_points, steps: Optional[int] = 100):
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
    indices = [i - 1 for i, value in enumerate(symbolic_path) if value == ","]
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
    indices = [i - 1 for i, value in enumerate(path) if value == ","]
    for i in indices:
        mask[i] = True
    return mask


def get_number_of_valence_bands(band_data: BandData, energy_zero: Optional[float] = 0.0) -> int:
    """ For a system with a band gap, find the number of valence bands.

    :param band_data: Object containing band structure data
    :param energy_zero: Zero of the energy, under which all valence band states lie.
    Convention taken to be 0.0 eV.
    :return: Number of valence bands.
    """
    max_energy_per_band = np.empty(shape=band_data.n_bands)
    for ib in range(0, band_data.n_bands):
        max_energy_per_band[ib] = np.amax(band_data.bands[:, ib]) * ha_to_ev
    n_valence = len(max_energy_per_band[max_energy_per_band <= energy_zero])
    return n_valence


def get_gw_bandedge_indices(path: str) -> BandIndices:
    """ Get GW band indices from GW_INFO.OUT, for exciting Nitrogen.

    Note, the strings may differ in Oxygen.

    :param path: Path to file `GW_INFO.OUT`
    :return (vbm, cbm): VBM and CBm indices.
    """
    file_name = os.path.join(path, "GW_INFO.OUT")
    try:
        with open(file=file_name) as fid:
            file_string = fid.read()
    except FileNotFoundError:
        raise FileNotFoundError(f'{file_name} cannot be found')

    # In both cases, take the last match, which corresponds to the GW indice.
    vb_max_str = re.findall(r'\s*Band index of VBM: .*$', file_string, flags=re.MULTILINE)[-1]
    cb_min_str = re.findall(r'\s*Band index of CBM: .*$', file_string, flags=re.MULTILINE)[-1]

    vbm = int(vb_max_str.split()[-1])
    cbm = int(cb_min_str.split()[-1])

    return BandIndices(VBM=vbm, CBm=cbm)


def get_gw_bandedge_k_indices(path: str) -> List[PointIndex]:
    """ Parse k-points at band edges.

    Valid for Nitrogen.

    Fundamental BandGap (eV):                 1.4935
         at k(VBM) =    0.000   0.000   0.000 ik =     1
            k(CBM) =    0.000   0.500   0.500 ik =     3

    :param path:
    :return:
    """
    file_name = os.path.join(path, "GW_INFO.OUT")
    try:
        with open(file=file_name) as fid:
            file_string = fid.read()
    except FileNotFoundError:
        raise FileNotFoundError(f'{file_name} cannot be found')

    # Check if k-points are defined w.r.t. VBM and CBM
    # This implies an indirect gap.
    line = re.findall(r'^.*k\(VBM\) = .*$', file_string, flags=re.MULTILINE)
    direct_gap = line == []

    if direct_gap:
        # Last line match will correspond to GW
        line = re.findall(r'^.*at k\s*= .*$', file_string, flags=re.MULTILINE)[-1]
        ik_vbm = int(line.split()[-1])
        k_vbm = [float(x) for x in line.split()[3:6]]
        ik_cbm = ik_vbm
        k_cbm = k_vbm

    else:
        # Last line match will correspond to GW
        vbm_line = re.findall(r'^.*k\(VBM\) = .*$', file_string, flags=re.MULTILINE)[-1]
        cbm_line = re.findall(r'^.*k\(CBM\) = .*$', file_string, flags=re.MULTILINE)[-1]

        ik_vbm = int(vbm_line.split()[-1])
        k_vbm = [float(x) for x in vbm_line.split()[3:6]]

        ik_cbm = int(cbm_line.split()[-1])
        k_cbm = [float(x) for x in cbm_line.split()[2:5]]

    return [PointIndex(k_vbm, ik_vbm), PointIndex(k_cbm, ik_cbm)]


def return_zero_if_no_file(func):
    """ Decorate `band_gaps` and return zeros if the file cannot be found
    (as the parser will throw an immediate exception).
    """
    def modified_func(path):
        full_file = os.path.join(path, 'EVALQP.DAT')
        if not os.path.isfile(full_file):
            print(f'{full_file} does not exist')
            return {'fundamental': 0.0, 'gamma': 0.0}
        return func(path)
    return modified_func


@return_zero_if_no_file
def band_gaps(path) -> dict:
    """ Parse eigenvalues, return band gaps in Ha.

    Written for exciting Nitrogen.

    :param path: Path to file.
    :return: Dict of band gaps.
    """
    eigenvalues: EigenValues = gw_eigenvalue_parser(path, NitrogenEvalQPColumns.E_GW)

    # Indirect (or smallest)
    band_indices = get_gw_bandedge_indices(path)
    k_valence, k_conduction = get_gw_bandedge_k_indices(path)
    fundamental_gap = eigenvalues.band_gap(band_indices, k_indices=[k_valence.index, k_conduction.index])

    # Direct at Gamma
    ik = eigenvalues.get_index([0., 0., 0.])
    direct_gamma_gap = eigenvalues.band_gap(band_indices, k_indices=[ik, ik])

    return {'fundamental': fundamental_gap, 'gamma': direct_gamma_gap}

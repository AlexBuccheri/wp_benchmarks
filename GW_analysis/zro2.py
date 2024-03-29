""" Post-process ZrO2 G0W0 results.
"""
import os
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.axes import Axes
import numpy as np
from typing import Tuple

from src.conversions import ha_to_ev
from src.parsers import band_gaps


def band_gaps_vs_q(root: str, q_points: list) -> dict:
    """ Get fundamental and direct Gamma gaps as a function of q.

    :param path: Path to file.
    :param q_points: List of q-points.
    :return: Fundamental and Gamma-Gamma gaps.
    """
    fundamental = np.empty(shape=len(q_points))
    gamma_gamma = np.empty(shape=len(q_points))

    for iq, q_point in enumerate(q_points):
        q_str = "".join(str(x) for x in q_point)
        path = os.path.join(root, q_str)
        gaps: dict = band_gaps(path)
        fundamental[iq] = gaps['fundamental']
        gamma_gamma[iq] = gaps['gamma']

    return {'fundamental': fundamental, 'gamma_gamma': gamma_gamma}


def print_data(n_empty, q_points, gaps):
    print(f'Band gap (eV) using {n_empty} % of empty states')

    print('# q-grid, Fundamental Gap, Gamma - Gamma Gap')
    for i in range(0, len(q_points)):
        if np.isclose(gaps['fundamental'][i], 0.):
            continue
        print(q_points[i], gaps['fundamental'][i] * ha_to_ev, gaps['gamma_gamma'][i] * ha_to_ev)


def initialise_bandgap_plot(title: str, y_limits) -> Tuple[Figure, Axes]:
    fig, ax = plt.subplots()
    ax.set_title(title)
    fig.set_size_inches(14, 10)
    plt.rcParams.update({'font.size': 16})
    ax.set_xlabel('q-points', fontsize=16)
    ax.set_ylabel('Quasiparticle Band Gap (eV)', fontsize=16)
    ax.set_ylim(y_limits[0], y_limits[1])
    return fig, ax


if __name__ == "__main__":
    # Notebook Results directory
    root = 'GW_results/zro2/rgkmax8'

    q_points = [[2, 2, 2], [4, 4, 4], [6, 6, 6]]
    total_q_points = [np.prod(q_point) for q_point in q_points]
    n_empty = 100

    # Print data
    gaps = band_gaps_vs_q(root, q_points)
    print_data(n_empty, q_points, gaps)

    # Plot data
    fig, ax = initialise_bandgap_plot('Convergence in Indirect Gap of ZrO2', (5.3, 5.4))
    ax.plot(total_q_points, gaps['fundamental'] * ha_to_ev, 'o', markersize=10, label=str(n_empty))
    plt.show()

    # Must be a smarter way of doing this
    fig, ax = initialise_bandgap_plot('Convergence in Gamma-Gamma Gap of ZrO2', (5.8, 6.0))
    ax.plot(total_q_points, gaps['gamma_gamma'] * ha_to_ev, 'o', markersize=10, label=str(n_empty))
    plt.show()

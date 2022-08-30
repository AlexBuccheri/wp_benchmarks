""" Post-process ZnO G0W0 results.
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
    :return: Fundamental gap
    """
    fundamental = np.empty(shape=len(q_points))

    for iq, q_point in enumerate(q_points):
        q_str = "".join(str(x) for x in q_point)
        path = os.path.join(root, q_str)
        gaps: dict = band_gaps(path)
        fundamental[iq] = gaps['fundamental']

    return {'fundamental': fundamental}


def initialise_gamma_gap_plot() -> Tuple[Figure, Axes]:
    fig, ax = plt.subplots()
    ax.set_title('Convergence in Gamma-Gamma Gap of ZnO')
    fig.set_size_inches(14, 10)
    plt.rcParams.update({'font.size': 16})
    ax.set_xlabel('q-points', fontsize=16)
    ax.set_ylabel('Quasiparticle Band Gap (eV)', fontsize=16)
    ax.set_ylim(1.8, 3.0)
    return fig, ax


def print_data(n_empty, q_points, gaps):
    print(f'Band gap (eV) using {n_empty} % of empty states')
    print('# q-grid, Fundamental Gap (direct)')

    NULL = 0.
    for i in range(0, len(q_points)):
        if np.isclose(gaps['fundamental'][i], NULL):
            continue
        print(q_points[i], gaps['fundamental'][i] * ha_to_ev)


if __name__ == "__main__":
    # Notebook Results directory
    root = 'GW_results/results/zno'
    subdirectory = '_percent_empty'

    q_points = [[2, 2, 1], [4, 4, 3], [6, 6, 4], [8, 8, 5]]
    total_q_points = [np.prod(q_point) for q_point in q_points]
    n_empties = [25, 50, 100]

    # Print data
    for n_empty in n_empties:
        gaps = band_gaps_vs_q(os.path.join(root, str(n_empty) + subdirectory), q_points)
        print_data(n_empty, q_points, gaps)

    # Plot data
    fig, ax = initialise_gamma_gap_plot()
    for n_empty in n_empties:
        gaps = band_gaps_vs_q(os.path.join(root, str(n_empty) + subdirectory), q_points)
        ax.plot(total_q_points, gaps['fundamental'] * ha_to_ev, 'o', markersize=10, label=str(n_empty))

    ax.legend(title='Number of Emtpy States (%)')
    plt.show()

    save_plots = False
    if save_plots:
        plt.savefig('zno_bs.jpeg', dpi=300, facecolor='w', edgecolor='w',
                    orientation='portrait', transparent=True, bbox_inches=None, pad_inches=0.1)


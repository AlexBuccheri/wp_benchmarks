""" Post-process TiO2 G0W0 results.
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

    :param root: Path to file.
    :param q_points: List of q-points.
    :return: Fundamental gap
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


def initialise_bandgap_plot(title: str, y_limits) -> Tuple[Figure, Axes]:
    fig, ax = plt.subplots()
    ax.set_title(title)
    fig.set_size_inches(14, 10)
    plt.rcParams.update({'font.size': 16})
    ax.set_xlabel('q-points', fontsize=16)
    ax.set_ylabel('Quasiparticle Band Gap (eV)', fontsize=16)
    ax.set_ylim(y_limits[0], y_limits[1])
    return fig, ax


def print_data(n_empty, q_points, gaps):
    print(f'Band gap (eV) using {n_empty} % of empty states')
    print('# q-grid, Fundamental Gap, Gamma - Gamma Gap')
    NULL = 0.

    for i in range(0, len(q_points)):
        if np.isclose(gaps['fundamental'][i], NULL):
            continue
        print(q_points[i], gaps['fundamental'][i] * ha_to_ev, gaps['gamma_gamma'][i] * ha_to_ev)


if __name__ == "__main__":
    # Notebook Results directory
    root = 'GW_results/tio2'
    subdirectory = '_percent_empty'

    q_points = [[2, 2, 3], [4, 4, 6], [6, 6, 10]]
    total_q_points = [np.prod(q_point) for q_point in q_points]
    n_empties = [10, 25, 50]   # 100

    # Print data
    for n_empty in n_empties:
        path = os.path.join(root, str(n_empty) + subdirectory)
        gaps = band_gaps_vs_q(path, q_points)
        print_data(n_empty, q_points, gaps)

    # Fewer frequency points
    print('Reduced frrequency points from 32 to 16'
          'Error looks large')
    path = os.path.join(root, '10_percent_empty_freq16')
    gaps = band_gaps_vs_q(path, q_points)
    print_data('10', q_points, gaps)

    # Plot data
    print("As the q-grid becomes more dense, TiO2 switches from direct gap to indirect gap")

    fig, ax = initialise_bandgap_plot('Convergence in Fundamental Gap of TiO2', (3.25, 3.5))
    for n_empty in n_empties:
        gaps = band_gaps_vs_q(os.path.join(root, str(n_empty) + subdirectory), q_points)
        ax.plot(total_q_points, gaps['fundamental'] * ha_to_ev, 'o', markersize=10, label=str(n_empty))

    ax.legend(title='Number of Emtpy States (%)')
    plt.show()

    fig, ax = initialise_bandgap_plot('Convergence in Gamma-Gamma Gap of TiO2', (3.25, 3.50))
    for n_empty in n_empties:
        gaps = band_gaps_vs_q(os.path.join(root, str(n_empty) + subdirectory), q_points)
        ax.plot(total_q_points, gaps['gamma_gamma'] * ha_to_ev, 'o', markersize=10, label=str(n_empty))

    ax.legend(title='Number of Emtpy States (%)')
    plt.show()

    save_plots = False
    if save_plots:
        plt.savefig('tio2_bs.jpeg', dpi=300, facecolor='w', edgecolor='w',
                    orientation='portrait', transparent=True, bbox_inches=None, pad_inches=0.1)

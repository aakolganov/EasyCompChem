"""
A script to plot reaction energy profiles with distinct colors and approximate/missing values.
"""

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
import os


def parse_energy_value(value):
    """
    Parse energy value from string or numeric input.
    Returns a tuple (energy: float or None, is_approx: bool)
    """
    if isinstance(value, str):
        if value.strip() == '???':
            return (None, False)
        elif value.strip().startswith('~'):
            # Approximate value
            approx_val = value.strip()[1:]
            if approx_val == '???':
                # Approximate missing data (~???), treat as missing
                return (None, True)
            try:
                val = float(approx_val)
                return (val, True)
            except:
                raise ValueError(f"Cannot parse approximate value: {value}")
        else:
            # Try to parse directly as float
            try:
                val = float(value)
                return (val, False)
            except:
                raise ValueError(f"Cannot parse energy value: {value}")
    else:
        # Numeric value
        return (float(value), False)

def plot_reaction_energy_profile(energy_lists, compound_labels,
                                          series_labels=None,
                                          title="Reaction Energy Profile",
                                          xlabel="Reaction Coordinate",
                                          ylabel="Energy (kJ/mol)",
                                          figsize=(14, 8),
                                          colors=None,
                                          markers=None,
                                          save_path=None,
                                          dpi=300,
                                          show_plot=True):
    """
    Plot reaction energy profiles with distinct colors and approximate/missing values.
    """

    def parse_energy_value(value):
        if isinstance(value, str):
            if value.strip() == '???':
                return (None, False)
            elif value.strip().startswith('~'):
                approx_val = value.strip()[1:]
                if approx_val == '???':
                    return (None, True)
                try:
                    val = float(approx_val)
                    return (val, True)
                except:
                    raise ValueError(f"Cannot parse approximate value: {value}")
            else:
                try:
                    val = float(value)
                    return (val, False)
                except:
                    raise ValueError(f"Cannot parse energy value: {value}")
        else:
            return (float(value), False)

    if colors is None:
        colors = [
            '#1f77b4',  # Strong blue
            '#ff7f0e',  # Bright orange
            '#2ca02c',  # Green
            '#d62728',  # Red
            '#9467bd',  # Purple
            '#8c564b',  # Brown
            '#e377c2',  # Pink
            '#7f7f7f',  # Gray
            '#bcbd22',  # Olive
            '#17becf',  # Cyan
            '#ff1493',  # Deep pink
            '#32cd32'   # Lime green
        ]

    fig, ax = plt.subplots(figsize=figsize)
    x_coords = np.arange(len(compound_labels))

    if series_labels is None:
        series_labels = [f'Series {i+1}' for i in range(len(energy_lists))]

    if markers is None:
        markers = ['o', 's', '^', 'D', 'v', '<', '>', 'p', '*', 'h'] * 10

    for i, energies_raw in enumerate(energy_lists):
        energies = []
        approx_flags = []
        missing_flags = []

        for val in energies_raw:
            energy, approx = parse_energy_value(val)
            energies.append(energy)
            approx_flags.append(approx)
            missing_flags.append(energy is None)

        energies_corrected = energies.copy()
        for idx, e in enumerate(energies):
            if e is None:
                neighbors = []
                for left_idx in range(idx-1, -1, -1):
                    if energies[left_idx] is not None:
                        neighbors.append(energies[left_idx])
                        break
                for right_idx in range(idx+1, len(energies)):
                    if energies[right_idx] is not None:
                        neighbors.append(energies[right_idx])
                        break
                if neighbors:
                    fill_value = max(neighbors) + 20
                else:
                    fill_value = 20
                energies_corrected[idx] = fill_value

        # Use distinct color for each series
        series_color = colors[i % len(colors)]

        # Plot with thicker lines for better visibility
        ax.plot(x_coords, energies_corrected, color=series_color,
                linestyle='-', linewidth=3, alpha=0.9, label=series_labels[i])

        for j, energy_val in enumerate(energies_corrected):
            if missing_flags[j]:
                ax.plot(x_coords[j], energy_val, marker='X', markersize=12,
                        color='red', markerfacecolor='red',
                        markeredgecolor='darkred', markeredgewidth=2)
            elif approx_flags[j]:
                ax.plot(x_coords[j], energy_val, marker=markers[i % len(markers)],
                        markersize=11, color=series_color, markerfacecolor='white',
                        markeredgewidth=3, markeredgecolor=series_color)
            else:
                ax.plot(x_coords[j], energy_val, marker=markers[i % len(markers)],
                        markersize=9, color=series_color, markerfacecolor=series_color,
                        markeredgecolor='white', markeredgewidth=1)

        for j, val in enumerate(energies_raw):
            if str(val).strip() in ['???', '~???']:
                label_val = '???'
                label_color = 'red'
                font_weight = 'bold'
            elif str(val).startswith('~'):
                if str(val) == '~???':
                    label_val = '???'
                    label_color = 'red'
                    font_weight = 'bold'
                else:
                    label_val = f'~{energies_corrected[j]:.1f}'
                    label_color = '#FF6600'  # Distinct orange for approx
                    font_weight = 'bold'
            else:
                label_val = f'{energies_corrected[j]:.1f}'
                label_color = 'black'
                font_weight = 'normal'

            ax.annotate(label_val, (x_coords[j], energies_corrected[j]),
                        textcoords='offset points',
                        xytext=(0, 15), ha='center', fontsize=10,
                        color=label_color, fontweight=font_weight,
                        bbox=dict(boxstyle="round,pad=0.3",
                                facecolor='white',
                                edgecolor=label_color,
                                alpha=0.95))

    ax.set_xticks(x_coords)
    ax.set_xticklabels(compound_labels, rotation=45, ha='right', fontsize=11)
    ax.set_xlabel(xlabel, fontsize=13, fontweight='bold')
    ax.set_ylabel(ylabel, fontsize=13, fontweight='bold')
    ax.set_title(title, fontsize=15, fontweight='bold', pad=25)
    ax.grid(True, alpha=0.3, linestyle='--')

    # Main legend with better styling
    main_legend = ax.legend(loc='center left', bbox_to_anchor=(1.02, 0.7),
                            frameon=True, fancybox=True, shadow=True,
                            fontsize=11, title='Pathway', title_fontsize=12,
                            borderpad=1, columnspacing=1, handletextpad=0.8)
    ax.add_artist(main_legend)

    legend_elements = [
        Line2D([0], [0], marker='o', color='gray', linestyle='None',
               markerfacecolor='gray', markersize=8, label='Completed'),
        Line2D([0], [0], marker='o', color='gray', linestyle='None',
               markerfacecolor='white', markeredgewidth=2, markersize=8, label='N (~)'),
        Line2D([0], [0], marker='X', color='red', linestyle='None',
               markersize=8, label='Missing (???)')
    ]

    legend2 = ax.legend(handles=legend_elements, loc='center left', bbox_to_anchor=(1.02, 0.25),
                       frameon=True, fancybox=True, shadow=True,
                       fontsize=10, title='Calculation Status', title_fontsize=11,
                       borderpad=1, columnspacing=1, handletextpad=0.8)

    plt.tight_layout()
    plt.subplots_adjust(right=0.78)  # More space for legends

    saved_files = []
    if save_path:
        directory = os.path.dirname(save_path) if os.path.dirname(save_path) else '.'
        os.makedirs(directory, exist_ok=True)

        png_file = f"{save_path}.png"
        plt.savefig(png_file, dpi=dpi, bbox_inches='tight',
                   facecolor='white', edgecolor='none', format='png')
        saved_files.append(png_file)

        svg_file = f"{save_path}.svg"
        plt.savefig(svg_file, bbox_inches='tight',
                   facecolor='white', edgecolor='none', format='svg')
        saved_files.append(svg_file)

        print(f"Plot saved successfully!")
        print(f"PNG: {png_file} ({dpi} DPI)")
        print(f"SVG: {svg_file} (vector format for Adobe Illustrator)")

    if show_plot:
        plt.show()

    if save_path:
        return fig, ax, saved_files
    else:
        return fig, ax

#example usage

if __name__ == "__main__":
    # Example 1: Mixed data types
    compound_labels = ['RuMACHO', 'RuMACHO+CO2', 'RuMACHO-HCOO',
                       'RuMACHO-HCOO+H2', 'TS(H2 insertion)', 'RuMACHO-H2+HCOO',
                       'RuMACHO-H2+HCOO+amine', 'TS2', 'RuMACHO+HCOOH+amine',
                       '2HCOOH+amine', 'TS3','HCOOH+tetrahedr. int.', 'TS4',
                       'Formamide+HCOOH']

    # Series with exact, approximate, and missing values
    energy_series_1 = [0, -6, -57,
                       -66, '~-1', -52,
                       -156, -103, -142,
                       -169, -113, -120, 45,
                       -83]



    energy_series_2 = [0, -26, -5,
                       -66, '~-1', -52,
                       -15, -104    , -16,
                       -16, -10, -19, '~44',
                       -10]      # Mix of types      # Approximate missing

    fig = plot_reaction_energy_profile(
        [energy_series_1, energy_series_2],
        compound_labels,
        series_labels=['e2e', 'e2'],
        title='CO2 hydrogenation Energy Profile',
        ylabel='Relative Electronic Energy (kJ/mol)',
        save_path='figures/reaction_profile'
    )

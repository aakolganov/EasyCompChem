"""
A script to plot interactive linear regression plots with hover features
"""

import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import linregress
import plotly.express as px
import plotly.io as pio
import numpy

def plot_linear_correlation(
    data,
    x_feature,
    y_feature,
    hover_feature=None,
    graph_name=None,
    html_name=None,
    x_label=None,       # Optional matplotlib x-axis label
    y_label=None        # Optional matplotlib y-axis label
):
    # Get data
    columns = [x_feature, y_feature]
    if hover_feature:
        columns.append(hover_feature)

    df_clean = data[columns].dropna()
    x = df_clean[x_feature]
    y = df_clean[y_feature]

    # Linear Regression
    slope, intercept, r_value, p_value, std_err = linregress(x, y)

    # === Matplotlib Plot ===
    plt.figure(figsize=(6, 4.5))
    plt.scatter(x, y, s=40, color='#1f77b4', alpha=0.7, edgecolors='k', linewidths=0.5)
    equation = rf'$y = {slope:.3f}x {"+" if intercept >= 0 else "-"} {abs(intercept):.3f}$'
    r_squared = rf'$R^2 = {r_value**2:.3f}$'
    label = equation + '\n' + r_squared

    plt.plot(x, intercept + slope * x, 'r-', linewidth=2, label=label)

    ax = plt.gca()  # Get current axes
    ax.grid(False)  # No grid, but keep ticks and axes

        # Keep ticks and axes
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    ax.tick_params(bottom=True, left=True, direction='out')
    # Make sure the axis lines (spines) are visible
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color('black')       # Optional: match your style
        spine.set_linewidth(1.0)       # Optional: control thickness


    # Set grey background for the plot area
    ax.set_facecolor('#ffffff')  # Light grey — change if needed

    # Optional: also change the outer (figure) background if desired
    plt.gcf().set_facecolor('#ffffff')




    # Use custom axis labels if provided
    plt.xlabel(x_label if x_label else x_feature, fontsize=14)
    plt.ylabel(y_label if y_label else y_feature, fontsize=14)

    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.legend(frameon=True, fontsize=12)
    # Turn off all grid lines explicitly
    #ax.grid(False)                      # Disable major grid
    #ax.grid(visible=False, which='minor')  # Disable minor grid
    plt.tight_layout()

    if graph_name:
        plt.savefig(graph_name, dpi=300, bbox_inches='tight')
    plt.show()

    # === Plotly Plot ===
    fig = px.scatter(
        data,
        x=x_feature,
        y=y_feature,
        hover_name=hover_feature,
        title=f'Interactive correlation between {x_feature} and {y_feature}',
        color=hover_feature,
        template='plotly_white'
    )

    fig.update_traces(marker=dict(size=15, line=dict(width=1, color='DarkSlateGrey')))

    fig.update_layout(
        font=dict(size=14),
        title_font=dict(size=18),
        xaxis_title=x_feature,
        yaxis_title=y_feature,
        legend_title=hover_feature if hover_feature else "Legend",
        margin=dict(l=60, r=60, t=60, b=60),
        hoverlabel=dict(bgcolor="white", font_size=12)
    )

    fig.show()

    if html_name:
        fig.write_html(html_name)

#example
if __name__ == '__main__':
    data = pd.read_csv(
        'cleaned_data_man.csv',
        encoding='latin1'  # or try 'cp1252'
    )


    plot_linear_correlation(
        data,
        x_feature="dPARA_z",
        y_feature = 'Formation Electronic Energy',
        hover_feature="Molecule",
        graph_name="Graphs/cleaned_df/dPARA_z_vs_en.png",
        html_name="Graphs/cleaned_df/dPARA_z_vs_en.html",
        x_label=r'$\Delta\delta^{\mathrm{para}}_{33}$ ,ppm',  # For matplotlib
        y_label=r'$\Delta E_{\mathrm{formation}}$ (kJ mol$^{-1}$)'
    )

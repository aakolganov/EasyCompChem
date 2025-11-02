"""
This module processes a dataset by applying various transformations to numeric
columns and computes pairwise Pearson correlation coefficients and R² values
between original and transformed data columns. Additionally, it visualizes the
results through heatmaps.

The primary goal is to analyze the relationship between variable transformations
and their mutual dependencies. The module uses correlation matrices, R² matrices,
and visual heatmaps to illustrate these relationships.

Dependencies:
- pandas
- seaborn
- matplotlib.pyplot

The output includes correlation matrices, R² matrices, and saved heatmaps
for each transformation applied.
"""

import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt

data_dropped = pd.read_csv('data.csv')



# Select only numeric columns
#numeric_df = data_dropped.select_dtypes(include=['number'])

# Define transformation functions.
# Note: For logarithmic, we use log(|x|+1) to avoid issues with negatives or zero.
transformations = {
    '': lambda x: x,           # Identity transformation: regular linear correlation
    #'quadratic': lambda x: x**2,
    #'cubic': lambda x: x**3,
    #'inverse_square': lambda x: np.where(x == 0, np.nan, 1/(x**2)),
    #'inverse_cube': lambda x: np.where(x == 0, np.nan, 1/(x**3)),
    #'logarithmic': lambda x: np.log(np.abs(x) + 1)
}
# Loop over each transformation.
for name, func in transformations.items():
    # Apply the transformation to every numeric column
    transformed_df = data_dropped.apply(func)

    # Create an output folder for scatter plots for this transformation.
    folder_name = f"scatter_plots/scatter_plots_{name}"
    #os.makedirs(folder_name, exist_ok=True)

    # Create an empty DataFrame to store the one-sided correlations.
    # Here, for each pair (i, j), we correlate the original column i with the transformed version of column j.
    corr_matrix = pd.DataFrame(index=data_dropped.columns, columns=data_dropped.columns, dtype=float)

    # Loop over every pair of original and transformed columns.
    for orig_col in data_dropped.columns:
        for trans_col in data_dropped.columns:
            # Compute Pearson correlation between original column and transformed column.
            corr_val = transformed_df[orig_col].corr(transformed_df[trans_col])
            corr_matrix.loc[orig_col, trans_col] = corr_val

            # # Generate scatter plot for this pair.
            # plt.figure(figsize=(8, 6))
            # plt.scatter(numeric_df[orig_col], transformed_df[trans_col], alpha=0.6)
            # plt.xlabel(orig_col)
            # plt.ylabel(f"{trans_col} ({name} transformed)")
            # plt.title(f"{orig_col} vs. {trans_col}\nPearson r = {corr_val:.2f}")
            #
            # # Save the scatter plot in the corresponding folder.
            # plot_filename = os.path.join(folder_name, f"{orig_col}_vs_{trans_col}.png")
            # plt.savefig(plot_filename)
            # plt.close()

    # Compute the R² matrix as the square of the correlation coefficients.
    r2_matrix = corr_matrix ** 2

    # Print the matrices to the console.
    print(f"\nCorrelation matrix {name}")
    print(corr_matrix)
    print(f"\nR² matrix {name}")
    print(r2_matrix)

    # Create and save the R² heatmap for this transformation.
    plt.figure(figsize=(18, 18))
    sns.heatmap(r2_matrix, annot=True, cmap='coolwarm', fmt=".2f")
    plt.title(f'Pairwise R² Matrix Heatmap - {name.capitalize()} Transformation\n(Original vs. Transformed)')

    heatmap_filename = f"R2_heatmap_{name}.png"
    plt.savefig(heatmap_filename)
    plt.show()  # Display the heatmap plot
    plt.close()
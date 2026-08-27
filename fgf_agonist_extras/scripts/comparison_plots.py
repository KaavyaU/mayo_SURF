import matplotlib.pyplot as plt
import pandas


# Obtain the original files
original_files = []
runs = 15
for x in range(1,runs + 1):
    original_files.append(f"kaavya_data/kaavya_data/runs/20260622_run{x}/preds/preds_evaluation.csv")

# Obtain filtered file
filtered_file = "filtered_preds_evaluation.csv"


# Convert original csv files to a dataframe
original_frames = []
for file in original_files:
    df = pandas.read_csv(file)
    original_frames.append(df)
original_df = pandas.concat(original_frames, ignore_index = True)

# Convert filtered csv file to a dataframe
filtered_df = pandas.read_csv(filtered_file)

# Create histograms
metrics = ["binder_length", "dG", "dSASA", "dG_dSASA_ratio", "shape_complementarity", "interface_hydrophobicity",
           "binder_helix_pct", "binder_sheet_pct", "binder_loop_pct", "surface_hydrophobicity"]
for metric in metrics:
    original_vals = original_df[metric].to_numpy()
    filtered_vals = filtered_df[metric].to_numpy()

    bins = 20
    plt.hist(original_vals, bins = bins, color = "tab:orange", label = "Predictions")
    plt.hist(filtered_vals, bins = bins, color = "tab:blue", label = "Filtered Predictions")                        

    plt.legend()
    plt.xlabel(metric)
    plt.ylabel('count')

    # Save histogram
    plt.savefig(f"comparison_histograms/{metric}_hist.png")

    # Clear histogram data
    plt.clf()

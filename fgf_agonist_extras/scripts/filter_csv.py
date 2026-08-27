import pandas

# Store csv files for all runs --> Update to match file structure
files = []
runs = 15
for x in range(1,runs + 1):
    files.append(f"kaavya_data/kaavya_data/runs/20260622_run{x}/preds/preds_evaluation.csv")

# Filter the csv files
filtered_frames = []
for file in files:
    df = pandas.read_csv(file)
    # Drop unnecessary columns
    df = df.drop(columns = ["pred_name", "design_name", "design_pdb"])
    # Drop rows with empty sequence (failed test)
    df = df.dropna(subset = ["sequence"])
    # Filter for best sequences
    df = df[df["shape_complementarity"] >= 0.65]
    df = df[df["n_hbonds"] >= 2]
    df = df[df["n_unsat_hbonds"] <= 1]
    df = df[df["surface_hydrophobicity"] <= 0.35]
    filtered_frames.append(df)

final_df = pandas.concat(filtered_frames, ignore_index = True)
final_df.to_csv("filtered_preds_evaluation.csv", index = False)

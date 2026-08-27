import pandas

# Obtain filtered file
filtered_file = "filtered_preds_evaluation.csv"

# Convert filtered csv file to a dataframe
df = pandas.read_csv(filtered_file)

# Get summary stats
summary = df.select_dtypes(include = 'number').agg(["mean", "median", "std", "max", "min"])

# Save to csv file
summary.to_csv("filtered_preds_stats.csv", index = True)
import pandas as pd

in_path  = "params_025.tsv"
out_path = "params_025_seed2to5.tsv"   # change name as you like

# Read (tab-separated). If your file has a commented header or weird spacing, tell me and I’ll adjust.
df = pd.read_csv(in_path, sep=r"\t", engine="python")

# Drop first 30 rows (seed 1 block), keep the rest
df2 = df.iloc[30:].copy()

# Shift idx down by 30 (assumes idx is numeric)
df2["#idx"] = df2["#idx"] - 30

# (Optional) also reset the dataframe row index (not the 'idx' column)
df2 = df2.reset_index(drop=True)

# Save back as TSV
df2.to_csv(out_path, sep="\t", index=False)

print(f"Wrote {len(df2)} rows to {out_path}")
print(df2.head())

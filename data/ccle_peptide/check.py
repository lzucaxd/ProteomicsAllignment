import pandas as pd

# -----------------------------
# Load CCLE normalized protein table
# -----------------------------
file = "Table_S2_Protein_Quant_Normalized.xlsx"
sheet = "Normalized Protein Expression"

print("Loading data...")
df = pd.read_excel(file, sheet_name=sheet)

# -----------------------------
# Identify column types
# -----------------------------
meta_cols = {
    "Protein_Id",
    "Gene_Symbol",
    "Description",
    "Group_ID",
    "Uniprot",
    "Uniprot_Acc"
}

peptide_cols = [c for c in df.columns if str(c).endswith("_Peptides")]

sample_cols = [
    c for c in df.columns
    if c not in meta_cols and c not in peptide_cols
]

data = df[sample_cols]

print("\nDetected columns")
print("Total columns:", len(df.columns))
print("Metadata columns:", len(meta_cols))
print("Peptide-count columns:", len(peptide_cols))
print("Sample columns:", len(sample_cols))

# -----------------------------
# Matrix shape
# -----------------------------
print("\nMatrix shape (proteins × cell lines):", data.shape)

# -----------------------------
# Overall missingness
# -----------------------------
total_values = data.size
missing_values = data.isna().sum().sum()

overall_missing = missing_values / total_values * 100

print("\nOverall Missingness")
print("Total values:", total_values)
print("Missing values:", missing_values)
print("Missing %:", round(overall_missing, 2))

# -----------------------------
# Missingness per protein (gene)
# -----------------------------
missing_per_protein = data.isna().mean(axis=1)

print("\nPer-protein missingness")
print("Median:", round(missing_per_protein.median()*100,2), "%")
print("Mean:", round(missing_per_protein.mean()*100,2), "%")
print("Min:", round(missing_per_protein.min()*100,2), "%")
print("Max:", round(missing_per_protein.max()*100,2), "%")

# -----------------------------
# Missingness per cell line
# -----------------------------
missing_per_cell_line = data.isna().mean(axis=0)

print("\nPer-cell line missingness")
print("Median:", round(missing_per_cell_line.median()*100,2), "%")
print("Mean:", round(missing_per_cell_line.mean()*100,2), "%")
print("Min:", round(missing_per_cell_line.min()*100,2), "%")
print("Max:", round(missing_per_cell_line.max()*100,2), "%")

# -----------------------------
# Save detailed outputs
# -----------------------------
protein_missing_table = pd.DataFrame({
    "GeneSymbol": df["Gene_Symbol"],
    "Protein_Id": df["Protein_Id"],
    "missing_fraction": missing_per_protein
})

cell_line_missing_table = pd.DataFrame({
    "cell_line": missing_per_cell_line.index,
    "missing_fraction": missing_per_cell_line.values
})

protein_missing_table.to_csv("ccle_missing_per_gene.csv", index=False)
cell_line_missing_table.to_csv("ccle_missing_per_cell_line.csv", index=False)

print("\nSaved:")
print("ccle_missing_per_gene.csv")
print("ccle_missing_per_cell_line.csv")

# -----------------------------
# Optional: quick histograms
# -----------------------------
try:
    import matplotlib.pyplot as plt

    plt.figure()
    missing_per_protein.hist(bins=50)
    plt.title("Missingness per protein")
    plt.xlabel("Fraction missing")
    plt.ylabel("Count")
    plt.savefig("missing_per_protein_hist.png")
    plt.close()

    plt.figure()
    missing_per_cell_line.hist(bins=50)
    plt.title("Missingness per cell line")
    plt.xlabel("Fraction missing")
    plt.ylabel("Count")
    plt.savefig("missing_per_cell_line_hist.png")
    plt.close()

    print("Saved histogram plots")

except:
    print("Matplotlib not installed — skipping plots")
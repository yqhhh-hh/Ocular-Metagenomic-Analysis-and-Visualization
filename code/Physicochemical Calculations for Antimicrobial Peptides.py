from Bio.SeqUtils.ProtParam import ProteinAnalysis
import pandas as pd

# === 1. Reading files ===
file_path = "/Users/yqhhh/Desktop/Book3.csv"
df = pd.read_csv(file_path, encoding="utf-8")

# === 2. sequence column name ===
print(df.columns)
seq_col = "Sequence"

# === 3. Cleaning sequence ===
df = df.dropna(subset=[seq_col])

def clean_seq(seq):
    seq = str(seq).upper()
    seq = seq.replace("*", "")
    seq = seq.replace("X", "")
    seq = seq.replace(" ", "")
    seq = seq.replace("\t", "")
    seq = seq.replace("\n", "")
    seq = seq.strip()
    return seq

df[seq_col] = df[seq_col].apply(clean_seq)

# Remove empty sequences after cleaning
df = df[df[seq_col] != ""]

# === 4. Calculation features ===
def analyze_peptide(seq):
    analysed_seq = ProteinAnalysis(seq)

    aa_comp = analysed_seq.amino_acids_percent

    result = {
        "length": len(seq),
        "charge": analysed_seq.charge_at_pH(7.0),
        "hydrophobicity": analysed_seq.gravy()
    }

    for aa, value in aa_comp.items():
        result[f"AA_{aa}"] = value

    return result

features = df[seq_col].apply(analyze_peptide)
features_df = pd.DataFrame(features.tolist())

# === 5. Merge results ===
final_df = pd.concat([df.reset_index(drop=True), features_df], axis=1)

# === 6. Save ===
out_path = "/Users/yqhhh/Desktop/0.5_0.6amp_features.csv"
final_df.to_csv(out_path, index=False)

print("Analysis complete, output file:", out_path)
print(final_df.head())
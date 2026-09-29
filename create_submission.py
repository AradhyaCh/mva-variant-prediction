import pandas as pd

# Required Track 1 CSV schema
columns = [
    "proband_id", "chrom_1", "pos_1", "ref_1", "alt_1",
    "chrom_2", "pos_2", "ref_2", "alt_2", "epcr",
    "finding_type", "notes",
]

proband_id = "PROBAND01"

candidates = [
    {
        "proband_id": proband_id,
        "chrom_1": "15",
        "pos_1": 40209701,
        "ref_1": "T",
        "alt_1": "G",
        "chrom_2": "15",
        "pos_2": 40220612,
        "ref_2": "T",
        "alt_2": "G",
        "epcr": 0.95,
        "finding_type": "primary",
        "notes": "Compound heterozygous BUB1B variants (stop_gained + missense) causing MVA Type 1",
    },
    {
        "proband_id": proband_id,
        "chrom_1": "11",
        "pos_1": 96092210,
        "ref_1": "TTGCTGCTGC",
        "alt_1": "T",
        "chrom_2": "",
        "pos_2": "",
        "ref_2": "",
        "alt_2": "",
        "epcr": 0.85,
        "finding_type": "primary",
        "notes": "Homozygous inframe deletion in CEP57 causing MVA Type 2",
    }
]

df_sub = pd.DataFrame(candidates, columns=columns)

# Write directly to the required Hugging Face filename format
submission_path = r"C:\Users\admin vit\Desktop\Dataset\Yobrrrr_vep_filter.csv"
df_sub.to_csv(submission_path, index=False)

print("=== SUBMISSION PREVIEW ===")
print(df_sub.to_string(index=False))
print(f"\nSaved fixed submission file to: {submission_path}")

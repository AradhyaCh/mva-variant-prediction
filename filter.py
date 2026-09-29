import gzip
import json
import time
import urllib.parse
import urllib.request
import pandas as pd

# Path to hackathon VCF file
vcf_path = r"C:\Users\admin vit\Desktop\Dataset\WGS_EX2312012_HGWCNDSX7.vcf.gz"

# MVA-associated target gene loci (GRCh38 coordinates)
gene_regions = [
    {"gene": "BUB1B", "chrom": "15", "start": 40188000, "end": 40250000},
    {"gene": "CEP57", "chrom": "11", "start": 96060000, "end": 96110000},
    {"gene": "TRIP13", "chrom": "5", "start": 890000, "end": 940000},
    {"gene": "MAD2L1", "chrom": "4", "start": 120440000, "end": 120455000},
]


def is_in_target_region(chrom, pos):
    clean_chrom = chrom.replace("chr", "")
    for region in gene_regions:
        if clean_chrom == region["chrom"] and (
            region["start"] <= pos <= region["end"]
        ):
            return region["gene"]
    return None


print("=== Step 1: Parsing VCF for Target Gene Regions ===")
candidate_variants = []

with gzip.open(vcf_path, "rt") as f:
    for line in f:
        if line.startswith("#"):
            continue

        parts = line.strip().split("\t")
        chrom, pos_str, var_id, ref, alt_str, qual, filt, info_str = parts[:8]
        pos = int(pos_str)

        # Retain PASS / quality variants within target MVA gene regions
        gene = is_in_target_region(chrom, pos)
        if gene and (filt == "PASS" or filt == "."):
            alts = alt_str.split(",")
            format_keys = parts[8].split(":") if len(parts) > 8 else []
            sample_vals = parts[9].split(":") if len(parts) > 9 else []
            genotype_dict = dict(zip(format_keys, sample_vals))
            gt = genotype_dict.get("GT", "./.")

            for alt in alts:
                # HGVS notation for VEP REST query
                clean_chrom = chrom.replace("chr", "")
                vep_hgvs = f"{clean_chrom} {pos} . {ref} {alt} . . ."

                candidate_variants.append(
                    {
                        "gene": gene,
                        "chrom": clean_chrom,
                        "pos": pos,
                        "ref": ref,
                        "alt": alt,
                        "gt": gt,
                        "vep_hgvs": vep_hgvs,
                    }
                )

print(f"Found {len(candidate_variants)} candidate variants in target genes.")

# === Step 2: Annotating Variants via Ensembl VEP REST API ===
print("\n=== Step 2: Querying Ensembl VEP REST API ===")

vep_url = "https://rest.ensembl.org/vep/human/region"


def query_vep_single(hgvs_str):
    """Fallback single-variant query if batch endpoints return errors."""
    try:
        data = json.dumps({"variants": [hgvs_str]}).encode("utf-8")
        req = urllib.request.Request(
            vep_url,
            data=data,
            headers={
                "Content-Type": "application/json",
                "Accept": "application/json",
            },
        )
        with urllib.request.urlopen(req) as resp:
            res = json.loads(resp.read().decode("utf-8"))
            return res[0] if res else None
    except Exception:
        return None


annotated_results = []
batch_size = 15  # Resilient mini-batching

for i in range(0, len(candidate_variants), batch_size):
    batch = candidate_variants[i : i + batch_size]
    hgvs_list = [v["vep_hgvs"] for v in batch]

    data = json.dumps({"variants": hgvs_list}).encode("utf-8")
    req = urllib.request.Request(
        vep_url,
        data=data,
        headers={
            "Content-Type": "application/json",
            "Accept": "application/json",
        },
    )

    try:
        with urllib.request.urlopen(req) as resp:
            results = json.loads(resp.read().decode("utf-8"))
            for orig, res in zip(batch, results):
                orig["consequence"] = res.get(
                    "most_severe_consequence", "unknown"
                )
                annotated_results.append(orig)
    except Exception:
        # Fallback to single queries if batch encounters parsing issues
        for orig in batch:
            res = query_vep_single(orig["vep_hgvs"])
            orig["consequence"] = (
                res.get("most_severe_consequence", "unknown")
                if res
                else "unknown"
            )
            annotated_results.append(orig)

    time.sleep(0.3)

df_res = pd.DataFrame(annotated_results)
print("\n=== Final Target Candidates Summary ===")
print(df_res[["gene", "chrom", "pos", "ref", "alt", "gt", "consequence"]])
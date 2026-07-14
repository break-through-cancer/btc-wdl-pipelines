from cyvcf2 import VCF
import pandas as pd

def get_af(v):
    af = None

    # Try FORMAT/AF first
    try:
        val = v.format("AF")
        if val is not None:
            af = float(val[0][0])
    except Exception:
        pass

    # Try INFO/AF
    if af is None:
        try:
            val = v.INFO.get("AF")
            if isinstance(val, tuple):
                af = float(val[0])
            elif val is not None:
                af = float(val)
        except Exception:
            pass

    # Try FORMAT/AD if AF unavailable
    if af is None:
        try:
            ad = v.format("AD")
            if ad is not None:
                ref_depth = ad[0][0]
                alt_depth = ad[0][1]
                total = ref_depth + alt_depth
                if total > 0:
                    af = alt_depth / total
        except Exception:
            pass

    return af


def load_shared_vcf(path, label):
    rows = []

    for v in VCF(path):
        key = f"{v.CHROM}:{v.POS}:{v.REF}:{','.join(v.ALT)}"
        rows.append({
            "variant_id": key,
            f"AF_{label}": get_af(v),
            f"QUAL_{label}": v.QUAL,
            f"FILTER_{label}": v.FILTER or "PASS",
        })

    return pd.DataFrame(rows)


df1 = load_shared_vcf("isec_results/0002.vcf", "file1")
df2 = load_shared_vcf("isec_results/0003.vcf", "file2")

merged = df1.merge(df2, on="variant_id")

merged["AF_diff"] = merged["AF_file1"] - merged["AF_file2"]
merged["abs_AF_diff"] = merged["AF_diff"].abs()

print(merged.describe())

merged.to_csv("shared_af_concordance.tsv", sep="\t", index=False)
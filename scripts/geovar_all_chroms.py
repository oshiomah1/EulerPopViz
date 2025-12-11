#!/usr/bin/env python3
import os, re, glob, gzip, argparse
import numpy as np
import pandas as pd
from tqdm import tqdm
from geovar import read_pop_panel, GeoVar, vcf_to_freq_table
#2  4 5  8  10 11 
# ================== CONFIG ==================
DATA_PATH = "/quobyte/bmhenngrp/from-lssc0"
VCF_DIR = f"{DATA_PATH}/data/genomes/CAAPA_freeze2_PHASED_common_rare/snps"
POP_PANEL = f"{DATA_PATH}/projects/CAAPA2_functional_annotation/geovar/data/cleaned_population_panel.txt"
RESULTS_DIR = f"{DATA_PATH}/projects/CAAPA2_functional_annotation/geovar/results/freq_mat_files"
os.makedirs(RESULTS_DIR, exist_ok=True)
BINS = [(0,0), (0,0.01), (0.01,0.05), (0.05,0.1), (0.1,1.0)]
# ============================================


def parse_args():
    p = argparse.ArgumentParser(description="Generate GeoVar codes for one or all chromosomes.")
    p.add_argument("--chr", dest="one_chr", help="Process only this chromosome (e.g. 1..22, X, Y)")
    return p.parse_args()

def find_chr_vcfs(vcf_dir, chrs):
    """
    Selects only the main per-chromosome VCFs:
    chrN.20240612_Freeze2.common_rare_phased.snps.vcf.gz
    Skips any '.nochr' or '.h3a' files.
    """
    out = {}
    for chr_ in chrs:
        pattern = os.path.join(vcf_dir, f"chr{chr_}.20240612_Freeze2.common_rare_phased.snps.vcf.gz")
        matches = glob.glob(pattern)
        if matches:
            out[chr_] = matches[0]
        else:
            print(f"[WARN] chr{chr_}: main .snps.vcf.gz file not found (skipping)")
    return out


def geovar_codes_streaming_fixed(geovar_obj, freq_mat_file_gz):
    """Stream large frequency matrices to avoid memory overflow."""
    assert geovar_obj.bins is not None
    geovar_codes = []
    test_bins = np.array([x[1] for x in geovar_obj.bins])
    with gzip.open(freq_mat_file_gz, 'rt') as f:
        header = f.readline()
        pops = np.array(header.split()[6:])
        geovar_obj.pops = pops
        for line in tqdm(f, desc=f"Streaming {os.path.basename(freq_mat_file_gz)}"):
            maf_vector = np.array(line.split()[6:]).astype(np.float64)
            cur_geovar = np.digitize(maf_vector, test_bins, right=True)
            cur_geovar_code = "".join([str(i) for i in cur_geovar])
            geovar_codes.append(cur_geovar_code)
    geovar_obj.geovar_codes = np.array(geovar_codes)
    geovar_obj.n_variants = geovar_obj.geovar_codes.size
    geovar_obj.n_populations = geovar_obj.pops.size


def process_chromosome(chr_, vcf, pop_df):
    """Main per-chromosome routine."""
    freq_out_gz = os.path.join(RESULTS_DIR, f"chr{chr_}.freq_mat.tsv.gz")
    codes_out_txt = os.path.join(RESULTS_DIR, f"caapa_chr{chr_}.txt")

    # Skip if already done
    if os.path.exists(codes_out_txt):
        print(f"[SKIP] chr{chr_}: Output exists → {codes_out_txt}")
        return

    print(f"[INFO] chr{chr_}: Writing allele frequency table...")
    tmp_plain = os.path.join(RESULTS_DIR, f"chr{chr_}.freq_mat.tsv")

    af_df = vcf_to_freq_table(
        vcf,
        pop_df=pop_df,
        outfile=tmp_plain,
        minor_allele=True
    )

    print(f"[INFO] chr{chr_}: Compressing freq table → {freq_out_gz}")
    with open(tmp_plain, "rb") as fin, gzip.open(freq_out_gz, "wb") as fout:
        for chunk in iter(lambda: fin.read(1024 * 1024), b""):
            fout.write(chunk)
    os.remove(tmp_plain)

    # GeoVar coding
    geovar = GeoVar(bins=BINS)
    geovar_codes_streaming_fixed(geovar_obj=geovar, freq_mat_file_gz=freq_out_gz)

    # Count and write
    codes, counts = geovar.count_geovar_codes()
    df_out = pd.DataFrame({"codes": codes, "counts": counts})
    df_out.to_csv(codes_out_txt, sep=' ', header=False, index=False)
    print(f"[OK] chr{chr_}: Done → {codes_out_txt}")


def main():
    args = parse_args()

    # Load pop panel once
    pop_df = read_pop_panel(POP_PANEL)

    # Determine which chromosomes to run
    CHRS_ALL = [str(i) for i in range(1, 23)] + ["X", "Y"]
    if args.one_chr:
        CHRS = [args.one_chr]
    else:
        CHRS = CHRS_ALL

    chr_vcfs = find_chr_vcfs(VCF_DIR, CHRS)
    if not chr_vcfs:
        raise FileNotFoundError(f"No matching VCFs found in {VCF_DIR}")

    for chr_, vcf in sorted(chr_vcfs.items(), key=lambda kv: (kv[0] not in map(str, range(1,23)), kv[0])):
        print(f"\n[INFO] Processing chr{chr_}: {vcf}")
        process_chromosome(chr_, vcf, pop_df)

    print("\n[DONE] All requested chromosomes processed.")
    print(f"Results in: {RESULTS_DIR}")


if __name__ == "__main__":
    main()

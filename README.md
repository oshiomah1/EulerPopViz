




# Part 1 – Creating the Frequency Matrix with `geovarr.py`

This document describes how to generate a **population frequency matrix** and **GeoVar counts file** from a VCF using the `geovarr.py` script.

---

## 1. Inputs

### 1.1 VCF file

- Input must be a **VCF** file.  
- **BCF is not supported** (convert to VCF first if needed).
- The file should contain the variants for all samples you want to analyze.

### 1.2 Population panel

- A **two-column text file**, separated by a single space.
- Must contain a header line with the columns:

  ```text
  sample pop

### Allele frequency Binning scheme 

In the script, the GeoVar bins object is typically defined as 5 allele frequency bins:

geovar = GeoVar(bins=[(0, 0), (0, 0.01),(0.01, 0.05),(0.05, 0.1),0.1, 1.0)])

## 2. Running geovarr.py

geovarr.py takes your VCF and population panel and produces:

- A frequency matrix file (zipped downstream)

- A GeoVar counts file.

The script expects you to set a few parameters inside the Python file (or via a config, depending on how you've structured it). At minimum, you’ll need to specify:


data_path = "/path/to/data/" # Path to the top level directory containing the VCF (optional, depending on your script)

population_panel = "/path/to/pop_panel.txt" # Path to the population panel

vcf_file = "/path/to/input.vcf" # Path to the input VCF file
 
outfile & freq_file = "/path/to/output/freq_matrix.tsv  # or .txt.gz depending on your script # Path for the output frequency matrix file

freq_mat_file= "/path/to/output/freq_matrix.tsv.gz # Path for the zipped output frequency matrix file

outfile = "/path/to/output/geovar_counts.txt" #  path for other outputs like GeoVar counts


## 3. Understanding Outputs from geovarr.py

#### The frequency matrix file (usually compressed, e.g. freq_matrix.tsv.gz) has columns like: 

| CHR   | SNP     | A1 | A2 | MAC | MAF        | African-American | Afro-Caribbean        | Afro-South_American |
|-------|---------|----|----|-----|------------|------------------|-----------------------|---------------------|
| chr21 | 5031162 | C  | A  | 6   | 0.00086856 | 0.0              | 0.0019305             |         0.13        |

### Column meanings

- CHR – Chromosome.

- SNP – Variant position / ID (depends on how your script defines it).

- A1 / A2 – Alleles (e.g., reference and alternate alleles).

- MAC – Minor Allele Count.

- MAF – Minor Allele Frequency (across all samples).

- Every column after MAF is a population columns – One column per population in your population panel, each containing the frequency of the minor allele in that population.

#### The GeoVar Counts File

The GeoVar counts file summarizes patterns of allele frequencies across populations, binned by frequency ranges.

It looks like this:

```
0000000000001 189451
0000000000010 14886
0000000000011 4740
```

In the GeoVar counts file, Column 1 is a string of digits, one digit per population. The number of digits equals the number of populations in your population panel.

The position of a digit corresponds to a particular population (in the same order as in the panel / frequency matrix). (#adjust this later for accuracy)

The value of each digit (0–4) indicates which frequency bin that population falls into for a given variant.

Interpretation of 0000000000011 4740:

In column 1, there are as many digits as populations (here, 13 digits → 13 populations).

For each variant included in this pattern:

- Populations at positions 1–11 have digit 0 → bin (0, 0) (no minor allele observed).

- Populations at positions 12 and 13 have digit 1 → bin (0, 0.01) (very low frequency).

- The second column (4740) means that 4,740 variants share this exact pattern of allele frequency bins across populations.


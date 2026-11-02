# HBsAg Loss GWAS: Analysis Code

![Language](https://img.shields.io/badge/language-Bash%20%7C%20R%20%7C%20Perl-blue)
![Status](https://img.shields.io/badge/status-research%20code-orange)
![License](https://img.shields.io/badge/license-not%20specified-lightgrey)

Analysis code for a genome-wide association study (GWAS) of genetic factors that influence **HBsAg loss** in chronic hepatitis B (CHB) patients treated with pegylated interferon.

> **中文概览**：本仓库是一项慢性乙肝（CHB）患者聚乙二醇干扰素治疗后 HBsAg 清除相关的全基因组关联分析（GWAS）代码，包括测序比对与变异检测、SNP 质控与关联分析、eQTL 分析以及转录组分析。各分析模块位于独立文件夹中。**本仓库只包含脚本**，原始测序数据和大部分中间输入文件未包含在内（见 [输入数据](#-input-data-not-included)）。

---

## 📑 Table of Contents

- [Overview](#-overview)
- [Pipeline](#-pipeline)
- [Repository Structure](#-repository-structure)
- [Requirements](#-requirements)
- [Installation](#-installation)
- [Usage](#-usage)
- [Input Data (Not Included)](#-input-data-not-included)
- [Known Issues](#-known-issues)
- [Citation](#-citation)

---

## 🔬 Overview

| Stage | Folder | Purpose |
|---|---|---|
| 1. Mapping & variant calling | `Mapping_and_varicant_calling/` | BWA-MEM alignment, duplicate marking, freebayes variant calling |
| 2. SNP QC & association | `SNP_quality_control_and_association_analysis/` | PLINK quality control, population structure (MDS), logistic association, plots |
| 3. eQTL analysis | `eQTL/` | GTEx v7 liver eQTL lookup, liftOver to hg38, eQTL/GWAS overlay, LocusZoom & UCSC data |
| 4. Transcriptome analysis | `transcriptome_analysis_of_HepG2_HepAD38_LM3/` | edgeR differential expression in HepG2, HepAD38 and LM3 cells after TP53BP2 knockdown |

Each folder contains a short readme that maps the scripts to the figures in the paper.

---

## 🧭 Pipeline

```
 FASTQ reads
     │
     ▼
 ┌──────────────────────────┐
 │ 1. BWA-MEM → GATK        │  Mapping_and_Variant_calling.sh
 │    MarkDuplicates        │
 │    freebayes             │
 └────────────┬─────────────┘
              ▼  merge.vcf.gz
 ┌──────────────────────────┐
 │ 2. PLINK QC              │  SNP_quality_control and association analysis.sh
 │    MDS · logistic GWAS   │  MDS_own.R · MDS_merged.R · Manhattan_plot_QQ_plot.R
 └────────────┬─────────────┘
              ▼  assoc_results.assoc
 ┌──────────────────────────┐
 │ 3. eQTL overlay          │  GTEx_v7_merge_GWAS.sh → eQTL_*.R
 └──────────────────────────┘
 ┌──────────────────────────┐
 │ 4. edgeR transcriptome   │  transcriptome_analysis.R
 └──────────────────────────┘
```

---

## 🗂 Repository Structure

```
GWAS_TP53BP2/
├── README.md
├── Mapping_and_varicant_calling/
│   └── Mapping_and_Variant_calling.sh
├── SNP_quality_control_and_association_analysis/
│   ├── SNP_quality_control and association analysis.sh   # main QC + association pipeline
│   ├── SNP_quality_control.sh                             # early filtering (older version)
│   ├── MDS_own.R
│   ├── MDS_merged.R
│   ├── Manhattan_plot_QQ_plot.R
│   ├── Global_SNP_distribution_map.R
│   └── rs7519753_world.txt
├── eQTL/
│   ├── GTEx_v7_merge_GWAS.sh
│   ├── eQTL_protein_gene.R
│   ├── eQTL_GWAS_PPplot.R
│   ├── eQTL_locusZoom_data_prepair.R
│   ├── remap_chip_seq.with_enhancer.R
│   ├── sig.SNP.chr1_223.72_223.92.R
│   ├── id_extract_gtex.pl
│   └── protein_coding.csv
└── transcriptome_analysis_of_HepG2_HepAD38_LM3/
    ├── transcriptome_analysis.R
    └── transcriptome_analysis.readme
```

---

## 🛠 Requirements

**Command-line tools**

- BWA, SAMtools, HTSlib (`bgzip`, `tabix`), BCFtools
- GATK 4 (`SortSam`, `MarkDuplicates`)
- freebayes (including `freebayes-parallel` and `fasta_generate_regions.py`)
- PLINK 1.9 and PLINK 2.0
- UCSC `liftOver`, `wget`, Perl, awk

**R packages**

- CRAN: `tidyverse`, `dplyr`, `magrittr`, `ggplot2`, `ggrepel`, `qqman`, `sf`, `raster`, `sp`, `spData`, `spDataLarge`, `tmaptools`, `tmap`, `devtools`, `VennDiagram`, `dbplyr`, `stringr`, `pheatmap`, `gridExtra`, `factoextra`
- Bioconductor: `edgeR`, `org.Hs.eg.db`, `clusterProfiler`

---

## 📦 Installation

**1. Clone the repository**

```bash
git clone https://github.com/guanguiwensy/GWAS_TP53BP2.git
cd GWAS_TP53BP2
```

**2. Install command-line tools** (conda/mamba recommended; verify package names on Bioconda)

```bash
mamba create -n hbsag_gwas -c conda-forge -c bioconda \
    bwa samtools bcftools htslib gatk4 freebayes plink plink2 ucsc-liftover wget perl
mamba activate hbsag_gwas
```

**3. Install R packages** (R ≥ 4.x recommended)

```r
install.packages(c("tidyverse", "dplyr", "magrittr", "ggplot2", "ggrepel", "qqman",
                   "sf", "raster", "sp", "spData", "spDataLarge", "tmaptools", "tmap",
                   "devtools", "VennDiagram", "dbplyr", "stringr", "pheatmap",
                   "gridExtra", "factoextra"))

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("edgeR", "org.Hs.eg.db", "clusterProfiler"))
```

> `spDataLarge` is not on CRAN. `Global_SNP_distribution_map.R` installs it from GitHub.

> 💡 Software versions were not recorded in the original repository. Please record the versions you use, e.g. `bwa 2>&1 | grep Version` and `plink --version`.

---

## ▶️ Usage

Run the stages in order. Each stage reads files produced by the previous one.

### Stage 1 — Mapping and variant calling

```bash
cd Mapping_and_varicant_calling
bash Mapping_and_Variant_calling.sh \
    -t 16 \
    -i /path/to/bwa_index/hg38.fa \
    -v /path/to/known_variants.vcf.gz \
    -f /path/to/hg38.fa \
    -a /path/to/hg38.fa.fai
```

| Flag | Description |
|:---:|---|
| `-t` | Number of threads |
| `-i` | BWA index prefix |
| `-v` | Input variant VCF passed to freebayes (`-@`) |
| `-f` | Reference FASTA |
| `-a` | Reference `.fai` index (`samtools faidx`) |

Expects paired-end reads named `S<n>.*1.fq.gz` / `S<n>.*2.fq.gz` for samples 1–105 in the working directory. Output: `merge.vcf.gz` (+ index).

### Stage 2 — SNP QC and association analysis

The main pipeline is `SNP_quality_control and association analysis.sh`. It is a list of commands rather than a parameterised script, so set the hard-coded reference paths first.

```bash
cd SNP_quality_control_and_association_analysis
cp "SNP_quality_control and association analysis.sh" run_qc.sh
# edit run_qc.sh: set dbSNP and 1000 Genomes paths for your system
bash run_qc.sh
```

Steps:

1. Keep biallelic SNPs with MAF ≥ 0.01; annotate rsIDs from dbSNP
2. Missingness filters (`--geno`, `--mind`), sex check, autosomes only
3. MAF ≥ 0.05 and Hardy–Weinberg filter (`--hwe 1e-10`)
4. Heterozygosity outlier removal; relatedness filter (PI-HAT > 0.2)
5. Merge with 1000 Genomes Phase 3 (hg38) and compute MDS components
6. Logistic regression with 10 MDS covariates (`--logistic`)

Plots:

```bash
Rscript MDS_own.R
Rscript MDS_merged.R
Rscript Manhattan_plot_QQ_plot.R
Rscript Global_SNP_distribution_map.R   # needs internet access for map data
```

### Stage 3 — eQTL analysis

```bash
cd eQTL
bash GTEx_v7_merge_GWAS.sh                 # GTEx v7 liver eQTLs + hg19→hg38 chain
Rscript eQTL_GWAS_PPplot.R                 # Figure 2A
Rscript eQTL_locusZoom_data_prepair.R      # LocusZoom input (Figure 2D)
Rscript remap_chip_seq.with_enhancer.R     # ReMap peaks annotated with enhancers
```

### Stage 4 — Transcriptome analysis

```bash
cd transcriptome_analysis_of_HepG2_HepAD38_LM3
Rscript transcriptome_analysis.R
```

Outputs: `diffSig_*.csv`, `pca.pdf`, `heatmap.pdf`, `venn.pdf`, `vol.pdf`.

---

## 📥 Input Data (Not Included)

| Stage | Required inputs |
|---|---|
| 1 | Paired-end FASTQ files; reference FASTA with BWA index and `.fai`; known-variant VCF |
| 2 | dbSNP VCF (`dbsnp_146.hg38.vcf.gz`); 1000 Genomes Phase 3 hg38 PLINK set (`1kg.all_hg38`); `sample_keep`; `inversion.txt`; helper scripts `MAF_check.R`, `hwe.R`, `check_heterozygosity_rate.R`, `heterozygosity_outliers_list.R` (called by the pipeline but not committed) |
| 3 | `results_as_0.001.txt` (GWAS summary); LD table `rs7519753_2000000_sig.LD.r2.csv`; ReMap and enhancer position tables. GTEx `Liver.allpairs.txt.gz` is downloaded by the script |
| 4 | `symbol.protein_coding.hepg2.txt`, `symbol.protein_coding.hepad38.txt`, `normalizeExp.g2.txt`, `normalizeExp.ad38.txt`, `lm3_edgerOut.tsv`, `lm3.normalized.csv`, `Type_I_interferon_pathway_genes.txt` |

> ⚠️ Human genetic and clinical data may be subject to access restrictions. Obtain it through the appropriate channels.

---

## ⚠️ Known Issues

These issues were found by reading the code. **Nothing was run end to end**, so each item should be verified before running the pipeline.

### Mapping and variant calling
- `getopts 'h:t:i:v:f:a:'` makes `-h` take an argument, so `-h` does not print usage.
- `-R "@RG\tID:S$i_$k..."` expands the empty variable `$i_`; read-group IDs are not unique across samples. Should be `${i}_$k`.
- `sample=*.merge.marked.bam` does not match MarkDuplicates output (`S<n>.merge.bam.marked.bam`).
- Sample range (`seq 1 105`), memory (`-Xmx60G`) and threads are hard-coded.

### SNP QC and association
- **Lines 4–5:** a line break after `--exclude-uncalled\` joins with `- O z` into `--exclude-uncalled- O z`, an invalid option.
- **Line 2:** `bgzip merge.vcf` fails because the mapping step already produced `merge.vcf.gz`.
- **Line 195:** `--annotate assoc_results.assoc.fisher` reads a file nothing creates, so it fails as written.
- The eQTL scripts read `assoc_results.assoc.fisher.txt` and `assoc_results.assoc.fisher0.0001.txt`, which the pipeline does not produce.
- **Line 67:** the sample removed for low call rate is hard-coded (`echo "1-55 1"`).
- **Line 148:** the 1000 Genomes panel is fetched via `ftp://`, which may no longer be served.
- **Lines 7 and 92:** absolute paths to `/home/guanguiwen/...` are hard-coded.
- `SNP_quality_control.sh` uses an undefined `$anno` and duplicates early steps with different file names; treat it as an older version.

### Population structure
- The pipeline labels Asian samples `ASN`, but `MDS_merged.R` tests for `"EAS"`, so Asian samples are not plotted even though the legend lists them.

### eQTL
- **`GTEx_v7_merge_GWAS.sh` line 11:** `$id=Liver.allpairs.txt` is not an assignment in bash (fails with `command not found`). Should be `id=...`. Later commands depend on it.
- **`eQTL_GWAS_PPplot.R` line 17:** `hg38id` is built from `gene_id` rather than the chromosome, so the merge with GWAS is likely to match nothing.
- `remap_chip_seq.with_enhancer.R` compares coordinates without checking chromosome.
- `eQTL_locusZoom_data_prepair.R` writes `rs7519753_2000000_sig.r2.csv` but reads `rs7519753_2000000_sig.LD.r2.csv`; it selects columns by fixed position; and a trailing comma in `mutate()` may cause an error.
- `sig.SNP.chr1_223.72_223.92.R` relies on `results_as` already existing in the workspace.

### Transcriptome analysis
- **Line 286:** `pdf(width = 1000, height = 1000)` is in inches, not pixels, producing an enormous PDF.
- **Lines 312–313:** the file appears to be GBK-encoded; the label `INF-α` shows as `INF-��` in UTF-8 sessions. Save as UTF-8 or use `"INF-alpha"`.
- `gene_logFC()` takes a third argument that is never passed; the `setdiff(..., y = "TP53BP2_DOWN")` step has no effect.
- `merge_gene` (line 262) is computed but unused.

### Global map
- `tm_*` functions come from `tmap`, which is not loaded (only `tmaptools` is).
- `install.packages()` runs on every execution.
- Row indices `15:18` are hard-coded.

### General
- No license file, and no environment file (conda YAML or renv lockfile), so software versions cannot be recovered from the repository.
- `SNP_quality_control and association analysis.sh.bak` is an older copy of the main pipeline; it is kept as-is.

---

## 📚 Citation

If you use this code, please cite the associated article. *(Citation details to be added.)*

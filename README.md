# Code for the HBsAg loss GWAS study

This repository contains the analysis code for a genome-wide association study (GWAS) of the genetic factors that affect HBsAg loss in chronic hepatitis B (CHB) patients treated with pegylated interferon. The analysis covers read mapping, variant calling, SNP quality control, association analysis, eQTL analysis and transcriptome analysis. Each stage lives in its own folder.

**The repository contains scripts only.** Raw sequencing data, reference genomes, and most intermediate inputs are not included (see [Input data](#input-data-not-included)).

## Repository layout

| Folder | Stage | Main files |
|---|---|---|
| `Mapping_and_varicant_calling/` | Alignment (BWA-MEM), duplicate marking, joint variant calling (freebayes) | `Mapping_and_Variant_calling.sh` |
| `SNP_quality_control_and_association_analysis/` | PLINK QC, population structure (MDS/PCA), association test, Manhattan/QQ plots, MAF map | `SNP_quality_control and association analysis.sh` (main pipeline), `SNP_quality_control.sh` (early filtering only), `MDS_own.R`, `MDS_merged.R`, `Manhattan_plot_QQ_plot.R`, `Global_SNP_distribution_map.R` |
| `eQTL/` | GTEx v7 liver eQTL lookup, liftOver to hg38, eQTL/GWAS overlay plots, LocusZoom and UCSC data preparation | `GTEx_v7_merge_GWAS.sh`, `eQTL_protein_gene.R`, `eQTL_GWAS_PPplot.R`, `eQTL_locusZoom_data_prepair.R`, `remap_chip_seq.with_enhancer.R`, `id_extract_gtex.pl`, `protein_coding.csv` |
| `transcriptome_analysis_of_HepG2_HepAD38_LM3/` | Differential expression (edgeR) in HepG2, HepAD38 and LM3 cells after siTP53BP2/shTP53BP2 knockdown, Venn, PCA and heatmaps | `transcriptome_analysis.R` |

Each folder has a short readme (`*.readme` / `*_readme.txt`) that maps the scripts to the figures in the paper.

## Languages and tools

- **Bash** for pipeline drivers, **awk/sed** for text processing, **Perl** for one helper script, **R** for statistics and plotting.
- **Command-line tools:** BWA, SAMtools, GATK (4.x, `SortSam` and `MarkDuplicates`), freebayes (including `freebayes-parallel` and `fasta_generate_regions.py`), BCFtools, HTSlib (`bgzip`, `tabix`), PLINK 1.9, PLINK 2.0, UCSC `liftOver`, wget.
- **R packages:**
  - CRAN: `tidyverse`, `dplyr`, `magrittr`, `ggplot2`, `ggrepel`, `qqman`, `sf`, `raster`, `sp`, `spData`, `spDataLarge`, `tmaptools`, `tmap`, `devtools`, `VennDiagram`, `dbplyr`, `stringr`, `pheatmap`, `gridExtra`, `factoextra`
  - Bioconductor: `edgeR`, `org.Hs.eg.db`, `clusterProfiler`

## Installation

These steps assume Linux (or macOS for the R parts) with a working `bash`.

### 1. Clone the repository

```bash
git clone https://github.com/guanguiwensy/GWAS_TP53BP2.git
cd GWAS_TP53BP2
```

### 2. Install the command-line tools

The easiest route is a conda/mamba environment from Bioconda and conda-forge. Package names are the usual ones, but check them on Bioconda before you rely on them.

```bash
mamba create -n hbsag_gwas -c conda-forge -c bioconda \
    bwa samtools bcftools htslib gatk4 freebayes plink plink2 ucsc-liftover wget perl
mamba activate hbsag_gwas
```

Versions were not recorded in the repository. Record the versions you actually use (`bwa 2>&1 | grep Version`, `plink --version`, and so on) when you reproduce the analysis.

### 3. Install the R packages

Start R (version 4.x recommended) and run:

```r
install.packages(c("tidyverse", "dplyr", "magrittr", "ggplot2", "ggrepel", "qqman",
                   "sf", "raster", "sp", "spData", "spDataLarge", "tmaptools", "tmap",
                   "devtools", "VennDiagram", "dbplyr", "stringr", "pheatmap",
                   "gridExtra", "factoextra"))

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("edgeR", "org.Hs.eg.db", "clusterProfiler"))
```

`spDataLarge` is not on CRAN. `Global_SNP_distribution_map.R` installs it from GitHub with `devtools::install_github("Nowosad/spDataLarge")`.

## Running the analysis

Run the stages in this order. Each stage writes files that the next stage reads.

### Stage 1: Mapping and variant calling

```bash
cd Mapping_and_varicant_calling
bash Mapping_and_Variant_calling.sh \
    -t 16 \
    -i /path/to/bwa_index/hg38.fa \
    -v /path/to/known_variants.vcf.gz \
    -f /path/to/hg38.fa \
    -a /path/to/hg38.fa.fai
```

| Flag | Meaning |
|---|---|
| `-t` | Threads |
| `-i` | BWA index prefix (the reference FASTA you indexed with `bwa index`) |
| `-v` | Variant VCF passed to freebayes as input variants (`-@`) |
| `-f` | Reference FASTA |
| `-a` | Reference `.fai` index (`samtools faidx`) |

The script expects paired-end FASTQ files named `S<n>.*1.fq.gz` and `S<n>.*2.fq.gz` for samples 1 to 105, in the working directory. It produces `merge.vcf.gz` and its index.

### Stage 2: SNP quality control and association analysis

Stage 2 is the main pipeline file, `SNP_quality_control_and_association_analysis/SNP_quality_control and association analysis.sh`. It is written as a list of commands rather than a script with arguments. Edit the hard-coded paths first (see [Known issues](#known-issues)), then run it from a directory that contains `merge.vcf.gz`:

```bash
cd SNP_quality_control_and_association_analysis
cp "SNP_quality_control and association analysis.sh" run_qc.sh
# edit run_qc.sh: set the dbSNP and 1000 Genomes paths for your system
bash run_qc.sh
```

It runs these steps in order:

1. Keep biallelic SNPs with minor allele frequency ≥ 0.01, annotate rsIDs from dbSNP, convert to PLINK format.
2. Sample and SNP missingness filters (`--geno`, `--mind`), sex check, autosomes only.
3. MAF filter (≥ 0.05) and Hardy-Weinberg filter (`--hwe 1e-10`).
4. Heterozygosity outlier removal, relatedness filter (PI-HAT > 0.2), founders only.
5. Merge with 1000 Genomes Phase 3 (hg38) and compute MDS components for the population-structure plots.
6. Logistic regression with the first 10 MDS components as covariates (`--covar covar_mds.txt --logistic`).

Then generate the plots:

```bash
Rscript MDS_own.R
Rscript MDS_merged.R
Rscript Manhattan_plot_QQ_plot.R
Rscript Global_SNP_distribution_map.R   # needs internet access for the map data
```

### Stage 3: eQTL analysis

```bash
cd eQTL
bash GTEx_v7_merge_GWAS.sh        # downloads GTEx v7 Liver eQTLs and the hg19→hg38 chain file
Rscript eQTL_GWAS_PPplot.R        # Figure 2A
Rscript eQTL_locusZoom_data_prepair.R   # LocusZoom input (Figure 2D)
Rscript remap_chip_seq.with_enhancer.R  # ReMap ChIP-seq peaks annotated with enhancers
```

`eQTL_protein_gene.R` keeps protein-coding genes; `GTEx_v7_merge_GWAS.sh` calls it at the end. `id_extract_gtex.pl` is called by the same script.

### Stage 4: Transcriptome analysis

```bash
cd transcriptome_analysis_of_HepG2_HepAD38_LM3
Rscript transcriptome_analysis.R
```

The script writes its tables (`diffSig_*.csv`) and figures (`pca.pdf`, `heatmap.pdf`, `venn.pdf`, `vol.pdf`) to the working directory.

## Input data (not included)

The following inputs are needed but are not in this repository:

- **Stage 1:** paired-end FASTQ files (`S1_*1.fq.gz`, ... `S105_*2.fq.gz`), the reference genome FASTA and its BWA index and `.fai`, and a known-variant VCF.
- **Stage 2:** dbSNP VCF (`dbsnp_146.hg38.vcf.gz`), the 1000 Genomes Phase 3 hg38 PLINK set (`1kg.all_hg38`), `sample_keep`, `inversion.txt` (PLINK inversion regions), and the helper R scripts `MAF_check.R`, `hwe.R`, `check_heterozygosity_rate.R`, `heterozygosity_outliers_list.R`, which the main pipeline calls but which are not committed.
- **Stage 3:** GTEx v7 `Liver.allpairs.txt.gz` (downloaded by the script), the `results_as_0.001.txt` GWAS summary file, the LD table `rs7519753_2000000_sig.LD.r2.csv`, ReMap and enhancer position tables.
- **Stage 4:** count tables `symbol.protein_coding.hepg2.txt`, `symbol.protein_coding.hepad38.txt`, `normalizeExp.g2.txt`, `normalizeExp.ad38.txt`, `lm3_edgerOut.tsv`, `lm3.normalized.csv`, and `Type_I_interferon_pathway_genes.txt`.

Human data may be subject to access restrictions. Obtain it through the appropriate channels.

## Known issues

These were found by reading the code. Nothing was run end to end, so treat each item as needing a check.

**Mapping and variant calling (`Mapping_and_Variant_calling.sh`)**

- `getopts 'h:t:i:v:f:a:'` makes `-h` take an argument, so `-h` does not print usage.
- `-R "@RG\tID:S$i_$k..."` expands `$i_` (an empty variable), so read-group IDs are not unique across samples. It should be `${i}_$k`.
- `sample=*.merge.marked.bam` does not match the files MarkDuplicates writes (`S<n>.merge.bam.marked.bam`), so freebayes gets an unmatched glob.
- Sample range is hard-coded (`seq 1 105`), and memory (`-Xmx60G`) and threads are fixed.

**SNP QC and association (`SNP_quality_control and association analysis.sh`)**

- Line 4–5: a line break inside `--exclude-uncalled\` followed by `- O z` joins into `--exclude-uncalled- O z`, an invalid bcftools option. Put `-O z` on the same line as the other options.
- Line 2: `bgzip merge.vcf` fails, because the mapping script already compressed `merge.vcf` to `merge.vcf.gz`. Use `merge.vcf.gz` as the input and drop this step, or run it only if the uncompressed file exists.
- Line 195: `plink --annotate assoc_results.assoc.fisher ...` uses a file that nothing creates, so it will fail as written. The `--assoc` run on line 187 produces `assoc_results.assoc`, which `Manhattan_plot_QQ_plot.R` reads.
- The `eQTL` scripts read `assoc_results.assoc.fisher.txt` and `assoc_results.assoc.fisher0.0001.txt`, which the pipeline does not create.
- Line 67: the sample removed for low call rate is hard-coded (`echo "1-55 1"`).
- Line 148: the 1000 Genomes panel is fetched over `ftp://`, which may no longer be served.
- Absolute paths to `/home/guanguiwen/...` are hard-coded (lines 7 and 92).
- `SNP_quality_control.sh` uses an undefined `$anno` and repeats the early steps of the main pipeline with different file names. Treat it as an older version.

**Population structure (`MDS_merged.R`, `Manhattan_plot_QQ_plot.R` and the pipeline)**

- The pipeline labels Asian populations as `ASN` (`sed 's/CHB/ASN/g'` etc.), but `MDS_merged.R` tests for `"EAS"`. Asian samples are therefore not plotted, though the legend shows them.

**eQTL (`eQTL/`)**

- `GTEx_v7_merge_GWAS.sh` line 11: `$id=Liver.allpairs.txt` is not an assignment in bash (it fails with `command not found` and leaves `$id` empty). It should be `id=Liver.allpairs.txt`. Every following command depends on it.
- `eQTL_GWAS_PPplot.R` line 17 builds `hg38id` from `gene_id` (column 2), not the chromosome, so the merge with the GWAS `hg38id` (`chr<CHR>:<BP>_<A2>_<A1>`) is likely to match nothing. Check this before you trust Figure 2A.
- `remap_chip_seq.with_enhancer.R` compares coordinates only and ignores chromosome, so peaks on one chromosome can match enhancers on another.
- `eQTL_locusZoom_data_prepair.R` writes `rs7519753_2000000_sig.r2.csv` but reads `rs7519753_2000000_sig.LD.r2.csv`, and it selects columns by fixed position (`[,c(11:16)]`). It also has a trailing comma in a `mutate()` call that may cause an error in R.
- `sig.SNP.chr1_223.72_223.92.R` uses `results_as` from the workspace. It is not defined in that file.

**Transcriptome analysis (`transcriptome_analysis.R`)**

- Line 286: `pdf("pca.pdf", width = 1000, height = 1000)` sets the page size in inches, not pixels, so the PDF will be enormous.
- Lines 312–313: the file appears to be GBK-encoded (`file` reports ISO-8859 text). The label `"INF-α"` is stored as the GBK bytes for α, which display as `INF-��` in UTF-8 sessions. Save the file as UTF-8, or write the label as `"INF-alpha"`.
- `gene_logFC()` has a third parameter that is never passed, and the `setdiff(..., y = "TP53BP2_DOWN")` step does nothing.
- `merge_gene` (line 262) is computed and not used.

**Global map (`Global_SNP_distribution_map.R`)**

- The `tm_*` functions come from `tmap`, which the script does not load. It loads only `tmaptools`.
- The script runs `install.packages()` every time it runs.
- Row indices `15:18` are hard-coded.

**General**

- The repository has no license file, and the scripts contain no setup or environment files (conda YAML, renv lockfile), so the software versions used in the paper cannot be recovered from the repository.
- `SNP_quality_control and association analysis.sh.bak` is an older copy of the main pipeline. It is kept as-is.

## Citation

If you use this code, please cite the associated article. Details are to be added.

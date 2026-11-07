# Reproduce the 1000 Genomes analyses

Run commands from the repository root. Use R 4.5 or later.

## 1. Packages and tools

Install BEDMatrix from the R console to run either analysis.

```r
install.packages("BEDMatrix")
```

To render either notebook, install `rmarkdown`, as well as `ggplot2` for rendering plots.

```r
install.packages(c("rmarkdown", "ggplot2"))
```

Install `EnsDb.Hsapiens.v75` to prepare the Part 1 annotation and `locuszoomr` to render the notebook's regional plot.

```r
install.packages("BiocManager")
BiocManager::install(c("locuszoomr", "EnsDb.Hsapiens.v75"), update = FALSE)
```

Part 2 uses PLINK 2 for LD clumping against the 1000 Genomes reference panel.
Download [PLINK 2](https://www.cog-genomics.org/plink/2.0/) and place `plink2` (`plink2.exe` for Windows) in a `bin/` folder at the repository root.

To convert raw 23andMe text, install `bcftools` (version 1.24).
Source releases and build instructions are on the [official BCFtools downloads page](https://www.htslib.org/download/).

## 2. Get the input data

Inputs data for the analyses should be placed in `data/raw/`.

| File | Source |
| --- | --- |
|  `LDREF.tar.bz2` | [FUSION download](https://data.broadinstitute.org/alkesgroup/FUSION/LDREF.tar.bz2) |
| `GD462.GeneQuantRPKM.50FN.samplename.resk10.txt.gz` | [Geuvadis on Zenodo](https://zenodo.org/records/15596289#files) |
| `Height.gwas.txt.gz` | [Height GWAS download](https://drive.google.com/file/d/1RWjk49QNZj9zvJHc9X_wyZ51fdy6xQjv/view) |

Leave the `.gz` files compressed. Extract `LDREF.tar.bz2` and place the extracted `LDREF` folder in `data/processed/`.

Check that `data/processed/LDREF/` contains the `.bed`, `.bim`, and `.fam` files for each of `1000G.EUR.1` through `1000G.EUR.22`.

## 3. Run Part 1

The [Bioconductor `EnsDb.Hsapiens.v75` package](https://bioconductor.org/packages/EnsDb.Hsapiens.v75/) supplies the Ensembl release 75 gene annotation in GRCh37 coordinates. Both the script and the notebook create `data/processed/gene_annot_ensembl_grch37.tsv.gz` automatically if it is missing.

Run the analysis across all 22 autosomes.

```powershell
# Run in PowerShell from the repository root.
Rscript --vanilla -e "source('R/part1.R'); run_part1('data/derived/part1_grch37')"
```

Check that `COMPLETE.txt`, `chromosome_summary.tsv`, and `top_hit_per_gene.tsv` appear in `data/derived/part1_grch37`. The per-chromosome association results are saved as `chr*_cis_eqtl.tsv.gz`.

### Reproduce the Part 1 notebook figures

Run `R/part1_runner_notebook.ipynb` with the notebook's working directory set to the repository root.

## 4. Run Part 2

Use a single 23andMe sample you have permission to analyze. Place the GRCh37 PLINK files in your local `data/processed/` with a shared prefix such as `selected_sample_autosomes_GRCh37` and the extensions `.bed`, `.bim`, and `.fam`.

If starting with a 23andMe text file, follow the
[bcftools conversion tutorial](https://samtools.github.io/bcftools/howtos/convert.html)
with an indexed GRCh37 reference FASTA. Use `-Ov -o selected_sample.GRCh37.vcf`
in the tutorial command to write a VCF. Set the VCF path below to your local
file and run the PLINK command.

```powershell
# Run in Powershell from the repository root
& .\bin\plink2.exe --vcf "path/to/selected_sample.GRCh37.vcf" --double-id --chr 1-22 --snps-only just-acgt --max-alleles 2 --make-bed --out "data/processed/selected_sample_autosomes_GRCh37"
```

In `R/part2.Rmd` , set your `.bim` path in `match_height_positions()` and your `.bed` path in `align_height_alleles()`.
In the notebook's final chunk, replace `Tiffany` with your sample ID in the
labels and output filename.

### Reproduce the Part 2 notebook figure

Run `part2_runner_notebook.ipynb`. 

## AI use

OpenAI's Codex was used to help draft and debug R code, prepare figures, and verify results. I directed the analysis and reviewed the changes. See [AI_USE.md](AI_USE.md) for details.

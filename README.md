# Reproduce the 1000 Genomes analyses

Run commands from the repository root. Use R 4.5 or later.

## 1. Packages and tools

Install BEDMatrix from the R console to run either analysis.

```r
install.packages("BEDMatrix")
```

To render either notebook, install `rmarkdown`. The Part 2 notebook also
uses `ggplot2`.

```r
install.packages(c("rmarkdown", "ggplot2"))
```

For the Part 1 notebook's regional plots, install the additional packages.

```r
install.packages("BiocManager")
BiocManager::install(c("locuszoomr", "EnsDb.Hsapiens.v75"), update = FALSE)
```

You can skip installing Pandoc if rendering in RStudio.
Rendering the R makrdown files requires Pandoc. To install Pandoc for PowerShell, download and run the Windows `.msi` installer linked in the [official Pandoc installation guide](https://pandoc.org/installing.html#windows).
The installer adds Pandoc to `PATH`. Open a new PowerShell window before rendering either notebook.

If you need to convert a VCF to PLINK files for Part 2, download [PLINK 2 for Windows](https://www.cog-genomics.org/plink/2.0/) and place `plink2.exe` in `bin/`.

## 2. Get the input data

Put the inputs for the analyses you want to run in `data/raw/`.

| File | Source |
| --- | --- |
|  `LDREF.tar.bz2` | [FUSION download](https://data.broadinstitute.org/alkesgroup/FUSION/LDREF.tar.bz2) |
| `GD462.GeneQuantRPKM.50FN.samplename.resk10.txt.gz` | [Geuvadis on Zenodo](https://zenodo.org/records/15596289#files) |
| `gene_annot.txt.gz` | Get the supplied annotation from the mentor |
| `Height.gwas.txt.gz` | [Height GWAS download](https://drive.google.com/file/d/1RWjk49QNZj9zvJHc9X_wyZ51fdy6xQjv/view) |

Leave the `.gz` files compressed. Extract `LDREF.tar.bz2` and place the extracted `LDREF` folder in `data/processed/`.

Check that `data/processed/LDREF/` contains the `.bed`, `.bim`, and `.fam` files for each of `1000G.EUR.1` through `1000G.EUR.22`.

## 3. Run Part 1

Run the analysis across all 22 autosomes.

```powershell
# Run in PowerShell from the repository root.
Rscript --vanilla -e "source('R/part1.R'); run_part1('data/derived/part1')"
```

We selected the output folder `data/derived/part1`. Check that `COMPLETE.txt`, `chromosome_summary.tsv`, and `top_hit_per_gene.tsv` appear in that folder. The per-chromosome association results are saved as `chr*_cis_eqtl.tsv.gz`.

### Reproduce the Part 1 notebook figures

Open [`R/part1.Rmd`](R/part1.Rmd) in RStudio and click **Knit**.

## 4. Run Part 2

Use a single 23andMe sample you have permission to analyze. Place the GRCh37 PLINK files in your local `data/processed/` with a shared prefix such as `selected_sample_autosomes_GRCh37` and the extensions `.bed`, `.bim`, and `.fam`.

If starting with a 23andMe text file, make a GRCh37 VCF using the [conversion tutorial](https://samtools.github.io/bcftools/howtos/convert.html).
Replace the VCF path below with your local file. Then, run the command.

```powershell
# Run in powershell from the repository root
& .\bin\plink2.exe --vcf "path/to/selected_sample.GRCh37.vcf" --double-id --chr 1-22 --snps-only just-acgt --max-alleles 2 --make-bed --out "data/processed/selected_sample_autosomes_GRCh37"
```

In `R/part2.Rmd` , set your `.bim` path in `match_height_positions()` and your `.bed` path in `align_height_alleles()`.
In the notebook's final chunk, replace `Tiffany` with your sample ID in the
labels and output filename.

### Reproduce the Part 2 notebook figure

Open [`R/part2.Rmd`](R/part2.Rmd) in RStudio and click **Knit**.

## AI use

OpenAI's Codex was used to help draft and debug R code, prepare figures, and verify results. I directed the analysis and reviewed the changes. See[AI_USE.md](AI_USE.md) for details.

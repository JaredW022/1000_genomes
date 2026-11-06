read_height_gwas <- function(file = "data/raw/Height.gwas.txt.gz") {
  gwas <- read.delim(file)

  # Keep SNPs with stronger evidence of association with height and reliable imputation quality
  keep <- gwas$P <= 0.001 & gwas$INFO > 0.8

  candidates <- gwas[
    which(keep),
    c("CHR", "BP", "SNP", "A1", "A2", "P", "OR")
  ]
  candidates$weight <- log(candidates$OR)

  candidates
}

read_height_bim <- function(file) {
  read.table(
    file,
    col.names = c("CHR", "SNP", "CM", "BP", "A1", "A2"),
    colClasses = c("integer", "character", "NULL", "integer",
                   "character", "character")
  )
}

# Returns one row per SNP containing GWAS effect, corresponding .bim, and reference SNP ID and alleles 
match_height_positions <- function(gwas) {
  tiffany <- read_height_bim("data/processed/Tiffany_autosomes_GRCh37.bim")

  reference_files <- file.path(
    "data/processed/LDREF",
    paste0("1000G.EUR.", 1:22, ".bim")
  )
  
  # Combines 22 Chromosom BIM files to SNP reference table
  reference <- NULL
  for (file in reference_files) {
    chromosome_table <- read_height_bim(file)
    reference <- rbind(reference, chromosome_table)
  }

  # Create chromosome-position keys to match GWAS weights with the
  # corresponding SNPs in Tiffany's and the 1000 Genomes genotype files.
  positions <- list(
    gwas = paste(gwas$CHR, gwas$BP, sep = ":"),
    tiffany = paste(tiffany$CHR, tiffany$BP, sep = ":"),
    reference = paste(reference$CHR, reference$BP, sep = ":")
  )

  # Drop repeated positions so match() cannot select the wrong SNP
  # when multiple variants share the same chromosome and base-pair position.
  unique_positions <- lapply(positions, function(x) {
    x[!duplicated(x) & !duplicated(x, fromLast = TRUE)]
  })
  shared <- Reduce(intersect, unique_positions)

  gwas_rows <- match(shared, positions$gwas)
  tiffany_rows <- match(shared, positions$tiffany)
  reference_rows <- match(shared, positions$reference)

  data.frame(
    gwas[gwas_rows, ],
    tiffany_bim_row = tiffany_rows,
    tiffany_a1 = tiffany$A1[tiffany_rows],
    tiffany_a2 = tiffany$A2[tiffany_rows],
    reference_snp = reference$SNP[reference_rows],
    reference_a1 = reference$A1[reference_rows],
    reference_a2 = reference$A2[reference_rows]
  )
}


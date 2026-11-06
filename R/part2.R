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

# Keep SNPs whose reference alleles match the GWAS and for which Tiffany
# has a nonmissing genotype we can use to count the GWAS effect allele.
align_height_alleles <- function(matches) {
  tiffany_genotypes <- BEDMatrix::BEDMatrix(
    "data/processed/Tiffany_autosomes_GRCh37.bed",
    simple_names = TRUE
  )
  # Allele order determined by the .bim file
  tiffany_a1_copies <- as.numeric(
    tiffany_genotypes[1, matches$tiffany_bim_row, drop = FALSE]
  )

  reference_matches <-
    (matches$reference_a1 == matches$A1 & matches$reference_a2 == matches$A2) |
    (matches$reference_a1 == matches$A2 & matches$reference_a2 == matches$A1)

  tiffany_matches <-
    (matches$tiffany_a1 == matches$A1 & matches$tiffany_a2 == matches$A2) |
    (matches$tiffany_a1 == matches$A2 & matches$tiffany_a2 == matches$A1)

  # A dot can mark an unlisted allele in a homozygous call.
  tiffany_homozygous <- matches$tiffany_a1 == "." &
    tiffany_a1_copies == 0 &
    (matches$tiffany_a2 == matches$A1 | matches$tiffany_a2 == matches$A2)

  keep <- reference_matches & !is.na(tiffany_a1_copies) &
    (tiffany_matches | tiffany_homozygous)

  aligned <- matches[keep, , drop = FALSE]
  a1_copies <- tiffany_a1_copies[keep]

  # If the effect allele is in Tiffany's A2 column, count copies as 2 minus A1.
  aligned$tiffany_effect_copies <- ifelse(
    aligned$tiffany_a2 == aligned$A1,
    2 - a1_copies,
    a1_copies
  )

  aligned
}

# Clump SNPs to keep only the most significant (lowest p-valiue)
# This is done to account for associations with nearby SNPs
# Returns only the SNPs which were selected through clumping (lead SNPs)
clump_height_chromosome <- function(aligned, chromosome, output_dir) {
  variants <- aligned[aligned$CHR == chromosome, , drop = FALSE]

  reference_prefix <- file.path(
    "data/processed/LDREF", paste0("1000G.EUR.", chromosome)
  )
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  output_prefix <- file.path(output_dir, paste0("chr", chromosome))
  clump_input <- paste0(output_prefix, "_clump_input.tsv")

  write.table(
    data.frame(ID = variants$reference_snp, P = variants$P),
    clump_input,
    sep = "\t", quote = FALSE, row.names = FALSE
  )

  # Clumping done with plink2's clump command
  status <- system2("bin/plink2.exe", c(
    "--bfile", shQuote(reference_prefix),
    "--clump", shQuote(clump_input),
    "--clump-p1", "0.001",
    "--clump-r2", "0.1",
    "--clump-kb", "250",
    "--out", shQuote(output_prefix)
  ))

  if (status != 0) {
    stop("PLINK clumping failed on chromosome ", chromosome)
  }

  clumps <- read.delim(paste0(output_prefix, ".clumps"))
  lead_rows <- match(clumps$ID, variants$reference_snp)
  variants[lead_rows, , drop = FALSE]
}
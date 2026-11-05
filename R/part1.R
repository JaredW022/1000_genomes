# Save tab-separated results without quotes or row numbers.
write_part1_table <- function(x, file) {
  write.table(x, file, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}

# Load expression and annotation data, match genes, and define cis windows.
read_part1_inputs <- function(config, include_expression = TRUE) {
  # Read the expression table, preserving sample IDs.
  expression <- read.delim(
    gzfile(config$expression_file),
    check.names = FALSE
  )

  # Sample IDs skip the first 4 columns, which are metadata
  sample_ids <- names(expression)[5:466]

  # Inspection mode keeps only gene metadata.
  if (!include_expression) {
    expression <- expression[, 1:4, drop = FALSE]
  }

  keep_rows <- !is.na(match(expression$Chr, config$chromosomes))
  expression <- expression[keep_rows, , drop = FALSE]

  # Strip version suffix from Ensembl gene IDs to match annotation file.
  expression$ensembl_id <- sub("[.][0-9]+$", "", expression$TargetID)

  annotation <- read.delim(gzfile(config$annotation_file), colClasses = "character")
  required <- c("ID", "CHR", "START", "STOP", "SYM", "TYPE")

  # Keep protein-coding genes on selected chromosomes and the required annotation columns.
  genes <- annotation[annotation$TYPE == "protein_coding" &
                        annotation$CHR %in% config$chromosomes, required, drop = FALSE]

  genes$expression_row <- match(genes$SYM, expression$ensembl_id)
  genes$has_expression <- !is.na(genes$expression_row)

  # TODO: Confirm annotation/BIM genome builds with the professor before the full run.
  genes$gene_start <- as.numeric(genes$START)
  genes$gene_end <- as.numeric(genes$STOP)

  # genes$window_start defines the -500kb cis-SNP for each gene. 
  genes$window_start <- pmax(1, genes$gene_start - config$cis_window_bp)

  # genes$window_end defines the +500kb cis-SNP for each gene.
  genes$window_end <- genes$gene_end + config$cis_window_bp

  list(expression = expression, genes = genes, sample_ids = sample_ids)
}


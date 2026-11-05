part1_expression_file <- "data/raw/GD462.GeneQuantRPKM.50FN.samplename.resk10.txt.gz"
part1_annotation_file <- "data/raw/gene_annot.txt.gz"
part1_genotype_dir <- "data/processed/LDREF"
part1_cis_window_bp <- 500000

# Save tab-separated results without quotes or row numbers.
write_part1_table <- function(x, file) {
  write.table(x, file, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}

# Load expression and annotation data, match genes, and define cis windows.
read_part1_inputs <- function(chromosomes = 1:22) {
  # Read the expression table, preserving sample IDs.
  expression <- read.delim(
    gzfile(part1_expression_file),
    check.names = FALSE
  )

  # Sample IDs skip the first 4 columns, which are metadata
  sample_ids <- names(expression)[5:466]

  keep_rows <- !is.na(match(expression$Chr, chromosomes))
  expression <- expression[keep_rows, , drop = FALSE]

  # Strip version suffix from Ensembl gene IDs to match annotation file.
  expression$ensembl_id <- sub("[.][0-9]+$", "", expression$TargetID)

  annotation <- read.delim(gzfile(part1_annotation_file), colClasses = "character")
  required <- c("ID", "CHR", "START", "STOP", "SYM", "TYPE")

  # Keep protein-coding genes on selected chromosomes and the required annotation columns.
  keep_genes <- annotation$TYPE == "protein_coding" & !is.na(match(annotation$CHR, chromosomes))
  genes <- annotation[keep_genes, required, drop = FALSE]

  genes$expression_row <- match(genes$SYM, expression$ensembl_id)
  genes$has_expression <- !is.na(genes$expression_row)

  # TODO: Confirm annotation/BIM genome builds with the professor before the full run.
  genes$gene_start <- as.numeric(genes$START)
  genes$gene_end <- as.numeric(genes$STOP)

  # genes$window_start defines the -500kb cis-SNP for each gene. 
  genes$window_start <- pmax(1, genes$gene_start - part1_cis_window_bp)

  # genes$window_end defines the +500kb cis-SNP for each gene.
  genes$window_end <- genes$gene_end + part1_cis_window_bp

  list(expression = expression, genes = genes, sample_ids = sample_ids)
}

fit_cis_regressions <- function(genotypes, expression, variants) {
  stopifnot(nrow(genotypes) == length(expression), ncol(genotypes) == nrow(variants)) # Not necessary
  
  X <- genotypes
  y <- expression
  observed <- rep(nrow(X), ncol(X))
  mean_by_snp <- colMeans(X, na.rm = TRUE)

  # Keep SNPs where the number of copies of the effect allele differs across individuals.
  has_variance <- apply(X, 2, var) > 0
  if (!any(has_variance)) return(NULL)
  # Build feature vector
  X <- X[, has_variance, drop = FALSE]
  variants <- variants[has_variance, , drop = FALSE]
  observed <- observed[has_variance]
  mean_by_snp <- mean_by_snp[has_variance]

  # Fit to linear regression model for each SNP
  # Extracts the coefficients, standard errors, and p-values.
  fits <- lapply(seq_len(ncol(X)), function(j) {
    fit <- lm(y ~ X[, j])
    coefficients <- summary(fit)$coefficients
    c(
      intercept = coefficients["(Intercept)", "Estimate"],
      beta = coefficients["X[, j]", "Estimate"],
      standard_error = coefficients["X[, j]", "Std. Error"],
      p_value = coefficients["X[, j]", "Pr(>|t|)"]
    )
  })

  # Builds matrix of estimates for each SNP
  estimates <- matrix(
    NA_real_,
    nrow = length(fits),
    ncol = 4,
    dimnames = list(NULL, c("intercept", "beta", "standard_error", "p_value"))
  )

  for (j in seq_along(fits)) {
    estimates[j, ] <- fits[[j]]
  }

  result <- data.frame(snp = variants$SNP, snp_bp = variants$BP,
                       effect_allele = variants$A1, other_allele = variants$A2,
                       beta = estimates[, "beta"],
                       standard_error = estimates[, "standard_error"],
                       p_value = estimates[, "p_value"],
                       intercept = estimates[, "intercept"],
                       n_samples = length(y), n_genotyped = observed,
                       mean_effect_allele_dosage = mean_by_snp, stringsAsFactors = FALSE)
  if (!nrow(result)) {
    return(NULL)
  }
  result
}


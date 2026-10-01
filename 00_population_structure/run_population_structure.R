#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L) y else x
value_after <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i)) return(default)
  if (i == length(args)) stop("Missing value after ", flag)
  args[[i + 1L]]
}

script_path <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
script_dir <- dirname(normalizePath(script_path, mustWork = TRUE))
root <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
vcf_file <- value_after("--vcf", file.path(root, "shared_data", "8pop_poolseq.vcf.gz"))
out_dir <- value_after("--outdir", file.path(script_dir, "results"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

required <- c("poolfstat", "FactoMineR", "factoextra", "pheatmap", "RColorBrewer", "ggplot2")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install required R packages: ", paste(missing, collapse = ", "))

labels <- c(
  AAI = "Dole7", AAK = "Starfood Holland", AAL = "Kingsect Belgium",
  AAQ = "Wild Serbian", AAU = "Pronutrix+ Belgium", AAW = "Papek Czechia",
  AAX = "Mystik Canada", ABE = "USA"
)

message("Reading eight-population unannotated VCF: ", vcf_file)
snp <- poolfstat::vcf2pooldata(vcf.file = vcf_file, poolsizes = rep(80, 8), verbose = TRUE)
if (!setequal(snp@poolnames, names(labels))) {
  stop("VCF samples do not match the expected eight IDs: ", paste(names(labels), collapse = ", "))
}
snp <- poolfstat::pooldata.subset(
  snp, pool.index = match(names(labels), snp@poolnames), min.maf = 0.01,
  cov.qthres.per.pool = c(0.10, 0.90), verbose = TRUE
)
snp@poolnames <- unname(labels)

freq <- snp@refallele.readcount / snp@readcoverage
colnames(freq) <- snp@poolnames
pca <- FactoMineR::PCA(t(freq), scale.unit = TRUE, graph = FALSE)
utils::write.csv(pca$ind$coord, file.path(out_dir, "PCA_coordinates.csv"))
for (axes in list(c(1, 2), c(3, 4))) {
  p <- factoextra::fviz_pca_ind(
    pca, axes = axes, col.ind = "cos2", gradient.cols = c("#EFC000", "#0073C2"),
    repel = TRUE, labelsize = 3
  ) + ggplot2::theme_classic(base_size = 8)
  ggplot2::ggsave(file.path(out_dir, sprintf("PCA_dim%d_%d.pdf", axes[1], axes[2])), p,
                  width = 3.2, height = 2.2)
}

fst <- poolfstat::compute.fstats(
  snp, nsnp.per.bjack.block = 0, computeDstat = FALSE, computeF3 = FALSE,
  computeF4 = FALSE, output.pairwise.fst = TRUE, output.pairwise.div = TRUE,
  computeQmat = FALSE, return.F2.blockjackknife.samples = FALSE,
  return.F4.blockjackknife.samples = FALSE, verbose = TRUE
)
utils::write.table(fst@pairwise.fst, file.path(out_dir, "pairwise_FST.tsv"),
                   sep = "\t", quote = FALSE, col.names = NA)
grDevices::pdf(file.path(out_dir, "pairwise_FST_heatmap.pdf"), width = 4.5, height = 4)
pheatmap::pheatmap(fst@pairwise.fst,
  fontsize = 7, col = grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, "Oranges"))(25))
grDevices::dev.off()

fstat <- poolfstat::compute.fstats(
  snp, nsnp.per.bjack.block = 20000, computeDstat = TRUE, computeQmat = TRUE,
  return.F2.blockjackknife.samples = TRUE, verbose = TRUE
)
he <- data.frame(Population = rownames(fstat@heterozygosities), fstat@heterozygosities,
                 row.names = NULL, check.names = FALSE)
utils::write.table(he, file.path(out_dir, "expected_heterozygosity.tsv"),
                   sep = "\t", quote = FALSE, row.names = FALSE)

farm_names <- unname(labels[c("AAI", "AAK", "AAL", "AAU", "AAW", "AAX")])
farm <- he[match(farm_names, he$Population), ]
wild <- he[he$Population == labels[["AAQ"]], ]
contrasts <- data.frame(
  Population = farm$Population, He_farm = farm$Estimate, He_wild = wild$Estimate,
  Delta_He = farm$Estimate - wild$Estimate,
  SE_delta = sqrt(farm$`bjack s.e.`^2 + wild$`bjack s.e.`^2)
)
contrasts$Z <- contrasts$Delta_He / contrasts$SE_delta
contrasts$P_raw <- 2 * stats::pnorm(-abs(contrasts$Z))
contrasts$P_Holm <- stats::p.adjust(contrasts$P_raw, method = "holm")
utils::write.table(contrasts, file.path(out_dir, "farm_vs_wild_heterozygosity.tsv"),
                   sep = "\t", quote = FALSE, row.names = FALSE)

message("Population-structure analysis complete: ", out_dir)

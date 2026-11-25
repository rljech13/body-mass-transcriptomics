#!/usr/bin/env Rscript
# monaLisa::calcBinnedMotifEnrR on one species.
# Two bins: story-set peaks vs the other proximal peaks.
# Default weights match Homer: GC bins and k-mers up to 3.
# Test is one-sided Fisher (enrichment in the set). Score threshold
# is the package default min.score = 10, not FIMO's 1e-4.

suppressPackageStartupMessages({
  library(monaLisa)
  library(Biostrings)
  library(TFBSTools)
  library(BiocParallel)
  library(SummarizedExperiment)
})

args <- commandArgs(trailingOnly = TRUE)
opt <- list()
i <- 1
while (i <= length(args)) {
  key <- sub("^--", "", args[[i]])
  opt[[key]] <- args[[i + 1]]
  i <- i + 2
}
stopifnot(all(c("meme", "fasta", "species", "geneset-dir", "sets", "out", "cores") %in% names(opt)))

set_tab <- read.delim(opt$sets, stringsAsFactors = FALSE)
sets <- setNames(set_tab$title, set_tab$set)

read_meme_pwms <- function(path) {
  lines <- readLines(path, warn = FALSE)
  bg <- c(A = 0.25, C = 0.25, G = 0.25, T = 0.25)
  motifs <- list()
  i <- 1
  n <- length(lines)
  while (i <= n) {
    if (!startsWith(lines[[i]], "MOTIF ")) {
      i <- i + 1
      next
    }
    parts <- strsplit(trimws(lines[[i]]), "\\s+")[[1]]
    id <- parts[[2]]
    alt <- if (length(parts) >= 3) parts[[3]] else id
    i <- i + 1
    while (i <= n && !grepl("letter-probability", lines[[i]], fixed = TRUE)) i <- i + 1
    if (i > n) stop("truncated motif ", id)
    hdr <- lines[[i]]
    nsites <- as.integer(sub(".*nsites=\\s*([0-9]+).*", "\\1", hdr))
    if (is.na(nsites) || nsites < 1) nsites <- 100L
    i <- i + 1
    rows <- list()
    while (i <= n && grepl("^\\s*[0-9]", lines[[i]])) {
      rows[[length(rows) + 1]] <- as.numeric(strsplit(trimws(lines[[i]]), "\\s+")[[1]])
      i <- i + 1
    }
    M <- do.call(rbind, rows)
    if (ncol(M) != 4) stop("bad row in ", id)
    prof <- t(M) * nsites
    rownames(prof) <- c("A", "C", "G", "T")
    pfm <- TFBSTools::PFMatrix(
      ID = id, name = alt, profileMatrix = prof,
      tags = list(nsites = as.character(nsites))
    )
    motifs[[length(motifs) + 1]] <- TFBSTools::toPWM(
      pfm, type = "log2probratio", pseudocounts = 0.8, bg = bg
    )
  }
  do.call(TFBSTools::PWMatrixList, motifs)
}

message("monaLisa ", as.character(packageVersion("monaLisa")),
        " reading ", opt$meme)
pwms <- read_meme_pwms(opt$meme)
message("matrices ", length(pwms))
seqs <- readDNAStringSet(opt$fasta)
gene <- sub("\\|.*", "", names(seqs))
message(opt$species, " sequences ", length(seqs))

bpp <- MulticoreParam(workers = as.integer(opt$cores), progressbar = FALSE)
rows <- list()
for (stem in names(sets)) {
  gs <- readLines(file.path(opt[["geneset-dir"]], paste0(stem, ".txt")))
  gs <- gs[nzchar(gs)]
  n_set <- sum(gene %in% gs)
  message(stem, " peaks ", n_set)
  if (n_set < 2) next
  bins <- factor(ifelse(gene %in% gs, "set", "rest"), levels = c("set", "rest"))
  se <- calcBinnedMotifEnrR(
    seqs = seqs,
    bins = bins,
    pwmL = pwms,
    background = "otherBins",
    test = "fisher",
    maxKmerSize = 3L,
    min.score = 10,
    BPPARAM = bpp,
    verbose = FALSE
  )
  if (!("set" %in% colnames(se))) stop("no set column for ", stem)
  nlp <- assay(se, "negLog10P")[, "set"]
  p <- 10^(-nlp)
  log2enr <- assay(se, "log2enr")[, "set"]
  fg <- assay(se, "sumForegroundWgtWithHits")[, "set"]
  bgw <- assay(se, "sumBackgroundWgtWithHits")[, "set"]
  rd <- as.data.frame(rowData(se))
  motif_id <- if ("motif.id" %in% names(rd)) rd$motif.id else rownames(se)
  motif_name <- if ("motif.name" %in% names(rd)) rd$motif.name else rownames(se)
  rows[[stem]] <- data.frame(
    set = stem,
    title = unname(sets[[stem]]),
    species = opt$species,
    motif_id = motif_id,
    tf = toupper(motif_name),
    p = as.numeric(p),
    log2enr = as.numeric(log2enr),
    fg_wgt = as.numeric(fg),
    bg_wgt = as.numeric(bgw),
    n_pos = n_set,
    n_neg = length(seqs) - n_set,
    stringsAsFactors = FALSE
  )
}

out <- do.call(rbind, rows)
dir.create(dirname(opt$out), recursive = TRUE, showWarnings = FALSE)
write.table(out, opt$out, sep = "\t", quote = FALSE, row.names = FALSE)
message("wrote ", opt$out, " rows ", nrow(out))

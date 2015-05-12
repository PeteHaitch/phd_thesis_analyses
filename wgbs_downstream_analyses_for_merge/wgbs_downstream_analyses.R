# Number methylation patterns observed at 4-tuples for Lister data

library(MethylationTuples)
x <- readRDS("Lister_4_tuples_strand_collapsed.rds")

# assayNames <- function(x) {
#   names(assays(x, withDimnames = FALSE))
# }

z <- lapply(colnames(x), function(sn, x, i) {
  Reduce("+", mclapply(assayNames(x), function(an, x, i) {
    zz <- assay(x, an) >= i
    zz[is.na(zz)] <- FALSE
    zz
  }, x = x[, sn], i = i, mc.cores = 3))
}, x = x, i = 2L)

sapply(z, table)

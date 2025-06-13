geneMatrix <- readr::read_tsv("~/projects/dcgna_t2d/workflow/auxiliary_data/geneMatrix.tsv_new.gz")
genesets <- readr::read_tsv("~/projects/dcgna_t2d/workflow/auxiliary_data/genesets_20220621.tsv")


readr::write_rds(geneMatrix,"inst/extdata/geneMatrix.rds")



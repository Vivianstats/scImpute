###############################################################################
## Zeisel data

## expression_mRNA_17-Aug-2014.txt
## can be downloaded from
## https://storage.googleapis.com/linnarsson-lab-www-blobs/blobs/cortex/expression_mRNA_17-Aug-2014.txt

data_dir = "./"
x = read.table(paste0(data_dir, "expression_mRNA_17-Aug-2014.txt"), skip = 11,
                header = FALSE, stringsAsFactors = FALSE,
                row.names = 1)
x = x[, -1]
cellnames = read.table(paste0(data_dir, "expression_mRNA_17-Aug-2014.txt"), 
                        skip = 7, nrows = 1, row.names = 1, 
                        stringsAsFactors = FALSE)
colnames(x) = cellnames[-1]

labels = read.table(paste0(data_dir, "expression_mRNA_17-Aug-2014.txt"), 
                        skip = 1, nrows = 1, row.names = 1, 
                        stringsAsFactors = FALSE)
labels = unlist(labels)
table(labels)


matching = c("1"="Interneurons", "2"="S1-Pyramidal", "3"="CA1-Pyramidal",
              "4"="Oligodendrocytes", "5"="Microglia", "6"="Endothelial",
              "7" = "Astrocytes", "8" = "Ependymal",
              "9"="Mural")
labels = matching[as.character(labels)]

saveRDS(x, paste0(data_dir, "rerun/zeisel_raw.rds"))
saveRDS(labels, paste0(data_dir, "rerun/zeisel_label9.rds"))


### level2 classes
labels = read.table(paste0(data_dir, "expression_mRNA_17-Aug-2014.txt"), 
                    skip = 9, nrows = 1, row.names = 1, 
                    stringsAsFactors = FALSE)
labels = unlist(labels)
table(labels)
saveRDS(labels, paste0(data_dir, "rerun/zeisel_label47.rds"))




write.csv(x, paste0(data_dir, "rerun/zeisel_raw.csv"), quote = FALSE)


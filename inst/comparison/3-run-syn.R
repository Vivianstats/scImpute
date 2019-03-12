data_dir = "./"
count_path = paste0(data_dir, "zeisel_samp.rds")

count = readRDS(count_path)
write.csv(count, file = paste0(data_dir, "zeisel_samp.csv"),
          row.names = TRUE, col.names = TRUE)
############ run scImpute
#library(devtools)
#install_github("Vivianstats/scImpute", ref = "895c262") #v0.0.3
library(scImpute)

count_path = paste0(data_dir, "zeisel_samp.csv")

out_dir = paste0(data_dir, "rerun/samp/")
dir.create(out_dir)
scimpute(count_path = count_path, infile = "csv", outfile = "csv",
         Kcluster = 9,
         out_dir = out_dir, drop_thre = 0.5, ncores = 36)

count = read.csv(paste0(out_dir, "scimpute_count.csv"), row.names = 1)
saveRDS(count, file = paste0(out_dir, "zeisel_scimpute_k9.rds"))


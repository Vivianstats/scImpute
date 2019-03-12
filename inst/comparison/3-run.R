data_dir = "./"
count_path = paste0(data_dir, "rerun/zeisel_raw.csv")

############ run scImpute
library(devtools)
install_github("Vivianstats/scImpute", ref = "895c262") #v0.0.3
library(scImpute)

out_dir = paste0(data_dir, "rerun/")
dir.create(out_dir)
scimpute(count_path = count_path, infile = "csv", outfile = "csv",
         Kcluster = 9,
         out_dir = out_dir, drop_thre = 0.5, ncores = 36)

count = read.csv(paste0(out_dir, "scimpute_count_k9.csv"), row.names = 1)
saveRDS(count, file = paste0(out_dir, "zeisel_scimpute_k9.rds"))

# ############ run SAVER
# devtools::install_github("mohuangx/SAVER", ref="b64a077") #v1.0.0

library(SAVER)
library(doParallel)

count_path = paste0(data_dir, "rerun/zeisel_raw.rds")

cl = makeCluster(35, outfile = "")
registerDoParallel(cl)

dat = readRDS(count_path)
out = saver(dat)

saveRDS(out, file = paste0(data_dir, "rerun/zeisel_saver.rds"))


# # ############ run MAGIC
# dat = read.csv(paste0(data_dir, "rerun/zeisel_magic.csv"))
# rownames(dat) = dat[, 1]
# dat = dat[,-1]
# saveRDS(dat, paste0(data_dir, "rerun/zeisel_magic.rds"))


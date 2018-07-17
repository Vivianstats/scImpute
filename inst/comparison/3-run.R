data_dir = "./"
count_path = paste0(data_dir, "rerun/zeisel_raw.rds")

############ run scImpute
library(devtools)
install_github("Vivianstats/scImpute", ref = "v0.0.8")
library(scImpute)

out_dir = "./rerun/"
dir.create(out_dir)
scimpute(count_path = count_path, infile = "rds", outfile = "rds",
         Kcluster = 7,
         out_dir = out_dir, drop_thre = 0.5, ncores = 35)


# ############ run SAVER
# devtools::install_github("mohuangx/SAVER")
library(SAVER)
library(doParallel)

cl = makeCluster(35, outfile = "")
registerDoParallel(cl)

saver1 = saver(dat, pred.genes = 1:4993, pred.genes.only = TRUE)
saveRDS(saver1, "~/Dropbox/Rpkgs-dev/scimpute_dev/diagnosis/SAVER-paper/data/rerun/zeisel_saver1.rds")
saver2 = saver(dat, pred.genes = 4994:9986, pred.genes.only = TRUE)
saveRDS(saver2, "~/Dropbox/Rpkgs-dev/scimpute_dev/diagnosis/SAVER-paper/data/rerun/zeisel_saver2.rds")
saver3 = saver(dat, pred.genes = 9987:14979, pred.genes.only = TRUE)
saveRDS(saver3, "~/Dropbox/Rpkgs-dev/scimpute_dev/diagnosis/SAVER-paper/data/rerun/zeisel_saver3.rds")
saver4 = saver(dat, pred.genes = 14980:19972, pred.genes.only = TRUE)
saveRDS(saver4, "~/Dropbox/Rpkgs-dev/scimpute_dev/diagnosis/SAVER-paper/data/rerun/zeisel_saver4.rds")

saver.all = combine.saver(list(saver1, saver2, saver3, saver4))
saveRDS(saver.all, "./rerun/zeisel_saver_combine.rds")

# ############ run MAGIC
dat = read.csv(paste0(data_dir, "rerun/zeisel_magic.csv"))
rownames(dat) = dat[, 1]
dat = dat[,-1]
saveRDS(dat, paste0(data_dir, "rerun/zeisel_magic.rds"))


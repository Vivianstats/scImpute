library(Rtsne)
library(ggplot2)
library(parallel)
library(ClusterR)
library(tidyr)
library(dplyr)
library(gridExtra)
library(grid)
#library(kernlab)

plot_dir = "./plots/"
data_dir = "./rerun/"


### data in Huang et al.
### downloaded from https://www.dropbox.com/sh/ri6fa3mbhvgapqk/AADwOzHfiCcLSqYnX9CTyd7_a?dl=0

### labels
labels = readRDS(paste0(data_dir, "rerun/zeisel_label9.rds"))
labels_samp = labels[match(colnames(count_samp), colnames(count_raw))]

info = readRDS(paste0(data_dir, "rerun/samp/fig2d_tsne.rds"))
ident = info[[2]]
labels_Huang = ident[[4]][[1]]
write.table(table(labels_samp, labels_Huang), file = paste0(data_dir, "rerun/samp/label-matrix.txt"))


#####################################################
### tSNE

methods = c("syn", "scImpute", "MAGIC", "SAVER")
name_appends = c("samp", "samp_scimpute", "samp_magic", "samp_saver")
names(name_appends) = methods

for(method in methods){
  set.seed(1234)
  print(method)
  dim = 2
  count_raw = readRDS(paste0(data_dir, "rerun/samp/zeisel_", name_appends[method], ".rds"))
  if(method == "SAVER") count_raw = count_raw$estimate
  count = log10(count_raw + 1)
  tsne = Rtsne(t(count), dims = dim)$Y
  saveRDS(tsne, file = paste0(data_dir, "rerun/samp/zeisel-", method, "-tsne", dim, ".rds"))
  gc()
}


### tSNE
data = lapply(methods, function(method){
  tsne = readRDS(file = paste0(data_dir, "rerun/samp/zeisel-", method,
                               "-tsne2", ".rds"))
  pdata = data.frame(tSNE1 = tsne[,1], tSNE2 = tsne[,2], type = labels_samp)
  pdata$method = method
  return(pdata)
})
data = Reduce(rbind, data)
data$method = factor(data$method, levels = c("syn", "scImpute", "MAGIC", "SAVER"))
gt = ggplot(data, aes(x = tSNE1, y = tSNE2, color = type)) +
  geom_point(alpha = 0.8, cex = 0.8) + facet_wrap(~method, nrow = 1) +
  theme_bw() +
  theme(strip.background = element_blank(),
        legend.position = "bottom",
        text = element_text(size=12))
ggsave(paste0(plot_dir,"Fig-samp.pdf"), gt, width = 11, height = 4)


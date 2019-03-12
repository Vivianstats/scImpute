library(Rtsne)
library(ggplot2)
library(parallel)
library(ClusterR)
library(tidyr)
library(dplyr)
library(gridExtra)
library(grid)
#library(kernlab)
###source("~/Dropbox/Rpkgs-dev/scimpute_dev/diagnosis/SAVER-paper/v2/comparison/supp.R")

plot_dir = "./plots/"
data_dir = "./rerun/"


#####################################################
### data characteristics
count_raw = readRDS(paste0(data_dir, "rerun/zeisel_raw.rds"))
### data in Huang et al.
### downloaded from https://www.dropbox.com/sh/ri6fa3mbhvgapqk/AADwOzHfiCcLSqYnX9CTyd7_a?dl=0
count_samp = readRDS(paste0(data_dir, "zeisel_samp.rds"))

datas =c("raw", "SAVER-paper")
summary = lapply(1:2, function(i){
  if(i == 1){count = count_raw}
  if(i == 2){count = count_samp}
  mean = log10(rowMeans(count)+1)
  sd = log10(apply(count, 1, sd)+1)
  zero = rowSums(count == 0)/ncol(count)
  da = data.frame(mean, sd, zero)
  da$data = datas[i]
  return(da)
})
summary = Reduce(rbind, summary)
g1 = ggplot(summary, aes(x = mean, y = sd, color = data)) +
  geom_point(alpha = 0.6, cex = 0.2) +
  xlab("log10(mean + 1)") + ylab("log10(sd + 1)") +
  scale_color_manual(values = c("#999999", "#CC0C00")) +
  # scale_y_continuous(labels = scaleFUN) +
  theme_bw() +
  theme(legend.position = "none",
        text = element_text(size=14))
g2 = ggplot(summary, aes(x = mean, y = zero, color = data)) +
  geom_point(alpha = 0.6, cex = 0.2) +
  xlab("log10(mean + 1)") + ylab("zero fraction") +
  scale_color_manual(values = c("#999999", "#CC0C00")) +
  theme_bw() +
  theme(legend.position = "none",
        text = element_text(size=14))
g = arrangeGrob(g1,g2,nrow = 2)
ggsave(paste0(plot_dir, "Fig2a.pdf"), g, width = 3, height = 5)

g3 = ggplot(summary, aes(x = mean, fill = data)) +
  geom_density(alpha = 0.6) + xlim(c(0,1.5)) +
  xlab("log10(mean + 1)") +
  scale_fill_manual(values = c("#999999", "#CC0C00")) +
  theme_bw() +
  theme(legend.position = "none",
        text = element_text(size=14))
g4 = ggplot(summary, aes(x = sd, fill = data)) +
  geom_density(alpha = 0.6) + xlim(c(0,1.5)) +
  xlab("log10(sd + 1)") +
  scale_fill_manual(values = c("#999999", "#CC0C00")) +
  theme_bw() +
  theme(legend.position = "none",
        text = element_text(size=14))
g5 = ggplot(summary, aes(x = zero, fill = data)) +
  geom_density(alpha = 0.6) +
  xlab("zero fraction") +
  scale_fill_manual(values = c("#999999", "#CC0C00")) +
  theme_bw() +
  theme(legend.position = "none",
        text = element_text(size=14))
gg = arrangeGrob(g3,g4,g5,nrow = 3)
ggsave(paste0(plot_dir, "Fig2b.pdf"), gg, width = 3, height = 5)



#####################################################
### tSNE

methods = c("Raw", "scImpute", "MAGIC", "SAVER")
name_appends = c("raw", "scimpute_k9", "magic", "saver")
names(name_appends) = methods
labels = readRDS(paste0(data_dir, "rerun/zeisel_label9.rds"))

# for(method in methods){
#   set.seed(1234)
#   print(method)
#   dim = 2
#   count_raw = readRDS(paste0(data_dir, "rerun/zeisel_", name_appends[method], ".rds"))
#   if(method == "SAVER") count_raw = count_raw$estimate
#   count = log10(count_raw + 1)
#   tsne = Rtsne(t(count), dims = dim)$Y
#   saveRDS(tsne, file = paste0(data_dir, "rerun/zeisel-", method, "-tsne", dim, ".rds"))
#   gc()
# }


### tSNE
data = lapply(methods, function(method){
  tsne = readRDS(file = paste0(data_dir, "rerun/zeisel-", method,
                               "-tsne2", ".rds"))
  pdata = data.frame(tSNE1 = tsne[,1], tSNE2 = tsne[,2], type = labels)
  pdata$method = method
  return(pdata)
})
data = Reduce(rbind, data)
data$method = factor(data$method, levels = c("Raw", "scImpute", "MAGIC", "SAVER"))
gt = ggplot(data, aes(x = tSNE1, y = tSNE2, color = type)) +
  geom_point(alpha = 0.8, cex = 0.8) + facet_wrap(~method, nrow = 1) +
  theme_bw() +
  theme(strip.background = element_blank(),
        legend.position = "bottom",
        text = element_text(size=12))
ggsave(paste0(plot_dir,"Fig2e.pdf"), gt, width = 11, height = 4)



##########################################################
### clustering

B = 100
J = 3005

for(kk in c(9,47)){
  print(kk)
  dim = 10
  temp_res = mclapply(1:B, function(b){
    set.seed(b)
    if(b %% 20 == 0) print(b)
    ind = sample(1:J, J, replace = TRUE)
    val = sapply(methods, function(method){
      mat = readRDS(paste0(data_dir, "rerun/zeisel-", method, "-pca.rds"))
      mat = mat[ind, 1:dim]
      truel = labels[ind]
      clusts = hclust(dist(mat), method = "median")
      Clabel = cutree(clusts, kk)
      Clabel = as.numeric(Clabel)
      v = sapply(c("jaccard_index", "adjusted_rand_index", "purity", "nmi"), function(x){
        external_validation(as.numeric(factor(truel)), Clabel, method = x)
      })
      return(v)
    })
    gc()
    mat = as.data.frame(val)
    colnames(mat) = methods
    mat$measure = c("Jaccard index", "adjusted Rand index", "purity", "nmi")
    mat = mat %>% gather(metric, value, -measure)
    return(mat)
  }, mc.cores = 36)
  da = Reduce(rbind, temp_res)
  

  val = sapply(methods, function(method){
    mat = readRDS(paste0(data_dir, "rerun/zeisel-", method, "-pca.rds"))
    mat = mat[, 1:dim]
    clusts = hclust(dist(mat), method = "median")
    Clabel = cutree(clusts, kk)
    Clabel = as.numeric(Clabel)
    v = sapply(c("jaccard_index", "adjusted_rand_index", "purity", "nmi"), function(x){
      external_validation(as.numeric(factor(labels)), Clabel, method = x)
    })
    return(v)
  })
  mat = as.data.frame(val)
  colnames(mat) = methods
  mat$measure = c("Jaccard index", "adjusted Rand index", "purity", "nmi")
  mat = mat %>% gather(metric, value, -measure)
  sd = sapply(1:nrow(mat), function(i){
    val = filter(da, measure == mat[i,"measure"] & metric == mat[i,"metric"])
    return(sd(val$value))
  })
  mat$sd = sd

  mat$metric = factor(mat$metric, levels = c("Raw", "scImpute", "MAGIC", "SAVER"))
  gc = ggplot(mat, aes(x = metric, y = value, fill = metric)) +
    geom_bar(stat = "identity", width = .7) + 
    facet_wrap(~measure, nrow = 4, scales = "free") +
    geom_errorbar(aes(ymin = value - sd, ymax = value+sd), width = .2) +
    theme_bw() + ylab("") + ylim(0,NA) +
    theme(strip.background = element_rect(colour="white", fill="white"),
          text = element_text(size=12),
          axis.text.x=element_text(size=8),
          axis.ticks.x=element_blank(),
          legend.position = "none") +
    scale_fill_manual(values = c("#999999", "#56B4E9", "#CC79A7", "#E69F00"))
  ggsave(paste0(plot_dir, "Fig2c-kk", kk, ".pdf"), gc, width = 2.5, height = 5)
}







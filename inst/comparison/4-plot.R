library(Rtsne)
library(ggplot2)
library(parallel)
library(kernlab)
library(ClusterR)
library(tidyr)
library(gridExtra)
library(grid)

plot_dir = "./plots/"
data_dir = "./rerun/"

#####################################################
### data characteristics
count_raw = readRDS(paste0(data_dir, "zeisel_raw.rds"))
### data in Huang et al.
### downloaded from https://www.dropbox.com/sh/ri6fa3mbhvgapqk/AADwOzHfiCcLSqYnX9CTyd7_a?dl=0
count_samp = readRDS("./zeisel_samp.rds")

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
  geom_point(alpha = 0.6, cex = 0.5) + 
  xlab("log10(mean + 1)") + ylab("log10(sd + 1)") +
  scale_color_manual(values = c("#999999", "#E69F00")) +
  scale_y_continuous(labels = scaleFUN) +
  theme_bw() +
  theme(legend.position = "none",
        text = element_text(size=14))
g2 = ggplot(summary, aes(x = mean, y = zero, color = data)) +
  geom_point(alpha = 0.6, cex = 0.5) + 
  xlab("log10(mean + 1)") + ylab("zero fraction") +
  scale_color_manual(values = c("#999999", "#E69F00")) +
  theme_bw() +
  theme(legend.position = "none",
        text = element_text(size=14))
g = arrangeGrob(g1,g2,nrow = 2)
g = arrangeGrob(g1,g2,nrow = 2)
ggsave(paste0(plot_dir, "Fig2a.pdf"), g, width = 4, height = 5)

#####################################################
### comparison

methods = c("Raw", "scImpute", "MAGIC", "SAVER")
name_appends = c("raw", "scimpute_k7", "magic", "saver_combine")
names(name_appends) = methods
dims = 2
labels = readRDS(paste0(data_dir, "zeisel_label.rds"))

# for(method in methods){
#   set.seed(1234)
#   print(method)
#   mclapply(dims, function(dim){
#     print(dim)
#     count_raw = readRDS(paste0(data_dir, "zeisel_", name_appends[method], ".rds"))
#     if(method == "SAVER") count_raw = count_raw$estimate
#     count = log10(count_raw + 1)
#     tsne = Rtsne(t(count), dims = dim)$Y
#     saveRDS(tsne, file = paste0(data_dir, "zeisel-", method, "-tsne", dim, ".rds"))
#     gc()
#   }, mc.cores = 2)
# }


### tSNE
data = lapply(methods, function(method){
  tsne = readRDS(file = paste0(data_dir, "zeisel-", method,
                               "-tsne2", ".rds"))
  pdata = data.frame(tSNE1 = tsne[,1], tSNE2 = tsne[,2], type = labels)
  pdata$method = method
  return(pdata)
})
data = Reduce(rbind, data)
data$method = factor(data$method, levels = c("Raw", "scImpute", "MAGIC", "SAVER"))
g3 = ggplot(data, aes(x = tSNE1, y = tSNE2, color = type)) +
  geom_point(alpha = 0.8, cex = 0.8) + facet_wrap(~method, nrow = 1) +
  theme_bw() +
  theme(strip.background = element_blank(),
        legend.position = "bottom",
        text = element_text(size=12)) 
ggsave(paste0(plot_dir,"Fig2b.pdf"), g3, width = 12, height = 4)


### evaluation
# dim = 2
# temp_res = sapply(methods, function(method){
#   set.seed(1234)
#   mat = readRDS( paste0(data_dir, "zeisel-", method, "-tsne", dim, ".rds"))
#   Clabel = kmeans(mat, centers = kk, nstart = 10)$cluster
#   saveRDS(Clabel,paste0(data_dir, "label-", method, "-tsne", dim, "-kk", kk, ".rds") )
#   gc()
#   return(0)
# })


mat = sapply(methods, function(method){
  set.seed(1234)
  Clabel = readRDS(paste0(data_dir, "label-", method, "-tsne2-kk7.rds") )
  Clabel = as.numeric(Clabel)
  sapply(c("jaccard_index", "adjusted_rand_index", "purity", "nmi"), function(x){
    external_validation(as.numeric(factor(labels)), Clabel, method = x)
  })
})
mat = data.frame(mat)
colnames(mat) = methods
mat$measure = c("Jaccard index", "adjusted Rand index", "purity", "nmi")

mat = mat %>% gather(key = method, value = "value",-(measure))
mat$method = factor(mat$method, levels = c("Raw", "scImpute", "MAGIC", "SAVER"))
g4 = ggplot(mat, aes(x = method, y = value, fill = method)) +
  geom_bar(stat = "identity", position=position_dodge()) +
  facet_wrap(~measure, nrow = 2) +
  theme_bw() + ylab("") + ylim(0,1) + xlab("") +
  theme(strip.background = element_rect(colour="white", fill="white"),
        text = element_text(size=12),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  scale_fill_manual(values = c("#999999", "#56B4E9", "#CC79A7", "#E69F00")) 
ggsave(paste0(plot_dir, "Fig2c.pdf"), g4, width = 6, height = 5)







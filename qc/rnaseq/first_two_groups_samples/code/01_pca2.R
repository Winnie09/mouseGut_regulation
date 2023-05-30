rm(list=ls())
## first group: TPM
tb1 = read.csv('/Users/wenpinhou/Dropbox/aravinda/data/rna/counts.raw.original_samples.csv', row.names = 1)
summary(tb1[,3])
summary(colSums(tb1))

meta = read.csv('/Users/wenpinhou/Dropbox/aravinda/data/meta/postnatal_meta_rnaseq_wh.csv')
rownames(meta) = meta[,2] 
str(meta)

## second group: counts
tb2 = read.csv('/Users/wenpinhou/Dropbox/aravinda/data/rna/counts_and_Rsq_for_new_data_wh.csv', row.names = 1)
str(tb2)


tb3 = read.csv('/Users/wenpinhou/Dropbox/aravinda/data/rna/HSCR_samples_wh.csv', row.names = 1)
tb3 = tb3[, 4:ncol(tb3)]
str(tb3)

gnome = readRDS('/Users/wenpinhou/Dropbox/resource/grcm38_geneid_genename_genelength.rds')
str(gnome)
int = intersect(rownames(tb2), gnome[,2])
tb2 = tb2[int, ]
tb3 = tb3[int, ]
tb1 = tb1[int, ]


tpm3 <- function(counts,len) {
  x <- counts/len
  return(t(t(x)*1e6/colSums(x)))
}


summary(tb1)

tb1.norm = log2(tpm3(tb1, gnome[match(int, gnome[,2]), 3]) + 1)
tb2.norm = log2(tpm3(tb2, gnome[match(int, gnome[,2]), 3]) + 1)
tb3.norm = log2(tpm3(tb3, gnome[match(int, gnome[,2]), 3]) + 1)
summary(colSums(tb2.norm))
summary(colSums(tb3.norm))
summary(colSums(tb1.norm))


tb = cbind(tb1.norm, tb2.norm)
colnames(tb)


tmp = as.vector(as.matrix(tb))
pdf('/Users/wenpinhou/Dropbox/aravinda/qc/first_two_groups_samples/plot/hist_expr_temporal.pdf', width = 5, height = 3.5)
hist(tmp[tmp!=0], breaks = 100, xlab = 'Gene expression (log2(TPM+1))', main = 'Histogram of non-zero gene expression values')
abline(v = 2, col = 'red')
dev.off()


expr = tb[rowSums(tb>2) > 1,] # 14862, 24
expr = expr[, c(21:24, 4:6, 15:16, 7:9, 17:18, 10:12, 19:20, 13:14, 1:3)]
time = meta[colnames(expr), 3]
str(expr)
colnames(expr)


tmp = as.vector(as.matrix(expr))
pdf('/Users/wenpinhou/Dropbox/aravinda/qc/first_two_groups_samples/plot/hist_expr_temporal_after_filtering.pdf', width = 5, height = 3.5)
hist(tmp[tmp!=0], breaks = 100, xlab = 'Gene expression (log2(TPM+1))', main = 'Histogram of non-zero gene expression values')
dev.off()
  
source('/Users/wenpinhou/Dropbox/resource/myfunc/01_function.R')
ls()
pr = PCA(expr)
str(pr) # 24, 9
saveRDS(pr, '/Users/wenpinhou/Dropbox/aravinda/qc/first_two_groups_samples/plot/pr.rds')

library(ggplot2)
library(RColorBrewer)
source('/Users/wenpinhou/Dropbox/resource/ggplot_theme.R')

setdiff(rownames(pr), meta[,2])

pdf('/Users/wenpinhou/Dropbox/aravinda/qc/first_two_groups_samples/plot/RNAseq_pca_timepoint.pdf', width = 2.8, height = 1.8)
pd = data.frame(pc1 = pr[,1], pc2 = pr[,2], sample = rownames(pr), time = meta[rownames(pr), 3])
ggplot(data = pd, aes(x = pc1, y = pc2, color = time)) + 
  geom_point(stroke = 0, size = 0.6) +
  xlab('Principal component 1') + ylab('Principal component 2') +
  ggtitle('timepoint') + 
  scale_color_brewer(palette = 'Dark2')
dev.off()


pdf('/Users/wenpinhou/Dropbox/aravinda/qc/first_two_groups_samples/plot/RNAseq_pca_samplegroups.pdf', width = 2.6, height = 1.8)
pd = data.frame(pc1 = pr[,1], pc2 = pr[,2], sample = rownames(pr), time = ifelse(rownames(pr) %in% colnames(tb1), 'group 1', 'group 2'))
ggplot(data = pd, aes(x = pc1, y = pc2, color = time)) + 
  geom_point(stroke = 0, size = 0.6) +
  xlab('Principal component 1') + ylab('Principal component 2') +
  ggtitle('Sample groups Rebecca sent')
dev.off()

colnames(meta)
pdf('/Users/wenpinhou/Dropbox/aravinda/qc/first_two_groups_samples/plot/RNAseq_pca_libprep.pdf', width = 2.6, height = 1.8)
pd = data.frame(pc1 = pr[,1], pc2 = pr[,2], sample = rownames(pr), time = meta[rownames(pr), 11])
ggplot(data = pd, aes(x = pc1, y = pc2, color = time)) + 
  geom_point(stroke = 0, size = 0.6) +
  xlab('Principal component 1') + ylab('Principal component 2') +
  ggtitle('submit for lib prep')
dev.off()


pdf('/Users/wenpinhou/Dropbox/aravinda/qc/first_two_groups_samples/plot/RNAseq_pca_pairedend.pdf', width = 2.6, height = 1.8)
pd = data.frame(pc1 = pr[,1], pc2 = pr[,2], sample = rownames(pr), time = meta[rownames(pr), 14])
ggplot(data = pd, aes(x = pc1, y = pc2, color = time)) + 
  geom_point(stroke = 0, size = 0.6) +
  xlab('Principal component 1') + ylab('Principal component 2') +
  ggtitle('paired end reads')
dev.off()




###### 
library(splines)
## assume within any one time point, the variance is the same as any other time points, i.e. across time points the variances are shared. 
pval = sapply(1:nrow(expr), function(i){
  v = as.numeric(as.vector(expr[i,]))
  
  mod <- lm(v~time)
  ano <- anova(mod)
  ano['Pr(>F)'][1,1]
})
hist(pval)


fdr = p.adjust(pval, method='fdr')
saveRDS(fdr, '/Users/wenpinhou/Dropbox/aravinda/qc/first_two_groups_samples/res/RNAseq_fdr.rds')

hist(fdr)
sum(fdr<0.05)
str(fdr)
names(pval)  <- names(fdr) <- rownames(expr)
par(mfrow=c(1,2))
plot(as.numeric(expr[which(fdr>0.5)[2],]))
plot(as.numeric(expr[which(fdr<0.05)[2],]))
dev.off()

expr.sig = expr[fdr < 0.05, ] ## 12786
str(expr.sig)

par(mfrow=c(1,1))
expr.sig.scale = t(apply(expr.sig, 1, scale))
str(expr.sig.scale)
colnames(expr.sig.scale) <- colnames(expr.sig)
  
d = dist(expr.sig.scale)
str(d)
hclu = hclust(d)
str(hclu)
clu = cutree(hclu, k = 10)
str(clu)
table(clu)
saveRDS(clu, '/Users/wenpinhou/Dropbox/aravinda/qc/first_two_groups_samples/res/dynamic_genes_cluster_10clu.rds')



#####
res.tb = data.frame(gene_name = names(clu), cluster = clu, 
                    pvalue = pval[names(clu)],
                    FDR = fdr[names(clu)],
                    stringsAsFactors = FALSE)
res.tb = res.tb[order(res.tb[,3]), ]
head(res.tb)
saveRDS(res.tb, '/Users/wenpinhou/Dropbox/aravinda/qc/first_two_groups_samples/res/dynamic_genes_10clu.rds')


write.table(res.tb, '/Users/wenpinhou/Dropbox/aravinda/qc/first_two_groups_samples/res/dynamic_genes_10clu.csv', row.names = F, sep = ',')


library(pheatmap)
library(RColorBrewer)

expr.sig.scale = expr.sig.scale[names(clu)[order(as.numeric(clu))], ]
n.bak = names(clu)
clu = as.character(clu)
names(clu) = n.bak
rowann = data.frame(cluster = clu[rownames(expr.sig.scale)])
rownames(rowann) = rownames(expr.sig.scale)
clu.col = colorRampPalette(rev(brewer.pal(11, 'Set3')))(length(unique(clu)))
names(clu.col) = unique(clu)
annotation_colors = list(cluster = clu.col)
pdf('/Users/wenpinhou/Dropbox/aravinda/qc/first_two_groups_samples/plot/hm_diffgene.pdf', width = 4.5, height = 4.5)
pheatmap(expr.sig.scale, show_rownames = F, cluster_cols = FALSE, cluster_rows = FALSE, annotation_row = rowann, annotation_colors = annotation_colors)
dev.off()


###################
## GO Enrichment
###################

diffgeneList <-
          sapply(sort(unique(clu)), function(i) {
            #########
            names(clu)[clu == i]
          })
str(diffgeneList)

library(topGO)
resList <- lapply(diffgeneList, function(diffgene) {
  allgene <- rownames(tb)
  gl = diffgene
  back = allgene
  ## gl <- sub(sep, '', diffgene)
  ## back <- sub(sep, '', allgene)
  geneList <- factor(as.integer(back %in% gl))
  names(geneList) <- back
  suppressMessages({
    GOdata <-
      new(
        "topGOdata",
        ontology = "BP",
        allGenes = geneList,
        geneSel = function(a) {
          a
        },
        annot = annFUN.org,
        mapping = "org.Mm.eg.db", ## org.Hs.eg.db
        ID = "Symbol"
      )
    resultFisher <-
      topGO::runTest(GOdata, algorithm = "classic", statistic = "fisher")
    sigres <-
      topGO::GenTable(
        GOdata,
        classicFisher = resultFisher,
        topNodes = length(resultFisher@score),
        orderBy = "classicFisher",
        numChar = 1000
      )
  })
  sigres$classicFisher[sigres$classicFisher == "< 1e-30"] <- 0
  sigres <- sigres[sigres$Annotated >= 10,]
  sigres$FDR <-
    stats::p.adjust(sigres$classicFisher, method = "fdr")
  ptcount <- 0
  fc <-
    ((sigres[, "Significant"] + ptcount) / (sum(GOdata@allScores[GOdata@feasible] ==
                                                  1) + ptcount)) / ((sigres[, "Annotated"] + ptcount) / (sum(GOdata@feasible) +
                                                                                                           ptcount))
  sigres <- data.frame(sigres, FC = fc)
  sigres <- sigres[order(sigres$FDR,-sigres$FC),]
})
names(resList) <- sort(unique(clu))
saveRDS(resList, '/Users/wenpinhou/Dropbox/aravinda/qc/first_two_groups_samples/plot/GOres.rds')


library(reshape2)
source('/Users/wenpinhou/Dropbox/trajectory_variability/package/Lamian/R/plotGOEnrich.R')
pdf('/Users/wenpinhou/Dropbox/aravinda/qc/first_two_groups_samples/plot/hm_diffgene_GO.pdf', width = 9, height = 10)
plotGOEnrich(goRes = resList, n = 10, fdr.cutoff = 0.05, fc.cutoff = 2)
dev.off()




str(expr.sig)
expr.sig.clu = rowsum(expr.sig, clu)
expr.sig.clu = expr.sig.clu[order(as.numeric(rownames(expr.sig.clu))), ]
rownames(expr.sig.clu)

tmp = t(sapply(seq(1,10), function(i){
  colSums(expr.sig[names(clu)[clu==i], ])
}))

str(expr.sig.clu)
str(tmp)

expr.sig.clu[2,]
tmp[2,]
expr.sig.clu[3,]
tmp[3,]
expr.sig.clu[4,]
tmp[4,]

expr.sig.clu.scale = t(apply(expr.sig.clu, 1, scale))

colnames(expr.sig.clu.scale) <- colnames(expr.sig.clu)
rownames(expr.sig.clu.scale) <- paste0('cluster', rownames(expr.sig.clu.scale))
str(expr.sig.clu.scale)
colnames(expr.sig.clu.scale) = time

pd = melt(expr.sig.clu.scale)
str(pd)

pdf('/Users/wenpinhou/Dropbox/aravinda/qc/first_two_groups_samples/plot/diffgene_cluster_pattern.pdf', width = 7, height = 4)
ggplot(data = pd, aes(x = Var2, y = value)) +
  geom_point(size = 0.5) +
  geom_smooth() + 
  facet_wrap(~Var1, nrow = 2) +
  xlab('Time points')  +
  ylab('Standardized mean log2(TPM + 1)') +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()


## ======================
## plot individual genes
## ======================
gs = c('Oct4', 'Sox2', 'Klf4', 'Myc', 'Pou5f1', 'Gata1', 'Gata2', 'Runx1', 'Bcl11a', 'Smad1', 'Pparg', 'Nr1h2')
gs = gs[toupper(gs) %in% toupper(rownames(expr.sig))]

pd <- sapply(gs, function(g){
  g = toupper(g)
  v = as.numeric(expr.sig[which(toupper(rownames(expr.sig)) == g), ])
})
str(pd)  
head(pd)
rownames(pd) = time
gs.fdr = sapply(colnames(pd), function(i) {
  round(res.tb[i, 4],3)
})
gs.fdr[gs.fdr == 0] = 0.001
colnames(pd) = paste0(colnames(pd), '(FDR<', gs.fdr, ')')

pd = reshape2::melt(pd)
head(pd)
colnames(pd) = c('time', 'gene', 'expression')

pdf('/Users/wenpinhou/Dropbox/aravinda/qc/first_two_groups_samples/plot/diffgene_dev.pdf', width = 7, height = 4)
ggplot(data = pd, aes(x = time, y = expression)) +
  geom_point(size = 0.5) +
  # geom_smooth(se = FALSE, method = "gam", formula = y ~ x) + 
  geom_smooth() + 
  facet_wrap(~gene, nrow = 2) +
  xlab('Time points')  +
  ylab('log2(TPM + 1)') +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()



## plot top significant genes
gs = rownames(res.tb)[1:20]
pd <- sapply(gs, function(g){
  g = toupper(g)
  v = as.numeric(expr.sig[which(toupper(rownames(expr.sig)) == g), ])
})
rownames(pd) = time
gs.fdr = sapply(colnames(pd), function(i) {
  round(res.tb[i, 4],3)
})
gs.fdr[gs.fdr == 0] = 0.001
colnames(pd) = paste0(colnames(pd), '(FDR<', gs.fdr, ')')

pd = reshape2::melt(pd)
head(pd)
colnames(pd) = c('time', 'gene', 'expression')
pdf('/Users/wenpinhou/Dropbox/aravinda/qc/first_two_groups_samples/plot/diffgene_top.pdf', width = 12, height = 6)
ggplot(data = pd, aes(x = time, y = expression)) +
  geom_point(size = 0.5) +
  geom_smooth(se = FALSE, method = "gam", formula = y ~ s(log(x))) + 
  facet_wrap(~gene, nrow = 3) +
  xlab('Time points')  +
  ylab('log2(TPM + 1)') +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()





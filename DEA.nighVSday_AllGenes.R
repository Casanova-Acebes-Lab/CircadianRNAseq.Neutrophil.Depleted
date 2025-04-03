
# Loading Libraries

library(limma)
library(edgeR)
library(Glimma)
library(gridExtra)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(R.utils)
library(clusterProfiler)
library(biomaRt)
library(cowplot)
library(DESeq2)
library(ggsci)
library(org.Mm.eg.db)


# Config file to define folder and file paths

config<-read.table(file="config.tsv", header = TRUE, fill=TRUE, sep = "\t")

keys <- config$key
values <- config$value
names(values) <- keys
###############
wkdir <- values["wkdir"]
datadir <- values["datadir"]  


source("functions.R") #Needed to call functions stored in file functions.R 

# Reading matrix and metadata

# Reading expression matrix
matrix <- read.table(paste0(datadir,"/normalized_counts.csv"),
                     sep = ",",
                     header = T,
                     check.names = F, 
                     row.names = 1)

matrix <- matrix[,-1]

# Metadata

metadata <- read.table(paste0(datadir,"/metadata_new.txt"),
                     header = T)


# New group for Day and Night samples 
metadata$cycle <- NA
metadata$cycle[1:length(metadata$cycle)/2] <- "Day"
metadata$cycle[19:length(metadata$cycle)] <- "Night"

# Factoring

metadata$condition <- as.factor(metadata$condition)
metadata$time <- as.factor(metadata$time)
metadata$cycle <- as.factor(metadata$cycle)

# Removing duplicated genes by variance. Keeping those with higher variance.

vari <- apply(matrix[1: ncol(matrix)-1 ],1, var)
matrix <- matrix[order(-vari), ]
matrix <- matrix[!duplicated(matrix$ext_gene), ]# check this

rownames(matrix) <- matrix$ext_gene
matrix<- matrix[1: ncol(matrix)-1 ]

# Writing deduplicated matrix to disk
write.table(matrix,paste0(datadir,"/RawMatrix.tsv"),
sep = "\t",
col.names = NA)

# Quality control

## Filtering genes with low expression
thresh <- matrix > 0.5
keep <- rowSums(thresh) >= 2
summary(keep)
counts.keep <- matrix[keep,]
dim(counts.keep)


## Deseq Object
dgeObj <- DGEList(counts.keep)
### have a look at dgeObj
dgeObj

## Checking library sizes
pdf(paste0(wkdir,"/plots/libSize.pdf"), width=12,height=12)
barplot(dgeObj$samples$lib.size, names=colnames(dgeObj), las=2, main="Barplot of library sizes")
abline(h=20e6, lty=2)
dev.off()


## Get log2 counts per million to use if necessary
logcounts <- cpm(dgeObj, log=TRUE)
### Check distributions of samples using boxplots
pdf(paste0(wkdir,"/plots/libNorm.pdf"), width=12,height=12)
boxplot(logcounts, 
        xlab="", 
        ylab="Log2 counts per million",
        las=2)
### Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(logcounts), col="blue")
dev.off()


### Writing log normalized matrix
write.table(logcounts,paste0(datadir,"/logcounts.tsv"),
sep = "\t",
col.names = NA)



###########  Differential expression Analysis with Limma Voom

## dge object
d0 <- DGEList(matrix)

# Design matrix

design <- model.matrix(~ 0 +  metadata$time)

# Limma low expressed genes filtering
keep <- filterByExpr(d0, design)
d <- d0[keep,,keep.lib.sizes=FALSE]

# Scale normalization
d <- calcNormFactors(d)


### Principal component Analisys

###colors
col.cell <- rep(c(rep(c("purple"), 3), rep(c("orange"), 3)),6)
col.cell.time <-rep(c("lightblue", "khaki", "brown", "lightgreen", "violet", "coral"), each = 6)

par(mfrow =c(1,1))
pdf(paste0(wkdir,"/plots/PCA.Condition.pdf"), width=12, height=12)
plotMDS(d, col=col.cell, main="Condition",pch = 19,gene.selection="common")
legend("topright", 
       legend=levels(metadata$condition),
       fill=c("orange", "purple"),
       col= coll.cell)
dev.off()


pdf(paste0(wkdir,"/plots/PCA.Time.pdf"), width=12, height=12)
par(mfrow =c(1,1))
plotMDS(d, col=col.cell.time, main="Time",pch = 19, gene.selection="common")
legend("topright",
       fill=c("lightblue", "khaki", "brown", "lightgreen", "violet", "coral"), 
       legend=levels(metadata$time))

dev.off()


pdf(paste0(wkdir,"/plots/PCA.sample.pdf"), width=12, height=12)
par(mfrow =c(1,1))
plotMDS(d, main="sample", gene.selection="common")

dev.off()


# Contrast 1. Wild type versus Depleted. Signature for Neutrophils

mm.WT.VS.DP <- model.matrix(~ 0 +  metadata$condition)

colnames(mm.WT.VS.DP) <- c("DP", "WT")

fit.DEA.WT.VS.DP <- differential.expression.voom(matrix, metadata, mm.WT.VS.DP) # calling function in functions.R


contr.DEA.WT.VS.DP <- makeContrasts(WT - DP,
                        levels = colnames(coef(fit.DEA.WT.VS.DP)))

tmp.DEA.WT.VS.DP <- contrasts.fit(fit.DEA.WT.VS.DP, contr.DEA.WT.VS.DP)
tmp.DEA.WT.VS.DP.voom <- eBayes(tmp.DEA.WT.VS.DP)
top.table.DEA.WT.VS.DP.voom <- topTable(tmp.DEA.WT.VS.DP.voom, sort.by = "P", n = Inf,adjust.method="fdr")


top.table.DEA.WT.VS.DP.voom<-as.data.frame(top.table.DEA.WT.VS.DP.voom)
top.table.DEA.WT.VS.DP.voom$diffexpressed <- "NO"
top.table.DEA.WT.VS.DP.voom$diffexpressed[top.table.DEA.WT.VS.DP.voom$adj.P.Val < 0.05 & top.table.DEA.WT.VS.DP.voom$logFC > 0] <- "UP"
top.table.DEA.WT.VS.DP.voom$diffexpressed[top.table.DEA.WT.VS.DP.voom$adj.P.Val < 0.05 & top.table.DEA.WT.VS.DP.voom$logFC < 0] <- "DOWN"

head(top.table.DEA.WT.VS.DP.voom, n=20)



# Contrast 1. Wild type versus Depleted. Signature for Neutrophils. Time correction added to linear model.
## Get better fit, I used this one.

mm.WT.VS.DP <- model.matrix(~ 0 + metadata$condition + metadata$time)

colnames(mm.WT.VS.DP) <- c("DP", "WT", "ZT5", "ZT9", "ZT13", "ZT17", "ZT21")

fit.DEA.WT.VS.DP <- differential.expression.voom(matrix, metadata, mm.WT.VS.DP)


contr.DEA.WT.VS.DP <- makeContrasts(WT - DP,
                        levels = colnames(coef(fit.DEA.WT.VS.DP)))

tmp.DEA.WT.VS.DP <- contrasts.fit(fit.DEA.WT.VS.DP, contr.DEA.WT.VS.DP)
tmp.DEA.WT.VS.DP.voom <- eBayes(tmp.DEA.WT.VS.DP)
top.table.DEA.WT.VS.DP.voom.corrected <- topTable(tmp.DEA.WT.VS.DP.voom, sort.by = "P", n = Inf,adjust.method="fdr")


top.table.DEA.WT.VS.DP.voom.corrected<-as.data.frame(top.table.DEA.WT.VS.DP.voom.corrected)
top.table.DEA.WT.VS.DP.voom.corrected$diffexpressed <- "NO"
top.table.DEA.WT.VS.DP.voom.corrected$diffexpressed[top.table.DEA.WT.VS.DP.voom.corrected$adj.P.Val < 0.05 & top.table.DEA.WT.VS.DP.voom.corrected$logFC > 0] <- "UP"
top.table.DEA.WT.VS.DP.voom.corrected$diffexpressed[top.table.DEA.WT.VS.DP.voom.corrected$adj.P.Val < 0.05 & top.table.DEA.WT.VS.DP.voom.corrected$logFC < 0] <- "DOWN"

head(top.table.DEA.WT.VS.DP.voom.corrected, n=20)


# Performance of Voom for this data looks better with Time correction
# We will keep the last result.

write.table(top.table.DEA.WT.VS.DP.voom,paste0(wkdir,"/results/DEG.WT.VS.DP.tsv"),
sep = "\t",
col.names = NA)



# Contrast 2. Night vs Day in WT and Depleted

## Isolating Wild Type

metadata.WT <- metadata %>% filter(grepl('WT', condition))

mm.Night.VS.Day.WT <- model.matrix(~ 0 +  metadata.WT$cycle) # Now we are testing DvsN

colnames(mm.Night.VS.Day.WT) <- c("Day", "Night")

matrix.WT <- matrix[, metadata.WT$sample]

all.equal(metadata.WT$sample, colnames(matrix.WT)) #Chek for the order, should be the same

fit.DEA.Night.VS.Day.WT <- differential.expression.voom(matrix.WT, metadata.WT, mm.Night.VS.Day.WT)

contr.DEA.Night.VS.Day.WT <- makeContrasts(Night - Day,
                                    levels = colnames(coef(fit.DEA.Night.VS.Day.WT)))

tmp.DEA.Night.VS.Day.WT <- contrasts.fit(fit.DEA.Night.VS.Day.WT, contr.DEA.Night.VS.Day.WT)
tmp.DEA.Night.VS.Day.WT.voom <- eBayes(tmp.DEA.Night.VS.Day.WT)
top.table.DEA.Night.VS.Day.WT.voom <- topTable(tmp.DEA.Night.VS.Day.WT.voom, sort.by = "P", n = Inf,adjust.method="fdr")


top.table.DEA.Night.VS.Day.WT.voom <- as.data.frame(top.table.DEA.Night.VS.Day.WT.voom)
top.table.DEA.Night.VS.Day.WT.voom$diffexpressed <- "NO"
top.table.DEA.Night.VS.Day.WT.voom$diffexpressed[top.table.DEA.Night.VS.Day.WT.voom$adj.P.Val < 0.05 & top.table.DEA.Night.VS.Day.WT.voom$logFC > 0] <- "UP"
top.table.DEA.Night.VS.Day.WT.voom$diffexpressed[top.table.DEA.Night.VS.Day.WT.voom$adj.P.Val < 0.05 & top.table.DEA.Night.VS.Day.WT.voom$logFC < 0] <- "DOWN"

head(top.table.DEA.Night.VS.Day.WT.voom, n=20)

write.table(top.table.DEA.Night.VS.Day.WT.voom,paste0(wkdir,"/results/DEG.Night.VS.Day.WT.tsv"),
            sep = "\t",
            col.names = NA)


# Isolating Depleted type

metadata.DP <- metadata %>% filter(grepl('depleted', condition))

mm.Night.VS.Day.DP <- model.matrix(~ 0 +  metadata.DP$cycle)

colnames(mm.Night.VS.Day.DP) <- c("Day", "Night")

matrix.DP <- matrix[, metadata.DP$sample]

all.equal(metadata.DP$sample, colnames(matrix.DP))

fit.DEA.Night.VS.Day.DP <- differential.expression.voom(matrix.DP, metadata.DP, mm.Night.VS.Day.DP)


contr.DEA.Night.VS.Day.DP <- makeContrasts(Night - Day,
                                           levels = colnames(coef(fit.DEA.Night.VS.Day.DP)))

tmp.DEA.Night.VS.Day.DP <- contrasts.fit(fit.DEA.Night.VS.Day.DP, contr.DEA.Night.VS.Day.DP)
tmp.DEA.Night.VS.Day.DP.voom <- eBayes(tmp.DEA.Night.VS.Day.DP)
top.table.DEA.Night.VS.Day.DP.voom <- topTable(tmp.DEA.Night.VS.Day.DP.voom, sort.by = "P", n = Inf,adjust.method="fdr")


top.table.DEA.Night.VS.Day.DP.voom <- as.data.frame(top.table.DEA.Night.VS.Day.DP.voom)
top.table.DEA.Night.VS.Day.DP.voom$diffexpressed <- "NO"
top.table.DEA.Night.VS.Day.DP.voom$diffexpressed[top.table.DEA.Night.VS.Day.DP.voom$adj.P.Val < 0.05 & top.table.DEA.Night.VS.Day.DP.voom$logFC > 0] <- "UP"
top.table.DEA.Night.VS.Day.DP.voom$diffexpressed[top.table.DEA.Night.VS.Day.DP.voom$adj.P.Val < 0.05 & top.table.DEA.Night.VS.Day.DP.voom$logFC < 0] <- "DOWN"

head(top.table.DEA.Night.VS.Day.DP.voom, n=20)

write.table(top.table.DEA.Night.VS.Day.DP.voom,paste0(wkdir,"/results/DEG.Night.VS.Day.DP.tsv"),
            sep = "\t",
            col.names = NA)


# Contrast 3. WT vs DP in night

metadata.WT.VS.DP.Night <- metadata %>% filter(grepl('Night', cycle))

mm.WT.VS.DP.Night  <- model.matrix(~ 0 +  metadata.WT.VS.DP.Night$condition)

colnames(mm.WT.VS.DP.Night) <- c("DP", "WT")

matrix.WT.VS.DP.Night <- matrix[, metadata.WT.VS.DP.Night$sample]

all.equal(metadata.WT.VS.DP.Night$sample, colnames(matrix.WT.VS.DP.Night))

fit.DEA.WT.VS.DP.Night <- differential.expression.voom(matrix.WT.VS.DP.Night, metadata.WT.VS.DP.Night, mm.WT.VS.DP.Night)


contr.DEA.WT.VS.DP.Night <- makeContrasts( WT - DP,
                                          levels = colnames(coef(fit.DEA.WT.VS.DP.Night )))

tmp.DEA.WT.VS.DP.Night  <- contrasts.fit(fit.DEA.WT.VS.DP.Night , contr.DEA.WT.VS.DP.Night )
tmp.DEA.WT.VS.DP.Night.voom <- eBayes(tmp.DEA.WT.VS.DP.Night )
top.table.DEA.WT.VS.DP.Night.voom <- topTable(tmp.DEA.WT.VS.DP.Night.voom, sort.by = "P", n = Inf,adjust.method="fdr")


top.table.DEA.WT.VS.DP.Night.voom <- as.data.frame(top.table.DEA.WT.VS.DP.Night.voom)
top.table.DEA.WT.VS.DP.Night.voom$diffexpressed <- "NO"
top.table.DEA.WT.VS.DP.Night.voom$diffexpressed[top.table.DEA.WT.VS.DP.Night.voom$adj.P.Val < 0.05 & top.table.DEA.WT.VS.DP.Night.voom$logFC > 0] <- "UP"
top.table.DEA.WT.VS.DP.Night.voom$diffexpressed[top.table.DEA.WT.VS.DP.Night.voom$adj.P.Val < 0.05 & top.table.DEA.WT.VS.DP.Night.voom$logFC < 0] <- "DOWN"

head(top.table.DEA.WT.VS.DP.Night.voom, n=20)

write.table(top.table.DEA.WT.VS.DP.Night.voom,paste0(wkdir,"/results/DEG.WT.VS.DP.Night.tsv"),
            sep = "\t",
            col.names = NA)



# Contrast 4. WT vs DP in day

metadata.WT.VS.DP.Day <- metadata %>% filter(grepl('Day', cycle))

mm.WT.VS.DP.Day  <- model.matrix(~ 0 +  metadata.WT.VS.DP.Day$condition)

colnames(mm.WT.VS.DP.Day) <- c("DP", "WT")

matrix.WT.VS.DP.Day <- matrix[, metadata.WT.VS.DP.Day$sample]

all.equal(metadata.WT.VS.DP.Day$sample, colnames(matrix.WT.VS.DP.Day))

fit.DEA.WT.VS.DP.Day <- differential.expression.voom(matrix.WT.VS.DP.Day, metadata.WT.VS.DP.Day, mm.WT.VS.DP.Day)


contr.DEA.WT.VS.DP.Day <- makeContrasts( WT - DP,
                                           levels = colnames(coef(fit.DEA.WT.VS.DP.Day )))

tmp.DEA.WT.VS.DP.Day  <- contrasts.fit(fit.DEA.WT.VS.DP.Day , contr.DEA.WT.VS.DP.Day)
tmp.DEA.WT.VS.DP.Day.voom <- eBayes(tmp.DEA.WT.VS.DP.Day)
top.table.DEA.WT.VS.DP.Day.voom <- topTable(tmp.DEA.WT.VS.DP.Day.voom, sort.by = "P", n = Inf,adjust.method="fdr")


top.table.DEA.WT.VS.DP.Day.voom <- as.data.frame(top.table.DEA.WT.VS.DP.Day.voom)
top.table.DEA.WT.VS.DP.Day.voom$diffexpressed <- "NO"
top.table.DEA.WT.VS.DP.Day.voom$diffexpressed[top.table.DEA.WT.VS.DP.Day.voom$adj.P.Val < 0.05 & top.table.DEA.WT.VS.DP.Day.voom$logFC > 0] <- "UP"
top.table.DEA.WT.VS.DP.Day.voom$diffexpressed[top.table.DEA.WT.VS.DP.Day.voom$adj.P.Val < 0.05 & top.table.DEA.WT.VS.DP.Day.voom$logFC < 0] <- "DOWN"

head(top.table.DEA.WT.VS.DP.Day.voom, n=20)

write.table(top.table.DEA.WT.VS.DP.Day.voom,paste0(wkdir,"/results/DEG.WT.VS.DP.Day.tsv"),
            sep = "\t",
            col.names = NA)



# Venn Diagram

WT.vs.DP <- row.names(top.table.DEA.WT.VS.DP.voom[top.table.DEA.WT.VS.DP.voom$adj.P.Val < 0.05 ,])

NI.vs.DA.WT <- row.names(top.table.DEA.Night.VS.Day.WT.voom[top.table.DEA.Night.VS.Day.WT.voom$adj.P.Val < 0.05 ,])

NI.vs.DA.DP <- row.names(top.table.DEA.Night.VS.Day.DP.voom[top.table.DEA.Night.VS.Day.DP.voom$adj.P.Val < 0.05 ,])

WT.vs.DP.NI <- row.names(top.table.DEA.WT.VS.DP.Night.voom[top.table.DEA.WT.VS.DP.Night.voom$adj.P.Val < 0.05 ,])

WT.vs.DP.DA <- row.names(top.table.DEA.WT.VS.DP.Day.voom[top.table.DEA.WT.VS.DP.Day.voom$adj.P.Val < 0.05 ,])


library(VennDiagram)

# Night vs Day
list_venn_Night<-list(WildType=NI.vs.DA.WT , Depleted=NI.vs.DA.DP)

v0 <- venn.diagram( list_venn_Night, filename=NULL, 
                    fill = c("orange", "gray84"),
                    alpha = 0.50,
                    col = "transparent")


pdf(paste0(wkdir,"/plots/venn.NvD.pdf"), width=12, height=12)
grid.arrange(gTree(children = v0),
             top="Night vs Day. Diffexpressed genes")
dev.off()


### Checking for overlaps
overlaps <- calculate.overlap(list_venn_Night)

core.genes.WT_DP <- overlaps$a3

### Unique genes
depleted.unique <- setdiff(NI.vs.DA.DP, core.genes.WT_DP)
wildType.unique <- setdiff(NI.vs.DA.WT, core.genes.WT_DP)



# WT vs DP
list_venn_Day<-list(Night=WT.vs.DP.NI , Day=WT.vs.DP.DA )

v0 <- venn.diagram( list_venn_Day, filename=NULL, 
                    fill = c("orange", "gray84"),
                    alpha = 0.50,
                    col = "transparent")


pdf(paste0(wkdir,"/plots/venn.WvD.pdf"), width=12, height=12)
grid.arrange(gTree(children = v0),
             top="Wild Type vs Depleted")
dev.off()


overlaps2 <- calculate.overlap(list_venn_Day)

core.genes.WT_DP2 <- overlaps2$a3

depleted.unique2 <- setdiff(WT.vs.DP.NI, core.genes.WT_DP2)
wildType.unique2 <- setdiff(WT.vs.DP.NI, core.genes.WT_DP2)



# Volcano plots

top.table.DEA.WT.VS.DP.voom.corrected$label <- NA
top.table.DEA.WT.VS.DP.voom.corrected$label[top.table.DEA.WT.VS.DP.voom.corrected$diffexpressed != "NO"] <- 
rownames(top.table.DEA.WT.VS.DP.voom.corrected)[top.table.DEA.WT.VS.DP.voom.corrected$diffexpressed != "NO"]



v.WT.VS.DP <-ggplot(data=top.table.DEA.WT.VS.DP.voom.corrected, aes(x=logFC, y=-log10(adj.P.Val), col=diffexpressed, label=label)) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel() +
  scale_color_manual(values=c("dodgerblue4", "dimgrey", "brown1")) + ggtitle("Wild Type Vs Neutrophil depleted. All times") +
  labs(x = "logFC", y= "-log10(adj.pval)")





top.table.DEA.Night.VS.Day.WT.voom$label <- NA
top.table.DEA.Night.VS.Day.WT.voom$label[top.table.DEA.Night.VS.Day.WT.voom$diffexpressed != "NO"] <- rownames(top.table.DEA.Night.VS.Day.WT.voom)[top.table.DEA.Night.VS.Day.WT.voom$diffexpressed != "NO"]


v.Night.VS.Day.WT <- ggplot(data=top.table.DEA.Night.VS.Day.WT.voom, aes(x=logFC, y=-log10(adj.P.Val), col=diffexpressed, label=label)) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel() +
  scale_color_manual(values=c("dodgerblue4", "dimgrey", "brown1")) + ggtitle("Night Vs Day. Wild Type") +
  labs(x = "logFC", y= "-log10(adj.pval)")




top.table.DEA.Night.VS.Day.DP.voom$label <- NA
top.table.DEA.Night.VS.Day.DP.voom$label[top.table.DEA.Night.VS.Day.DP.voom$diffexpressed != "NO"] <- rownames(top.table.DEA.Night.VS.Day.DP.voom)[top.table.DEA.Night.VS.Day.DP.voom$diffexpressed != "NO"]


v.Night.VS.Day.DP <- ggplot(data=top.table.DEA.Night.VS.Day.DP.voom, aes(x=logFC, y=-log10(adj.P.Val), col=diffexpressed, label=label)) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel() +
  scale_color_manual(values=c("dodgerblue4", "dimgrey", "brown1")) + ggtitle("Night Vs Day. Neutrophil depleted") +
  labs(x = "logFC", y= "-log10(adj.pval)")



top.table.DEA.WT.VS.DP.Night.voom$label <- NA
top.table.DEA.WT.VS.DP.Night.voom$label[top.table.DEA.WT.VS.DP.Night.voom$diffexpressed != "NO"] <- rownames(top.table.DEA.WT.VS.DP.Night.voom)[top.table.DEA.WT.VS.DP.Night.voom$diffexpressed != "NO"]


v.WT.VS.DP.Night <- ggplot(data=top.table.DEA.WT.VS.DP.Night.voom, aes(x=logFC, y=-log10(adj.P.Val), col=diffexpressed, label=label)) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel() +
  scale_color_manual(values=c("dodgerblue4", "dimgrey", "brown1")) + ggtitle("Wild Type Vs Depleted. Night") +
  labs(x = "logFC", y= "-log10(adj.pval)")



top.table.DEA.WT.VS.DP.Day.voom$label <- NA
top.table.DEA.WT.VS.DP.Day.voom$label[top.table.DEA.WT.VS.DP.Day.voom$diffexpressed != "NO"] <- rownames(top.table.DEA.WT.VS.DP.Day.voom)[top.table.DEA.WT.VS.DP.Day.voom$diffexpressed != "NO"]


v.WT.VS.DP.Day <- ggplot(data=top.table.DEA.WT.VS.DP.Day.voom, aes(x=logFC, y=-log10(adj.P.Val), col=diffexpressed, label=label)) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel() +
  scale_color_manual(values=c("dodgerblue4", "dimgrey", "brown1")) + ggtitle("Wild Type Vs Depleted. Day") +
  labs(x = "logFC", y= "-log10(adj.pval)")



pdf(paste0(wkdir,"/plots/volcanos.pdf"), width=25, height=10)
plot_grid(v.WT.VS.DP, v.Night.VS.Day.WT, v.Night.VS.Day.DP, v.WT.VS.DP.Night, v.WT.VS.DP.Day)
dev.off()


####################### Enrichment Analysis for unique diffexpressed genes in each group #############


library(R.utils)
R.utils::setOption("clusterProfiler.download.method","wget")

# For Depleted

depleted.unique.GSEA <- top.table.DEA.Night.VS.Day.DP.voom[row.names(top.table.DEA.Night.VS.Day.DP.voom) %in% depleted.unique,] 
depleted.unique.GSEA <- depleted.unique.GSEA[order(abs(depleted.unique.GSEA$logFC), decreasing = TRUE),]

depleted.unique.GSEA <- depleted.unique.GSEA[1]
depleted.unique.GSEA$logFC <- abs(depleted.unique.GSEA$logFC)

depleted.unique.GSEA.symbol <- depleted.unique.GSEA

x <- row.names(depleted.unique.GSEA)

# Cluster profiler works with entrezgene codes. Getting a new annotation

mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
IDs <- biomaRt::getBM(attributes=c("entrezgene_id", "mgi_symbol"), 
                      filters = "mgi_symbol", values= x, mart=mart)

 # Input for enrichment should be an ordered named vector qith entrez codes

depleted.unique.GSEA$mgi_symbol <- row.names(depleted.unique.GSEA)

depleted.unique.GSEA <- merge(depleted.unique.GSEA, IDs, by="mgi_symbol")
depleted.unique.GSEA <- na.omit(depleted.unique.GSEA)
row.names(depleted.unique.GSEA) <- depleted.unique.GSEA$entrezgene_id


depleted.unique.GSEA <- depleted.unique.GSEA[order(abs(depleted.unique.GSEA$logFC), decreasing = TRUE),]


depleted.unique.GSEA <- depleted.unique.GSEA[2]

foldchanges <- depleted.unique.GSEA$logFC
names(foldchanges) <- row.names(depleted.unique.GSEA)


kegg.unique.DP <- enrichKEGG(gene = row.names(depleted.unique.GSEA), 
                       organism = 'mmu',
                       pvalueCutoff = 0.05,
                       minGSSize = 15,
                       use_internal_data=F)

go.unique.DP <- enrichGO(gene = row.names(depleted.unique.GSEA),
                     OrgDb = org.Mm.eg.db,
                     ont = "BP",
                     pAdjustMethod ="BH",
                     pvalueCutoff = 0.05)


pdf(paste0(wkdir,"/plots/enrichKEGG.unique.DP.pdf"), width=12, height=12)
dotplot(kegg.unique.DP, showCategory=30) + ggtitle("Enrichment Analysis KEGG. Night vs Day. Neutrophil Depleted unique diffexpressed genes")
dev.off()


pdf(paste0(wkdir,"/plots/enrichGO.unique.DP.pdf"), width=12, height=12)
dotplot(go.unique.DP, showCategory=30) + ggtitle("Enrichment Analysis GO. Night vs Day. Neutrophil Depleted unique diffexpressed genes")
dev.off()



### For Wild Type


wt.unique.GSEA <- top.table.DEA.Night.VS.Day.WT.voom[row.names(top.table.DEA.Night.VS.Day.WT.voom) %in% wildType.unique,] 
wt.unique.GSEA <- wt.unique.GSEA[order(abs(wt.unique.GSEA$logFC), decreasing = TRUE),]

wt.unique.GSEA <- wt.unique.GSEA[1]
wt.unique.GSEA$logFC <- abs(wt.unique.GSEA$logFC)

wt.unique.GSEA.symbol <- wt.unique.GSEA

x <- row.names(wt.unique.GSEA)

mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
IDs <- biomaRt::getBM(attributes=c("entrezgene_id", "mgi_symbol"), 
                      filters = "mgi_symbol", values= x, mart=mart)


wt.unique.GSEA$mgi_symbol <- row.names(wt.unique.GSEA)

wt.unique.GSEA <- merge(wt.unique.GSEA, IDs, by="mgi_symbol")
wt.unique.GSEA <- na.omit(wt.unique.GSEA)
row.names(wt.unique.GSEA) <- wt.unique.GSEA$entrezgene_id


wt.unique.GSEA <- wt.unique.GSEA[order(abs(wt.unique.GSEA$logFC), decreasing = TRUE),]


wt.unique.GSEA <- wt.unique.GSEA[2]

foldchanges <- wt.unique.GSEA$logFC
names(foldchanges) <- row.names(wt.unique.GSEA)


kegg.unique.WT <- enrichKEGG(gene = row.names(wt.unique.GSEA), 
                       organism = 'mmu',
                       pvalueCutoff = 0.05,
                       minGSSize = 15,
                       use_internal_data=F)

go.unique.WT <- enrichGO(gene = row.names(wt.unique.GSEA),
                     OrgDb = org.Mm.eg.db,
                     ont = "BP",
                     pAdjustMethod ="BH",
                     pvalueCutoff = 0.05)


pdf(paste0(wkdir,"/plots/enrichKEGG.unique.WT.pdf"), width=12, height=12)
dotplot(kegg.unique.WT, showCategory=30) + ggtitle("Enrichment Analysis KEGG. Night vs Day. Wild Type unique diffexpressed genes")
dev.off()


pdf(paste0(wkdir,"/plots/enrichGO.unique.WT.pdf"), width=12, height=12)
dotplot(go.unique.WT, showCategory=30) + ggtitle("Enrichment Analysis GO. Night vs Day. Wild Type unique diffexpressed genes")
dev.off()



######### Gene set Enrichment Analysis ############
## I have to review the used metric. Final NES values are displayed in ABS values


# For depleted

depleted <- top.table.DEA.Night.VS.Day.DP.voom[order(abs(top.table.DEA.Night.VS.Day.DP.voom$logFC), 
decreasing = TRUE),]

depleted <- abs(depleted[1])


# Getting entrez
mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
IDs <- biomaRt::getBM(attributes=c("entrezgene_id", "mgi_symbol"), 
                      filters = "mgi_symbol", values= row.names(depleted), mart=mart)


depleted$mgi_symbol <- row.names(depleted)

depleted <- merge(depleted, IDs, by="mgi_symbol")
depleted <- na.omit(depleted)

which(duplicated(depleted$entrezgene_id))

depleted <-depleted[!duplicated(depleted$entrezgene_id), ]


row.names(depleted) <- depleted$entrezgene_id


depleted <- depleted[order(depleted$logFC, decreasing = TRUE),] # maybe a mixed metric here??. Check this


depleted <- depleted[2]

foldchanges <- depleted$logFC
names(foldchanges) <- row.names(depleted)


gseGO.NvD.Depleted <- gseGO(geneList = foldchanges, 
                     ont = "BP", 
                     OrgDb = "org.Mm.eg.db", 
                     minGSSize = 15, 
                     maxGSSize = 500, 
                     eps = 1e-10, 
                     nPermSimple = 10000, 
                     seed = TRUE,
                     pAdjustMethod ="BH",
                     pvalueCutoff = 0.05)

  gseGO.NvD.Depleted@result %>% 
  arrange(pvalue) %>% 
  head(8)


pdf(paste0(wkdir,"/plots/gseGO.NvD.Depleted.pdf"), width=12, height=12)
dotplot(gseGO.NvD.Depleted, showCategory=30, x="NES") + ggtitle("GSEA GO. Nigh vs Day. Neutrophil Depleted")
dev.off()

# For WT

wild.type <- top.table.DEA.Night.VS.Day.WT.voom[order(abs(top.table.DEA.Night.VS.Day.WT.voom$logFC), 
decreasing = TRUE),]

wild.type <- abs(wild.type[1])


mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
IDs <- biomaRt::getBM(attributes=c("entrezgene_id", "mgi_symbol"), 
                      filters = "mgi_symbol", values= row.names(wild.type), mart=mart)


wild.type$mgi_symbol <- row.names(wild.type)

wild.type <- merge(wild.type, IDs, by="mgi_symbol")
wild.type <- na.omit(wild.type)

which(duplicated(wild.type$entrezgene_id))

wild.type <-wild.type[!duplicated(wild.type$entrezgene_id), ]


row.names(wild.type) <- wild.type$entrezgene_id


wild.type <- wild.type[order(wild.type$logFC, decreasing = TRUE),]


wild.type <- wild.type[2]

foldchanges.WT <- wild.type$logFC
names(foldchanges.WT) <- row.names(wild.type)


gseGO.NvD.WT <- gseGO(geneList = foldchanges.WT, 
                     ont = "BP", 
                     OrgDb = "org.Mm.eg.db", 
                     minGSSize = 15, 
                     maxGSSize = 500, 
                     eps = 1e-10, 
                     nPermSimple = 10000, 
                     seed = TRUE)

  gseGO.NvD.WT@result %>% 
  arrange(pvalue) %>% 
  head(8)


pdf(paste0(wkdir,"/plots/gseGO.NvD.WT.pdf"), width=12, height=12)
dotplot(gseGO.NvD.WT, showCategory=30, x="NES") + ggtitle("GSEA Gene GO. Nigh vs Day. Wild Type")
dev.off()


### GSEA for WT vs Depleted all times. Corrected by time



WT.DT.C <- top.table.DEA.WT.VS.DP.voom.corrected[order(abs(top.table.DEA.WT.VS.DP.voom.corrected$logFC), 
decreasing = TRUE),]

WT.DT.C <- abs(WT.DT.C[1])

mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
IDs <- biomaRt::getBM(attributes=c("entrezgene_id", "mgi_symbol"), 
                      filters = "mgi_symbol", values= row.names(WT.DT.C),
                       mart=mart)


WT.DT.C$mgi_symbol <- row.names(WT.DT.C)

WT.DT.C <- merge(WT.DT.C, IDs, by="mgi_symbol")
WT.DT.C <- na.omit(WT.DT.C)

which(duplicated(WT.DT.C$entrezgene_id))

WT.DT.C <-WT.DT.C[!duplicated(WT.DT.C$entrezgene_id), ]


row.names(WT.DT.C) <- WT.DT.C$entrezgene_id


WT.DT.C <- WT.DT.C[order(WT.DT.C$logFC, decreasing = TRUE),]


WT.DT.C <- WT.DT.C[2]

foldchanges <- WT.DT.C$logFC
names(foldchanges) <- row.names(WT.DT.C)


gseGO.ALL.WT.DP.C <- gseGO(geneList = foldchanges, 
                     ont = "BP", 
                     OrgDb = "org.Mm.eg.db", 
                     minGSSize = 150, 
                     maxGSSize = 500, 
                     eps = 1e-10, 
                     nPermSimple = 10000, 
                     seed = TRUE,
                     pAdjustMethod ="BH",
                     pvalueCutoff = 0.05)

  gseGO.ALL.WT.DP.C@result %>% 
  arrange(pvalue) %>% 
  head(8)


pdf(paste0(wkdir,"/plots/gseGO.ALL.WT.DP.C.pdf"), width=12, height=12)
dotplot(gseGO.ALL.WT.DP.C, showCategory=30, x="NES") + ggtitle("GSEA GO. Wild Type vs Neutrophil Depleted. All times. Corrected")
dev.off()




### GSEA for WT vs Depleted Night


WT.DT.N <- top.table.DEA.WT.VS.DP.Night.voom[order(abs(top.table.DEA.WT.VS.DP.Night.voom$logFC), 
decreasing = TRUE),]

WT.DT.N <- abs(WT.DT.N[1])

mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
IDs <- biomaRt::getBM(attributes=c("entrezgene_id", "mgi_symbol"), 
                      filters = "mgi_symbol", values= row.names(WT.DT.N),
                       mart=mart)


WT.DT.N$mgi_symbol <- row.names(WT.DT.N)

WT.DT.N <- merge(WT.DT.N, IDs, by="mgi_symbol")
WT.DT.N <- na.omit(WT.DT.N)

which(duplicated(WT.DT.N$entrezgene_id))

WT.DT.N <-WT.DT.N[!duplicated(WT.DT.N$entrezgene_id), ]


row.names(WT.DT.N) <- WT.DT.N$entrezgene_id


WT.DT.N <- WT.DT.N[order(WT.DT.N$logFC, decreasing = TRUE),]


WT.DT.N <- WT.DT.N[2]

foldchanges <- WT.DT.N$logFC
names(foldchanges) <- row.names(WT.DT.N)


gseGO.WT.DT.N <- gseGO(geneList = foldchanges, 
                     ont = "BP", 
                     OrgDb = "org.Mm.eg.db", 
                     minGSSize = 150, 
                     maxGSSize = 500, 
                     eps = 1e-10, 
                     nPermSimple = 10000, 
                     seed = TRUE,
                     pAdjustMethod ="BH",
                     pvalueCutoff = 0.05)

  gseGO.WT.DT.N@result %>% 
  arrange(pvalue) %>% 
  head(8)

gseGO.WT.DT.N <- setReadable(gseGO.WT.DT.N, OrgDb = org.Mm.eg.db, keyType="ENTREZID")


pdf(paste0(wkdir,"/plots/gseGO.WT.DT.N.pdf"), width=12, height=12)
dotplot(gseGO.WT.DT.N, showCategory=30, x="NES") + ggtitle("GSEA GO. Wild Type vs Neutrophil Depleted. Night")
dev.off()



gseGO.WT.DT.N2 <- gseGO(geneList = foldchanges, 
                     ont = "BP", 
                     OrgDb = "org.Mm.eg.db", 
                     minGSSize = 15, 
                     maxGSSize = 500, 
                     eps = 1e-10, 
                     nPermSimple = 10000, 
                     seed = TRUE,
                     pAdjustMethod ="BH",
                     pvalueCutoff = 1)

  gseGO.WT.DT.N2@result %>% 
  arrange(pvalue) %>% 
  head(8)


gseGO.WT.DT.N2 <- setReadable(gseGO.WT.DT.N2, OrgDb = org.Mm.eg.db, keyType="ENTREZID")


### GSEA for WT vs Depleted Day


WT.DT.D <- top.table.DEA.WT.VS.DP.Day.voom[order(abs(top.table.DEA.WT.VS.DP.Day.voom$logFC), 
decreasing = TRUE),]

WT.DT.D <- abs(WT.DT.D[1])

mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
IDs <- biomaRt::getBM(attributes=c("entrezgene_id", "mgi_symbol"), 
                      filters = "mgi_symbol", values= row.names(WT.DT.D),
                       mart=mart)


WT.DT.D$mgi_symbol <- row.names(WT.DT.D)

WT.DT.D <- merge(WT.DT.D, IDs, by="mgi_symbol")
WT.DT.D <- na.omit(WT.DT.D)

which(duplicated(WT.DT.N$entrezgene_id))

WT.DT.D <-WT.DT.D[!duplicated(WT.DT.D$entrezgene_id), ]


row.names(WT.DT.D) <- WT.DT.D$entrezgene_id


WT.DT.D <- WT.DT.D[order(WT.DT.D$logFC, decreasing = TRUE),]


WT.DT.D <- WT.DT.D[2]

foldchanges <- WT.DT.D$logFC
names(foldchanges) <- row.names(WT.DT.D)


gseGO.WT.DT.D <- gseGO(geneList = foldchanges, 
                     ont = "BP", 
                     OrgDb = "org.Mm.eg.db", 
                     minGSSize = 150, 
                     maxGSSize = 500, 
                     eps = 1e-10, 
                     nPermSimple = 10000, 
                     seed = TRUE,
                     pAdjustMethod ="BH",
                     pvalueCutoff = 0.05)

  gseGO.WT.DT.D@result %>% 
  arrange(pvalue) %>% 
  head(8)
gseGO.WT.DT.D <- setReadable(gseGO.WT.DT.D, OrgDb = org.Mm.eg.db, keyType="ENTREZID")



pdf(paste0(wkdir,"/plots/gseGO.WT.DT.D.pdf"), width=12, height=12)
dotplot(gseGO.WT.DT.D, showCategory=30, x="NES") + ggtitle("GSEA GO. Wild Type vs Neutrophil Depleted. Day")
dev.off()



gseGO.WT.DT.D2 <- gseGO(geneList = foldchanges, 
                     ont = "BP", 
                     OrgDb = "org.Mm.eg.db", 
                     minGSSize = 150, 
                     maxGSSize = 500, 
                     eps = 1e-10, 
                     nPermSimple = 10000, 
                     seed = TRUE,
                     pAdjustMethod ="BH",
                     pvalueCutoff = 1)

  gseGO.WT.DT.D2@result %>% 
  arrange(pvalue) %>% 
  head(8)

gseGO.WT.DT.D2 <- setReadable(gseGO.WT.DT.D2, OrgDb = org.Mm.eg.db, keyType="ENTREZID")


######## Plotting things


# Night vs Day. Enrichment


# GO Night vs Day WT vs DP

gseGO.WT.DT.N2.DF <- as.data.frame(gseGO.WT.DT.N2)
gseGO.WT.DT.D2.DF <- as.data.frame(gseGO.WT.DT.D)


gseGO.WT.DT.N2.DF$group <- "Night"
gseGO.WT.DT.D2.DF$group <- "Day"

head(gseGO.WT.DT.N2.DF)
head(gseGO.WT.DT.D2.DF)




#### Pathway subsetting

pathways.Day <- c("GO:0072593", "GO:0050727", "GO:0045637", "GO:0007159", "GO:0050900",
"GO:0032944", "GO:0050870", "GO:0046651", "GO:0002274", "GO:0032943", "GO:1903706", "GO:0002694", "GO:0070661", "GO:0042113")


pathways.Night <- c("GO:0030198", "GO:0043062", "GO:0045229", "GO:0045765", "GO:1901342", 
"GO:0048514", "GO:0034764", "GO:0031589", "GO:0009617", "GO:0001525", "GO:0050673")


pathways <- c(pathways.Day, pathways.Night)


### Joining

gseGO.WT.DT.ND <- rbind (gseGO.WT.DT.N2.DF,gseGO.WT.DT.D2.DF)

gseGO.WT.DT.ND <- gseGO.WT.DT.ND[gseGO.WT.DT.ND$ID %in% pathways,]
 

 write.table(gseGO.WT.DT.ND,paste0(wkdir,"/results/gseGO.WT.ND.DT.tsv"),
sep = "\t",
col.names = NA)



library(ggh4x)

strip <- strip_themed(background_x = elem_list_rect(
       fill = c("gold", "dodgerblue4")))

pdf(paste0(wkdir,"/plots/gseGO.WT.DT.ND.JOIN.pdf"), width=12, height=12)


ggplot(gseGO.WT.DT.ND, aes(y=NES, x=reorder(Description,NES), fill=-p.adjust)) +
  geom_bar(stat='identity', position=position_dodge()) +
  theme(axis.text.x = element_text(size=rel(1.15)),axis.title = element_text(size=rel(1.15))) + facet_wrap2(~group, strip=strip) +   theme(
                                                                                     axis.text.x=element_blank(),
                                                                                     axis.ticks.x=element_blank())  + coord_flip() +
labs(title = "Enriched Pathways in Wild Type vs Neutophil Depleted contrast",
              subtitle = "Night vs Day differences", x='Pathway', y='Normalized Enrichment Score') +
scale_fill_gradientn( colours = c("blue", "lightcyan1", "moccasin", "red3"))  + theme_bw() +
theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA)) +
        theme(
      strip.text.x = element_text(
        size = 14, face = "bold"
        ))
        

dev.off()


dim(gseGO.WT.DT.N)

dim(gseGO.WT.DT.D)



##### Scaterplots


genes <- c("Csf3r", "S100a9", "Mki67", "Arntl", "Cry1", "Rorc", "Nr1d1", "Nr1d2", "Clock", "Dbp", "Top2a")



metadata$time <- paste0("ZT", metadata$time)

matrix.circadian <- matrix[genes,]

matrix.circadian<- matrix.circadian[, metadata$sample]

matrix.circadian <- t(matrix.circadian)

matrix.circadian <- as.data.frame(matrix.circadian)

matrix.circadian.1.3w <- colMeans(matrix.circadian[1:3,]) 
matrix.circadian.7.9w <- colMeans(matrix.circadian[7:9,]) 
matrix.circadian.13.15w <- colMeans(matrix.circadian[13:15,]) 
matrix.circadian.19.21w <- colMeans(matrix.circadian[19:21,])
matrix.circadian.25.27w <- colMeans(matrix.circadian[25:27,])
matrix.circadian.31.33w <- colMeans(matrix.circadian[31:33,])

matrix.circadian.4.6w <- colMeans(matrix.circadian[4:6,])
matrix.circadian.10.12w <- colMeans(matrix.circadian[10:12,])
matrix.circadian.16.18w <- colMeans(matrix.circadian[16:18,])
matrix.circadian.22.24w <- colMeans(matrix.circadian[22:24,])
matrix.circadian.28.30w <- colMeans(matrix.circadian[28:30,])
matrix.circadian.34.36w <- colMeans(matrix.circadian[34:36,])


data.circadian <- data.frame(matrix.circadian.1.3w, matrix.circadian.4.6w,
 matrix.circadian.7.9w, matrix.circadian.10.12w, 
 matrix.circadian.13.15w, matrix.circadian.16.18w,
 matrix.circadian.19.21w, matrix.circadian.22.24w ,
 matrix.circadian.25.27w, matrix.circadian.28.30w,
 matrix.circadian.31.33w, matrix.circadian.34.36w)    


data.circadian <- as.data.frame(t(data.circadian))

data.circadian$time <- c("ZT1","ZT1", "ZT5","ZT5", "ZT9","ZT9",
"ZT13","ZT13", "ZT17","ZT17", "ZT21","ZT21")


data.circadian$condition <- c("Wild Type", "Neutrophil Depleted", "Wild Type", "Neutrophil Depleted",
"Wild Type", "Neutrophil Depleted","Wild Type", "Neutrophil Depleted",
"Wild Type", "Neutrophil Depleted","Wild Type", "Neutrophil Depleted")

data.circadian$time <- factor(data.circadian$time)


library(ggplot2)
library(dplyr)


# Plot
g1 <-data.circadian %>%
  ggplot( aes(x=time, y=Csf3r, color=condition, group=condition)) +
    geom_line() +
    geom_point() +
labs(title="Csf3r Circadian behavior", y= "Average Normalized Expression") +
theme_classic() + scale_color_manual(values=c("orange","purple")) + scale_x_discrete(limit =c("ZT1","ZT5","ZT9",
"ZT13", "ZT17", "ZT21"))


# Plot
g2 <-data.circadian %>%
  ggplot( aes(x=time, y=S100a9, color=condition, group=condition)) +
    geom_line() +
    geom_point() +
labs(title="S100a9 Circadian behavior", y= "Average Normalized Expression") +
theme_classic() + scale_color_manual(values=c("orange","purple")) + scale_x_discrete(limit =c("ZT1","ZT5","ZT9",
"ZT13", "ZT17", "ZT21"))





# Plot
g3 <-data.circadian %>%
  ggplot( aes(x=time, y=Mki67, color=condition, group=condition)) +
    geom_line() +
    geom_point() +
labs(title="Mki67 Circadian behavior", y= "Average Normalized Expression") +
theme_classic() + scale_color_manual(values=c("orange","purple")) + scale_x_discrete(limit =c("ZT1","ZT5","ZT9",
"ZT13", "ZT17", "ZT21"))

# Plot
g4 <-data.circadian %>%
  ggplot( aes(x=time, y=Arntl, color=condition, group=condition)) +
    geom_line() +
    geom_point() +
labs(title="Arntl Circadian behavior", y= "Average Normalized Expression") +
theme_classic() + scale_color_manual(values=c("orange","purple")) + scale_x_discrete(limit =c("ZT1","ZT5","ZT9",
"ZT13", "ZT17", "ZT21"))




# Plot
g5 <-data.circadian %>%
  ggplot( aes(x=time, y=Cry1, color=condition, group=condition)) +
    geom_line() +
    geom_point() +
labs(title="Cry1 Circadian behavior", y= "Average Normalized Expression") +
theme_classic() + scale_color_manual(values=c("orange","purple")) + scale_x_discrete(limit =c("ZT1","ZT5","ZT9",
"ZT13", "ZT17", "ZT21"))


# Plot
g6 <-data.circadian %>%
  ggplot( aes(x=time, y=Rorc, color=condition, group=condition)) +
    geom_line() +
    geom_point() +
labs(title="Rorc Circadian behavior", y= "Average Normalized Expression") +
theme_classic() + scale_color_manual(values=c("orange","purple")) + scale_x_discrete(limit =c("ZT1","ZT5","ZT9",
"ZT13", "ZT17", "ZT21"))





# Plot
g7 <-data.circadian %>%
  ggplot( aes(x=time, y=Nr1d1, color=condition, group=condition)) +
    geom_line() +
    geom_point() +
labs(title="Nr1d1 Circadian behavior", y= "Average Normalized Expression") +
theme_classic() + scale_color_manual(values=c("orange","purple")) + scale_x_discrete(limit =c("ZT1","ZT5","ZT9",
"ZT13", "ZT17", "ZT21"))



# Plot
g8 <-data.circadian %>%
  ggplot( aes(x=time, y=Nr1d2, color=condition, group=condition)) +
    geom_line() +
    geom_point() +
labs(title="Nr1d2 Circadian behavior", y= "Average Normalized Expression") +
theme_classic() + scale_color_manual(values=c("orange","purple")) + scale_x_discrete(limit =c("ZT1","ZT5","ZT9",
"ZT13", "ZT17", "ZT21"))


# Plot
g9 <-data.circadian %>%
  ggplot( aes(x=time, y=Clock, color=condition, group=condition)) +
    geom_line() +
    geom_point() +
labs(title="Clock Circadian behavior", y= "Average Normalized Expression") +
theme_classic() + scale_color_manual(values=c("orange","purple")) + scale_x_discrete(limit =c("ZT1","ZT5","ZT9",
"ZT13", "ZT17", "ZT21"))



# Plot
g10 <-data.circadian %>%
  ggplot( aes(x=time, y=Dbp, color=condition, group=condition)) +
    geom_line() +
    geom_point() +
labs(title="Dbp Circadian behavior", y= "Average Normalized Expression") +
theme_classic() + scale_color_manual(values=c("orange","purple")) + scale_x_discrete(limit =c("ZT1","ZT5","ZT9",
"ZT13", "ZT17", "ZT21"))


# Plot
g10 <-data.circadian %>%
  ggplot( aes(x=time, y=Top2a, color=condition, group=condition)) +
    geom_line() +
    geom_point() +
labs(title="Top2a Circadian behavior", y= "Average Normalized Expression") +
theme_classic() + scale_color_manual(values=c("orange","purple")) + scale_x_discrete(limit =c("ZT1","ZT5","ZT9",
"ZT13", "ZT17", "ZT21"))



pdf(paste0(wkdir,"/plots/genes.circadian.pdf"), width=25, height=10)
plot_grid(g1,g2,g3,g4,g5,g6,g7,g8,g9,g10, ncol=4)
dev.off()




#### WT matrix for CircaN algorithm


metadata.WT <- metadata[metadata$condition == "WT",]

matrix.WT <- matrix[,colnames(matrix) %in% metadata.WT$sample]


write.table(metadata.WT,paste0(wkdir,"/results/metadata.Neutros.WT.tsv"),
sep = "\t",
col.names = NA)

write.table(matrix.WT,paste0(wkdir,"/results/matrix.Neutros.WT.tsv"),
sep = "\t",
col.names = NA)
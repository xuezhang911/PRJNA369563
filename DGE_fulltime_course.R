# load data 
data <- read.delim("/Users/xung0001/Desktop/drug.matrix",header = T)
sample_info <- read.delim("/Users/xung0001/Desktop/sample.txt",header = T)
# name the column
names(sample_info) <- c("GEO","sample","time")

sample_info$sample <- factor(c(rep("negative",3),rep("positive",3),rep("mock",15),rep("TGFb",15),rep("TGFb_retinoic",15),rep("TGFb_retinoic_rapamyin",15),
                               rep("TGFb_butyrate",15)),levels = c("negative", "positive","mock","TGFb","TGFb_retinoic","TGFb_retinoic_rapamyin","TGFb_butyrate"))

library(dplyr) # this allow us to use pipe function
sample_info<- sample_info[,-1,drop=FALSE]
sample_info %>% View() 
rownames(data) <- data$X
data<- data[,-1,drop=FALSE]
data1 <- data[,which(names(data) %in%row.names(sample_info))]
metadata <- sample_info
metadata$time <- as.factor(metadata$time)
levels(metadata$time)
# change the order of levels
metadata$time <- factor(metadata$time, levels = c("0h","2h","6h","24h","48h","144h") )
metadata$time
# round data 
countdata <- data1
countdata <- round(countdata,digits=0)
# make sure the colnames of countdata is insistent with metadata
countdata <- countdata[match(row.names(metadata),colnames(countdata))]
library(DESeq2)
row.names(countdata) <- substring(row.names(countdata),first = 1,last = 15)
#load count matrix 
treat <- factor(paste(metadata$sample,metadata$time,sep="."))
design <- model.matrix(~0+treat)
colnames(design) <- levels(treat)

dds <- DESeqDataSetFromMatrix(countData = countdata,
                              colData = metadata,
                              design = design)
dds <- DESeq(dds)
## step2 sample-level quality control 
#check if there are some genes not expressed 
h <- counts(dds)
sel <-rowSums(counts(dds)) == 0
sprintf("%s%% percent (%s) of %s genes are not expressed",
        round(sum(sel) * 100/ nrow(h),digits=1),sum(sel), nrow(h))
"40.6% percent (23447) of 57773 genes are not expressed"
# check the sequencing depth
library(tidyverse)
dat <- tibble(x=colnames(h),y=colSums(h))
#gglot2 visualize the sequencing depth for all the dataset Rplot1
x11(width = 10,height = 8)
par(mar=c(5,5,4,1)+.4)
p <- ggplot(dat,aes(x,y)) + 
  geom_col() + 
  scale_y_continuous(name="reads") +
  facet_grid(~ factor(x), scales = "free") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90,size=8,face="bold"),
        axis.title.x=element_blank())
p
pdf("time_course_sequencing_depth.pdf",width =10 ,height =10 )
plot(p)
dev.off()
p=NULL

# filter the dataset(30conditions)
dds <- dds[ rowSums(counts(dds)) > 30, ]
nrow(dds)
#normalization
# logtransformed normalization _method1
rld <- rlog(dds)
rld_mat <- assay(rld)
# vst normalization
vsd <- DESeq2::vst(dds)
#normalized counts 
vst <- assay(vsd)
vst %>%  view()
#normalization method view 
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)
normalized_counts %>%  view()%>% suppressMessages()
par( mfrow = c( 1, 3 ) )
plot(vst[rowSums(vst)>1,][,1:3], pch=10,cex=0.3)
plot(rld_mat[rowSums(rld_mat)>1,][,1:3], pch=10,cex=0.3)
plot(log10(counts(dds))[rowSums(counts(dds))>0,][,1:3], pch=10,cex=0.3)
# clustering
# heatmap
# three ways vst and rlog normalized_counts
vsd_cor <- cor(vst)  
rld_cor <- cor(rld_mat)  
j <- cor(counts(dds))
df <- as.data.frame(colData(dds)[, c("sample", "time")])
#plot with heatmap 
library(pheatmap)
pheatmap(rld_cor, cluster_rows=T, show_rownames=F, cluster_cols=T, annotation_col=df, fontsize=6)
pheatmap(vsd_cor, cluster_rows=T, show_rownames=F, cluster_cols=T, annotation_col=df, fontsize=6)
pheatmap(j, cluster_rows=T, show_rownames=F, cluster_cols=T, annotation_col=df, fontsize=6)
# PCA using r-log transformed value
rld_mat <- assay(rld)
pca <- prcomp(t(rld_mat))
summary(pca)
# we know PC1: 38%. PC2: 11%
p <- ggplot(df) + geom_point(aes(x=PC1, y=PC2, color = sample, shape=time))+ 
  xlab(paste0("PC1: ","38%", "variance")) +  ylab(paste0("PC2: ", "11%", "variance"))
p

# we can also discover other PCAs
# Create data frame with metadata and PC3 and PC4 values for input to ggplot
df <- cbind(metadata, pca$x)
ggplot(df) + geom_point(aes(x=PC3, y=PC4, color = sample,shape=time))+ 
  xlab(paste0("PC1: ","7%", "variance")) +  ylab(paste0("PC2: ", "5%", "variance"))
``
#Step3 before DGE analysis
## we check if the samples have big batch effect with sva package
library(sva)
treat <- factor(paste(metadata$sample,metadata$time,sep="."))
design <- model.matrix(~0+treat)
colnames(design) <- levels(treat)

dds <- DESeqDataSetFromMatrix(countData = countdata,
                              colData = metadata,
                              design = design)
dds <- DESeq(dds)
dat <- counts(dds, normalized=TRUE)
idx <- rowMeans(dat) > 1
dat <- dat[idx,]
mod <- model.matrix(~sample, colData(dds))
mod0 <- model.matrix(~ 1, colData(dds))
svseq <- svaseq(dat, mod, mod0, n.sv=2)
#visualize surrogate variables over samples 
par(mfrow=c(2,1),mar=c(3,5,3,1))
stripchart(svseq$sv[,1] ~ dds$sample,vertical=TRUE,main="SV1")
abline(h=0)
stripchart(svseq$sv[,2] ~ dds$sample,vertical=TRUE,main="SV2")
abline(h=0)
#it seems the batch effect is very small
## Step4 DGE analysis
treat <- factor(paste(metadata$sample,metadata$time,sep="."))
design <- model.matrix(~0+treat)
colnames(design) <- levels(treat)

dds <- DESeqDataSetFromMatrix(countData = countdata,
                              colData = metadata,
                              design = design)
dds <- DESeq(dds)

#Differentiatl expression(GLM model for multiple factor)
# Compute contrasts
library(limma)
# correlation 
corfit <- duplicateCorrelation(dat,design)
summary(corfit)
corfit$consensus.correlation
# fitting models
fit <- lmFit(dat,design,correlation = corfit$consensus.correlation)
# set up comparation groups 
contrastMatrix <- makeContrasts(negative.0h,positive.0h,
                             TGFb_vs_Mock_2h= TGFb.2h-mock.2h,
                             TGFb_butyrate_vs_Mock_2h=TGFb_butyrate.2h-mock.2h,
                             TGFb_retinoic_vs_Mock_2h=TGFb_retinoic.2h-mock.2h,
                             TGFb_retinoic_rapamyin_vs_Mock_2h=TGFb_retinoic_rapamyin.2h-mock.2h,
                             TGFb_vs_Mock_6h=TGFb.6h-mock.6h,
                             TGFb_butyrate_vs_Mock_6h=TGFb_butyrate.6h-mock.6h,
                             TGFb_retinoic_vs_Mock_6h= TGFb_retinoic.6h-mock.6h,
                             TGFb_retinoic_rapamyin_vs_Mock_6h=TGFb_retinoic_rapamyin.6h-mock.6h,
                             TGFb_vs_Mock_24h=TGFb.24h-mock.24h,
                             TGFb_butyrate_vs_Mock_24h=TGFb_butyrate.24h-mock.24h,
                             GFb_retinoic_vs_Mock_24h=TGFb_retinoic.24h-mock.24h,
                             TGFb_retinoic_rapamyin_vs_Mock_24h=TGFb_retinoic_rapamyin.24h-mock.24h,
                             TGFb_vs_Mock_48h=TGFb.48h-mock.48h,
                             TGFb_butyrate_vs_Mock_48h=TGFb_butyrate.48h-mock.48h,
                             GFb_retinoic_vs_Mock_48h=TGFb_retinoic.48h-mock.48h,
                             TGFb_retinoic_rapamyin_vs_Mock_48h=TGFb_retinoic_rapamyin.48h-mock.48h,
                             TGFb_vs_Mock_144h=TGFb.144h-mock.144h,
                             TGFb_butyrate_vs_Mock_144h=TGFb_butyrate.144h-mock.144h,
                             GFb_retinoic_vs_Mock_144h=TGFb_retinoic.144h-mock.144h,
                             TGFb_retinoic_rapamyin_vs_Mock_144h=TGFb_retinoic_rapamyin.144h-mock.144h,
                                levels = colnames(design))
# DGE analysis 
fit2 <- contrasts.fit(fit,contrastMatrix)
fit2 <- eBayes(fit2)
res <- topTable(fit2,adjust="BH",sort.by = "F",n=Inf,p.value=0.05)

# take one example :
res_TGFb_vs_Mock_6h=topTable(fit2,coef="TGFb_vs_Mock_6h",n=Inf,p.value=0.05)
# extract annotation
# annotating and exporting results
library("AnnotationDbi")
library("org.Hs.eg.db")
columns(org.Hs.eg.db)
res$symbol <- mapIds(org.Hs.eg.db,
                     keys=row.names(res),
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")
res$entrez <- mapIds(org.Hs.eg.db,
                     keys=row.names(res),
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first")

## check FOXP3 gene expression over time 
# create data frame
x <- na.omit(res[res$symbol=="FOXP3",])
x <- x[1:22]
x <- data.frame(t(log(t(x))))
x <- round(x,1)
x <- data.frame(x)
FOXP3 <-t(x)
TGFb <- c(NA,x$TGFb_vs_Mock_2h,x$TGFb_vs_Mock_6h,x$TGFb_vs_Mock_24h,x$TGFb_vs_Mock_48h,x$TGFb_vs_Mock_144h)
TGFb_b_butyrate <- c(NA,x$TGFb_butyrate_vs_Mock_2h,x$TGFb_butyrate_vs_Mock_6h,x$TGFb_butyrate_vs_Mock_24h,x$TGFb_butyrate_vs_Mock_48h,x$TGFb_butyrate_vs_Mock_144h)
TGFb_ATRAX <- c(NA,x$TGFb_retinoic_vs_Mock_2h,x$TGFb_retinoic_vs_Mock_6h,x$GFb_retinoic_vs_Mock_24h,x$GFb_retinoic_vs_Mock_48h,x$GFb_retinoic_vs_Mock_144h)
TGFb_ATRAX_Rapa <- c(NA,x$TGFb_retinoic_rapamyin_vs_Mock_2h,x$TGFb_retinoic_rapamyin_vs_Mock_6h,x$TGFb_retinoic_rapamyin_vs_Mock_24h,x$TGFb_retinoic_rapamyin_vs_Mock_48h,x$TGFb_retinoic_rapamyin_vs_Mock_144h)
nTreg <- c(x$positive.0h,rep(NA,5))
unstim <- c(x$negative.0h,rep(NA,5))
df <- c()
df$sample <- factor(c(rep("unstim",6),rep("nTreg",6),rep("TGFb",6),rep("TGFb_ATRAX",6),rep("TGFb_ATRAX_Rapa",6)),levels = c("unstim","nTreg","TGFb","TGFb_ATRAX","TGFb_ATRAX_Rapa"))
df$time <- factor(c(rep(c("0h","2h","6h","24h","48h","6d"),5)), levels = c("0h","2h","6h","24h","48h","6d") )
df$FOXP3 <- c(unstim,nTreg,TGFb,TGFb_ATRAX,TGFb_ATRAX_Rapa)              
                    df <- as.data.frame(df)
                    df$FOXP3 <- as.numeric(df$FOXP3)
    str(df$FOXP3)
  # visualization                          
                    p <- ggplot(df,aes(x=time, y=FOXP3, color = sample, shape=time)) + geom_point(na.rm=T,size=3)+geom_line(aes(group=sample),na.rm=T)   
                    f<- p + xlab("time") + ylab("Normalized coutnts")+ggtitle("FOXP3 mRNA")+theme_classic()+theme(
                      plot.title = element_text(size=14, face="bold.italic",hjust = 0.4),
                      axis.text.x = element_text(size = 12,color="black",face = "bold"),
                      axis.text.y = element_text(size = 12,color="black",face="bold"),
                      axis.title.x = element_text(color="black", size=14, face="bold"),
                      axis.title.y = element_text(color="black", size=14, face="bold")
                    )+scale_y_continuous(limits = c(0, 10))
f

# volcano to visualize the differential expressed genes under one comparation
res_TGFb_vs_Mock_6h=topTable(fit2,coef="TGFb_vs_Mock_6h",n=Inf,p.value=0.05)
resOrdered <- res_TGFb_vs_Mock_6h
resOrdered$symbol <- mapIds(org.Hs.eg.db,
                            keys=row.names(resOrdered),
                            column="SYMBOL",
                            keytype="ENSEMBL",
                            multiVals="first")
resOrdered$change <- as.factor(
  ifelse(
    resOrdered$P.Value<0.01 & abs(resOrdered$logFC)>0,
    ifelse(resOrdered$logFC>1, "Up", "Down"),
    "NoDiff"))

# exporting results
resOrdered <- resOrdered[order(resOrdered$P.Value,decreasing = F),]
resOrdered <- resOrdered[order(abs(resOrdered$logFC),decreasing = T),]
resOrdered <- as.data.frame(resOrdered)[1:100,]
write.csv(resOrdered, file="DEG_TGFb_mock_6h.csv")
# volcano
# Gernerating vocalno plot with ggplot2

library(ggplot2)

valcano <- ggplot(data=resOrdered, aes(x=logFC, y=-log10(P.Value), color=change)) + 
  geom_point(alpha=0.8, size=1) + 
  theme_bw(base_size=15) + 
  theme(
    panel.grid.minor=element_blank(),
    panel.grid.major=element_blank(),
    plot.title = element_text(size=14, face="bold",hjust = 0.4),
    axis.text = element_text(size = 12, face = "bold"),
    axis.title = element_text(colour="black",size = 14, face = "bold"),
    axis.text.x=element_text(hjust=0.5) ) + 
  ggtitle("TGFb_vs_Mock_24h") + 
  scale_color_manual(name="", values=c("red", "green", "black"), limits=c("Up", "Down", "NoDiff")) + 
  geom_vline(xintercept=c(-1, 1), lty=2, col="gray", lwd=0.5) + 
  geom_hline(yintercept=-log10(0.01), lty=2, col="gray", lwd=0.5)
valcano

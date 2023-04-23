# on uppmax folder cd /crex/proj/snic2021-23-14/Xue/PRJNA369563/rsem_ERCC2
rsem-generate-data-matrix *rsem.genes.results > output.matrix
cat output.matrix |sed -e "s/_rsem.genes.results//g" output.matrix >results.matrix
Scp xung0001@rackham.uppmax.uu.se:/crex/proj/snic2021-23-14/Xue/PRJNA369563/rsem_ERCC_2/results.matrix ./Desktop
# start to work on Rstudio
library(readr)
data <- read.delim("/Users/xung0001/Desktop/results.matrix",header = T)
# prepare sample table 

sample_info <- read.delim("/Users/xung0001/Desktop/sampleinfo.txt",header = F)
names(sample_info) <- c("GEO","sample","time")
sample_info$sample <- factor(c(rep("mock",6),rep("TGF-b",6),rep("TGFb_retinoic",6),rep("TGFb_retinoic_rapamyin",6),
rep("TGFb_butyrate",6),rep("positive",3)),levels = c("mock","TGF-b","TGFb_retinoic","TGFb_retinoic_rapamyin","TGFb_butyrate","positive"))
sample_info <- sample_info[order(sample_info$GEO),]
library(dplyr) # this allow us to use pipe function
row.names(sample_info) <- sample_info[,1]
sample_info<- sample_info[,-1,drop=FALSE]
sample_info %>% View() 
rownames(data) <- data$X
data<- data[,-1,drop=FALSE]
data1 <- data[,which(names(data) %in%row.names(sample_info))]
countdata <- data1
metadata <- sample_info
countdata <- round(countdata,digits=0) #make sure all the numbers have no digits
#Generate DGE
dds <- DESeqDataSetFromMatrix(countData = countdata,
 colData = metadata,
design = ~sample)
#check if there are some genes not expressed 
h <- counts(dds)
sel <-rowSums(counts(dds)) == 0
sprintf("%s%% percent (%s) of %s genes are not expressed",
     round(sum(sel) * 100/ nrow(h),digits=1),sum(sel), nrow(h))

# "14.1% percent (13) of 92 genes are not expressed"
 ## heatmap 
 select <- order(rowMeans(counts(dds, normalized=T)), decreasing=T)
 nt <- normTransform(dds)
 log2.norm.counts <- assay(nt)[select,]
 df <- as.data.frame(colData(dds)[, c("sample", "time")])
 pheatmap(log2.norm.counts, cluster_rows=T, show_rownames=F, cluster_cols=T, annotation_col=df, fontsize=6)

 # PCA analysis
  rld <- rlog(dds)
 pcaData <- plotPCA(rld, intgroup=c("sample", "time"), returnData=T)
 percentVar <- round(100*attr(pcaData, "percentVar"))
 pca <- ggplot(pcaData, aes(PC1, PC2, color=sample, shape=time)) + 
   geom_point(size=3) + 
   ggtitle("DESeq2 PCA") + 
   xlab(paste0("PC1: ", percentVar[1], "% variance")) + 
   ylab(paste0("PC2: ", percentVar[2], "% variance"))
 pca
 
# filter non-expressed genes 
dds <- dds[ rowSums(counts(dds)) > 1, ]
 dds <- DESeq(dds)
 res <- results(dds) # differential gene expression matrix 
 ## including single experiment data 
 resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)),by="row.names",sort=FALSE)
 resdata <- resdata[order(resdata$log2FoldChange,decreasing = T),]


# Gernerating vocalno plot with ggplot2
 library(ggplot2)
 resdata$change <- as.factor(
   ifelse(
     resdata$padj<0.01 & abs(resdata$log2FoldChange)>0,
     ifelse(resdata$log2FoldChange>1, "Up", "Down"),
     "NoDiff"
   )
 )
 valcano <- ggplot(data=resdata, aes(x=log2FoldChange, y=-log10(padj), color=change)) + 
   geom_point(alpha=0.8, size=1) + 
   theme_bw(base_size=15) + 
   theme(
     panel.grid.minor=element_blank(),
     panel.grid.major=element_blank()
   ) + 
   ggtitle("DESeq2 Valcano") + 
   scale_color_manual(name="", values=c("red", "green", "black"), limits=c("Up", "Down", "NoDiff")) + 
   geom_vline(xintercept=c(-1, 1), lty=2, col="gray", lwd=0.5) + 
   geom_hline(yintercept=-log10(0.01), lty=2, col="gray", lwd=0.5)
 valcano
# since there is no diffirentiated expressed genes in our spike-in data, we stop continue analyzing this data
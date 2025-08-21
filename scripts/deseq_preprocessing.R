


if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
#BiocManager::install("pcaExplorer")
#BiocManager::install("org.Hs.eg.db")
#BiocManager::install("DESeq2")
#library(org.Hs.eg.db)
#library(pcaExplorer)

check.packages <- function(pkg){
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if (length(new.pkg)) 
        install.packages(new.pkg, dependencies = TRUE)
    sapply(pkg, require, character.only = TRUE, warn.conflicts = FALSE, quietly=TRUE)
}

packages<-c("pcaExplorer", "tidyverse", "DESeq2", "dplyr", 
            "ggthemes", "RColorBrewer", "knitr", "ggpubr", 
            "gplots", "genefilter", "org.Hs.eg.db", "sva")

#check.packages(packages)


library('DESeq2')

setwd("./")
getwd()
options(show.error.locations = TRUE)


args = commandArgs(trailingOnly=TRUE)

countdata <- read.table(args[1], header = T, row.names=1)

# Convert to matrix
head(countdata)
#countdata %<>% dplyr::select(-1)
#countdata <- idplyr::select(-1)
coundata <- countdata[-c(1)]
countdata$Description <- NULL
head(countdata)
countdata <- as.matrix(countdata)

condition <- read.table(args[2], header = T)
condition$Treatment <- as.factor(condition$treatment)



control_group_name = args[3]
treatment_group_name = args[4]



library("sva")



(coldata <- data.frame(row.names=colnames(countdata[-1]), condition))
row.names(coldata) <- coldata$sample 

ddsFCseq <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~treatment)

ddsFCseq$sample
colData(ddsFCseq)



ddsFCseq_ko10 <- estimateSizeFactors(ddsFCseq)
ddsFCseq_ko10 <- ddsFCseq_ko10[ rowSums(counts(ddsFCseq_ko10, normalized = TRUE) >= 10) >= 3, ]  

ddsFCseq_ko10$treatment<- relevel(ddsFCseq_ko10$treatment, ref = control_group_name)

ddsFull <- DESeq(ddsFCseq_ko10)
#res <- results( ddsFull )
#head(res)
#boxplot(log10(assays(ddsFull)[["cooks"]]), range=0, las=2)

#plotDispEsts(ddsFull)



dat  <- counts(ddsFull, normalized = TRUE)
idx  <- rowMeans(dat) > 1
dat  <- dat[idx, ]
#mod  <- model.matrix(~ SeqRun+Treatment, colData(ddsFull))
mod  <- model.matrix(~Treatment, colData(ddsFull))
mod0 <- model.matrix(~   1, colData(ddsFull))
svseq <- svaseq(dat, mod, mod0, numSVmethod = "be")
#svseq$sv
plot(svseq$sv, pch=19, col="blue")
stripchart(svseq$sv[,1] ~ ddsFull$Treatment, vertical=TRUE)
#stripchart(svseq$sv[,2] ~ ddsFull$Treatment, vertical=TRUE)
#stripchart(svseq$sv[,3] ~ ddsFull$Treatment, vertical=TRUE)
#stripchart(svseq$sv[,4] ~ ddsFull$CellLine, vertical=TRUE)

ddssva <- ddsFull
ddssva$SV1 <- svseq$sv[,1]
#ddssva$SV2 <- svseq$sv[,2]
#ddssva$SV3 <- svseq$sv[,3]
#ddssva$SV4 <- svseq$sv[,4]
#ddssva$SV5 <- svseq$sv[,5]
design(ddssva) <- ~ SV1 + Treatment
#design(ddssva) <- ~ SV1 + SV2 + SV3 + Treatment




ddsFCseq_ko10 <- estimateSizeFactors(ddssva)
ddsFCseq_ko10 <- ddsFCseq_ko10[ rowSums(counts(ddsFCseq_ko10, normalized = TRUE) >= 10) >= 3, ]  

 ddsFCseq_ko10$Treatment<- relevel(ddsFCseq_ko10$Treatment, ref = control_group_name)

ddsFull2 <- DESeq(ddsFCseq_ko10)  # full stands for full dataset after filtering

rld <- rlog( ddsFull2 , blind=FALSE)
#head( assay(rld) )

sampleDists <- dist( t( assay(rld) ) )
#sampleDists

sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( rld$Sample_labs,
rld$Treatment, sep="-" )
colnames(sampleDistMatrix) <- NULL

#colours = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
#heatmap.2( sampleDistMatrix, trace="none", col=colours, labRow= rld$Treatment,
#           labCol= rld$Sample, cexRow = 0.8, cexCol = 0.8, srtCol=45, margins =c(6.5,8) )

#theme_set(theme_bw(24))

#mycols = c(control_group_name="#3265cc",
#           treatment_group_name="#d6604d")
## Treatment
#rld$Hour <- as.factor(rld$Hour)
#plotPCA(rld, intgroup = c("Treatment")) + 
#  geom_point(size = 4, alpha = 0.40) + 
#  scale_shape_manual(rld$Hour, values = c(18:22)) + 
  #ylim(-7,7) +
#  scale_colour_manual("Treatment", values = mycols, breaks = c(control_group_name, treatment_group_name))





  
#pcaobj <- stats::prcomp(t(assay(rld)))
#pcaExplorer::hi_loadings(pcaobj, topN = 10, whichpc = 1)
#pcaExplorer::hi_loadings(pcaobj, topN = 10, whichpc = 2)








annotation <- as.data.frame(read.table("all_readcounts.tsv",
                        header = T, row.names=1))
#annotation <- dplyr::select(annotation, 1)
#annotation <- select(annotation, 1)

annotation <- annotation[1]

colnames(annotation) <- c("gene_name")


#library('readr')


saveRDS(ddsFCseq_ko10, "svadata/ddsFCseq_ko10.RDS")
saveRDS(ddsFull, "svadata/ddsFull.RDS")
saveRDS(annotation, "svadata/annotation.RDS")
saveRDS(rld, "svadata/rld.RDS")

write.table(assay(rld), "svadata/rlog_transformed_reads.tsv", sep="\t", quote=FALSE)




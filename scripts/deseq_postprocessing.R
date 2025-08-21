check.packages <- function(pkg){
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if (length(new.pkg)) 
        install.packages(new.pkg, dependencies = TRUE)
    sapply(pkg, require, character.only = TRUE)
}

#BiocManager::install("pathview")

#BiocManager::install("gage")
#BiocManager::install("gageData")

# Usage example
packages<-c("pcaExplorer", "tidyverse", "DESeq2", "dplyr", 
            "ggthemes", "RColorBrewer", "knitr", "ggpubr", 
            "gplots", "genefilter", "org.Hs.eg.db", "knitr", "filesstrings", "fgsea")

#check.packages(packages)

library('dplyr')

library('pathview')
library('gage')
library('gageData')
library('ggplot2')
library('ggthemes')
library('DESeq2')

library('filesstrings')

#library("AnnotationDbi")
#library("org.Hs.eg.db")
mainDir <- "./"
subDir <- "./svadata"

#dir.create(file.path(subDir))
setwd(file.path(mainDir))


ddsFCseq <- readRDS("./svadata/ddsFCseq_ko10.RDS")
ddsFull <- readRDS("./svadata/ddsFull.RDS")
rld <- readRDS("./svadata/rld.RDS")
annotation <- readRDS( "svadata/annotation.RDS")




treatments <- factor(paste0(ddsFCseq$Treatment))
as.data.frame(levels(factor(treatments)))


args = commandArgs(trailingOnly=TRUE)

control_dwn <- args[1]
experimetal_up <- args[2]
comparison_name <- args[3]
#cell_type <- "polyA"
## experiment name
myexp <- paste(experimetal_up, "vs.", control_dwn, "_", comparison_name, sep='')


get_DE_genes <- function(x=diff_express,  sig = 0.05, lfc2 = 1) {
  res_tax <- merge(as.data.frame(diff_express), as.data.frame(counts(dds, normalized=T)), 
                   by='row.names', sort=F)
  names(res_tax)[1] <- 'GENE_ID'
  res_tax_sig = filter(res_tax, padj < sig & lfc2 < abs(log2FoldChange))
  res_tax_sig_exp1 <- res_tax_sig[order(res_tax_sig$padj),]
return(res_tax_sig)
}  



ddsFCseq2 <- ddsFCseq
ddsFCseq2$Treatment <- factor(paste0(ddsFCseq$Treatment))


ddsFCseq2 = ddsFCseq2[ , ddsFCseq2$Treatment == control_dwn | ddsFCseq2$Treatment == experimetal_up ]
#ddsFCseq2 = ddsFCseq2[ , ddsFCseq2$RNA == cell_type ]
colData(ddsFCseq2) <- colData(ddsFCseq2)[(1:3)]
as.data.frame(colData(ddsFCseq2))

#ddsFCseq2$Sample <- droplevels(ddsFCseq2$Sample)
#ddsFCseq2$Donor <- droplevels(ddsFCseq2$Donor)
ddsFCseq2$Treatment <- droplevels(ddsFCseq2$treatment)


design(ddsFCseq2) <- ~ treatment

ddsFCseq2$Treatment<- relevel(ddsFCseq2$treatment, ref = control_dwn)
dds <- DESeq(ddsFCseq2)


resultsNames(dds)
res_H4US <- results(dds, contrast=c("treatment", experimetal_up, control_dwn))  #order matters
rld <- rlog(dds)

#pcaExplorer::pcaplot(rld, intgroup = c("Treatment"), ntop = 500,
 #       pcX = 1, pcY = 2, title = "",
  #      ellipse = TRUE)


#pcaplot3d(rld, intgroup = c("Treatment"),ntop = 1000,
#        pcX = 1, pcY = 2, pcZ = 3)

#pcaobj <- prcomp(t(assay(rld)))
#hi_loadings(pcaobj, topN = 10, annotation = annotation)


## DEseq

diff_express=res_H4US
diff_express$GENE <- annotation$gene_name[match(row.names(diff_express), row.names(annotation))]

res_H4US_rep <- merge(as.data.frame(diff_express), as.data.frame(counts(dds, normalized=T)), 
                   by='row.names', sort=F)

csv_name <- paste0("./svadata/",myexp,".csv")
write.csv(as.data.frame(res_H4US_rep), file = csv_name, row.names = T, quote = F)


de_table <- get_DE_genes(x=diff_express,  sig = 0.05, lfc2 = 1)

#DT::datatable(de_table, filter = 'top',
#          extensions = 'Buttons', options = list(
#          dom = 'Bfrtip',
#          buttons = c('copy', 'csv', 'excel'),
#          pageLength = 5, autoWidth = T)
#              )

library('org.Hs.eg.db')

columns(org.Hs.eg.db)


diff_express$symbol = mapIds(org.Hs.eg.db,
                     keys=as.character(diff_express$GENE),
                     column="SYMBOL",
                     keytype="SYMBOL",
                     multiVals="first")
sum(is.na(diff_express$symbol))



diff_express$entrez = mapIds(org.Hs.eg.db,
                     keys=as.character(diff_express$GENE), 
                     column="ENTREZID",
                     keytype="SYMBOL",
                     multiVals="first")
sum(is.na(diff_express$entrez))


diff_express$name =   mapIds(org.Hs.eg.db,
                     keys=as.character(diff_express$GENE), 
                     column="GENENAME",
                     keytype="SYMBOL",
                     multiVals="first")
sum(is.na(diff_express$name))





#BiocManager::install("pathview")

data("kegg.sets.hs")
data("sigmet.idx.hs")
kegg.sets.hs = kegg.sets.hs[sigmet.idx.hs]
#head(kegg.sets.hs, 3)

foldchanges = diff_express$log2FoldChange
names(foldchanges) = diff_express$entrez
head(foldchanges)
## consider remove NAs
#foldchanges[!is.na(foldchanges)]
#sum(is.na(foldchanges))




# Get the results
keggres = gage(foldchanges, gsets=kegg.sets.hs, same.dir=TRUE)


dt<-as.data.frame(rbind(keggres$greater[1:6,], keggres$less[1:4,]))

ggplot(dt, aes(reorder(rownames(dt), stat.mean), stat.mean)) +
  geom_col(aes(fill=stat.mean > 0)) +
  coord_flip() +
  labs(x="KEEG Pathway", y="stat.mean()",
       title="") + 
  theme_minimal() + scale_fill_gdocs()

# save png  666 *273 inc

# Look at both up (greater), down (less), and statatistics.
#lapply(keggres, head)


### Top25 upregulated PATHWAYS

# Get the pathways
keggrespathways = data.frame(id=rownames(keggres$greater), keggres$greater) %>% 
  tbl_df() %>% 
  filter(row_number()<=25) %>% 
  .$id %>% 
  as.character()
keggrespathways



# Get the IDs.
keggresids = substr(keggrespathways, start=1, stop=8)
keggresids

# Define plotting function for applying later
plot_pathway = function(pid) pathview(gene.data=foldchanges, pathway.id=pid, species="hsa", new.signature=FALSE)

# plot multiple pathways (plots saved to disk and returns a throwaway list object)
tmp = sapply(keggresids, function(pid) pathview(gene.data=foldchanges, pathway.id=pid, species="hsa",
                                                low="#2166ac", mid="#bababa", high="#d73027", "*pathview.png"))


## Pathview analysis
#### do not forguet to change the folder path!!
pathdir <- paste0("./pathview_output_",myexp,"/")
dir.create(pathdir)
filenames <- Sys.glob(file.path(mainDir, "hsa*"))
file.move(filenames, pathdir, overwrite = TRUE)

filenames_pngs <- Sys.glob(file.path(pathdir, "*pathview.png"))
#knitr::include_graphics(filenames_pngs)




## FGSEA Analysis #####

#res2 <- as.tibble(diff_express) %>% 
#  dplyr::select(symbol, stat) %>% 
#  na.omit() %>% 
#  distinct() %>% 
#  group_by(symbol) %>% 
#  summarize(stat=mean(stat))
#res2


library(fgsea)

#ranks <- deframe(res2)
#head(ranks, 20)

# Load the pathways into a named list
pathways.hallmark <- fgsea::gmtPathways("h.all.v6.2.symbols.gmt")

# Look at them all if you want (uncomment)
# pathways.hallmark

# Show the first few pathways, and within those, show only the first few genes. 
#pathways.hallmark %>% 
#  head() %>% 
 # lapply(head)


#fgseaRes <- fgsea::fgsea(pathways=pathways.hallmark, stats=ranks, nperm=1000)


#fgseaResTidy <- fgseaRes %>%
#  as_tibble() %>%
#  arrange(desc(NES))

# Show in a nice table:
#fgseaResTidy %>% 
#  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
#  arrange(padj) %>% 
#  DT::datatable()


#fgseaResTidy %>%
#  filter(., NES > 2.2 | NES < -1.3) %>%
#ggplot(., aes(reorder(pathway, NES), NES)) +
#  geom_col(aes(fill=NES>2.2)) +
#  coord_flip() +
#  labs(x="Pathway", y="Normalized Enrichment Score",
#       title="Hallmark pathways NES from GSEA") + 
#  theme_minimal() + scale_fill_gdocs()


sessionInfo()

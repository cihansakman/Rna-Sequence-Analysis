library(clusterProfiler)
library(GSEABase)
library(org.Hs.eg.db)

#R Codes for finding ENTREZID
hs <- org.Hs.eg.db
#my.symbols <- c("CEP55","DSG3","TROAP","KIF14","PKM","UBE2T","GINS4","ERO1A","ITGB4","ARNTL2","SHCBP1","SEMA4B","CDCA8","AUNIP","OCIAD2","CCT6A","CCDC34","CASC9","TRIM59","GJB5") #Top 20 DEGs
#my.symbols <- c("ECE1","FCN3","CYP4B1","S100A6","MUC1","GLUL","ALOX5","SFTPA2","S100A4",
"SFTPA1","CTSD","HBB","SCGB1A1","CASP1","VWF","A2M","TNS2","STAT6","LTA4H","ALDH2") #Top 20 Variation Genes

#select(hs, keys = my.symbols,columns = c("ENTREZID", "SYMBOL"),keytype = "SYMBOL")


filename <- "c7.all.v7.1.entrez.gmt" 
gmtfile <- system.file(filename)
c6 <- read.gmt(gmtfile)
yourEntrezIdList<- c(55165,1830,10024,9928,5315,29089,84296,30001,3691,56938,79801,10509,55143,79000,132299,908,91057,101805492,286827,2709) #ENTREZID of DE genes

#yourEntrezIdList<- c(1889,8547,1580,6277,4582,2752,240,729238,6275,653509,1509,3043,7356,834,7450,2,23371,6778,4048,217) #ENTREZID of Annovar output

ImmunSigEnrich <- enricher(yourEntrezIdList, TERM2GENE=c6, pvalueCutoff = 0.05)
ImmunSigEnrich <- setReadable(ImmunSigEnrich, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
write.csv(ImmunSigEnrich,"MyImmunePathwayRelatedGenes.csv")
#write.csv(ImmunSigEnrich,"AnnovarMyImmunePathwayRelatedGenes.csv")

goEnrich<-enrichGO(gene= yourEntrezIdList,OrgDb= org.Hs.eg.db, ont= "ALL",pAdjustMethod="BH",pvalueCutoff = 0.05,readable= TRUE)
write.csv(goEnrich,"MyGORelatedGenes.csv")
#write.csv(goEnrich,"AnnovarMyGORelatedGenes.csv")

keggEnrich<-enrichKEGG(gene= yourEntrezIdList,organism= "hsa",pAdjustMethod="BH", pvalueCutoff = 0.1)
write.csv(keggEnrich,"MyKEGGRelatedGenes.csv")
#write.csv(keggEnrich,"AnnovarMyKEGGRelatedGenes.csv")
# Exit the R session


outfile="/media/data02/cihan2/enrichment/DE_GSEA_output.pdf"
#outfile="/media/data02/cihan2/enrichment/Annovar_GSEA_output.pdf"
pdf(file=outfile)
yourEntrezIdList <- sort(yourEntrezIdList, decreasing=TRUE)
library(enrichplot)
ego <- enrichGO(yourEntrezIdList, OrgDb = "org.Hs.eg.db", ont="ALL", readable=TRUE)
barplot(ego, showCategory=20)
dotplot(ego, showCategory=30)


barplot(keggEnrich, showCategory=20)
dotplot(keggEnrich, showCategory=30)

dev.off()

# Exit the R session
quit(save="no")



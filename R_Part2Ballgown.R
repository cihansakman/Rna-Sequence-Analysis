# Load libraries needed for this analysis
library(ballgown)
library(genefilter)
library(dplyr)
library(devtools)

# Define a path for the output PDF to be written
outfile="/media/data02/cihan2/enrichment/GSEA_output.pdf"

# Load phenotype data
pheno_data = read.csv("NORMAL_vs_LUSC.csv")

# Display the phenotype data
pheno_data

# Load the ballgown object from file
load('bg.rda')

# The load command, loads an R object from a file into memory in our R session. 
# You can use ls() to view the names of variables that have been loaded
ls()

# Print a summary of the ballgown object
bg

# Open a PDF file where we will save some plots. We will save all figures and then view the PDF at the end
pdf(file=outfile)

# Extract FPKM values from the 'bg' object
fpkm = texpr(bg,meas="FPKM")

# View the last several rows of the FPKM table
tail(fpkm)

# Transform the FPKM values by adding 1 and convert to a log2 scale
fpkm = log2(fpkm+1)

# View the last several rows of the transformed FPKM table
tail(fpkm)

# Create boxplots to display summary statistics for the FPKM values for each sample
boxplot(fpkm,col=as.numeric(pheno_data$type)+1,las=2,ylab='log2(FPKM+1)')


#function for finding the top column.
#Due to we have 227818 transcripts in bg we're going to search for the top Gene until 227818. When we found the top gene we're going to assign it's ID to top_gene_id variable. The ID of the top gene CEP55 is ENSG00000138180
for (val in 1: 227818)
{
    if (geneIDs(bg)[val] == "ENSG00000138180") {
	top_gene_id = val
	break
	}
}


# Display the transcript ID for a single row of data
ballgown::transcriptNames(bg)[top_gene_id]

# Display the gene name for a single row of data 
ballgown::geneNames(bg)[top_gene_id]

# Create a BoxPlot comparing the expression of a single gene for all replicates of both conditions
plot(fpkm[top_gene_id,] ~ as.factor(pheno_data$type), border=c(2,3), main=paste(ballgown::geneNames(bg)[top_gene_id],' : ', ballgown::transcriptNames(bg)[top_gene_id]),pch=19, xlab="Type", ylab='log2(FPKM+1)') 

# Add the FPKM values for each sample onto the plot 
#points(fpkm[top_gene_id,] ~ jitter(as.numeric(pheno_data$type)), col=as.numeric(pheno_data$type)+1, pch=16)

# Create a plot of transcript structures observed in each replicate and color transcripts by expression level
plotTranscripts(ballgown::geneIDs(bg)[top_gene_id], bg, main=c('CEP55 in sample SRR3474779'), sample=c('SRR3474779'))
plotTranscripts(ballgown::geneIDs(bg)[top_gene_id], bg, main=c('CEP55 in sample SRR3474795'), sample=c('SRR3474795'))
plotTranscripts(ballgown::geneIDs(bg)[top_gene_id], bg, main=c('CEP55 in sample SRR3474918'), sample=c('SRR3474918'))
plotTranscripts(ballgown::geneIDs(bg)[top_gene_id], bg, main=c('CEP55 in sample SRR3474972'), sample=c('SRR3474972'))
plotTranscripts(ballgown::geneIDs(bg)[top_gene_id], bg, main=c('CEP55 in sample SRR3475081'), sample=c('SRR3475081'))
plotTranscripts(ballgown::geneIDs(bg)[top_gene_id], bg, main=c('CEP55 in sample SRR3475089'), sample=c('SRR3475089'))
plotTranscripts(ballgown::geneIDs(bg)[top_gene_id], bg, main=c('CEP55 in sample SRR3475326'), sample=c('SRR3475326'))
plotTranscripts(ballgown::geneIDs(bg)[top_gene_id], bg, main=c('CEP55 in sample SRR3475332'), sample=c('SRR3475332'))
plotTranscripts(ballgown::geneIDs(bg)[top_gene_id], bg, main=c('CEP55 in sample SRR3475345'), sample=c('SRR3475345'))
plotTranscripts(ballgown::geneIDs(bg)[top_gene_id], bg, main=c('CEP55 in sample SRR3475362'), sample=c('SRR3475362'))
plotTranscripts(ballgown::geneIDs(bg)[top_gene_id], bg, main=c('CEP55 in sample SRR3475375'), sample=c('SRR3475375'))
plotTranscripts(ballgown::geneIDs(bg)[top_gene_id], bg, main=c('CEP55 in sample SRR3475377'), sample=c('SRR3475377'))

#plotMeans('TST',bg,groupvar="type",legend=FALSE)

# Close the PDF device where we have been saving our plots
dev.off()

# Exit the R session
quit(save="no")
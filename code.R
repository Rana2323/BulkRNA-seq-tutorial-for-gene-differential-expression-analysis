library(DESeq2)#each time you open your R

#loading the gene expression data 
data = as.matrix(read.csv("bladder_counts.csv", row.names=1))
#loading the phenotype data 
pheno = read.csv("bladder_sample_sheet.csv", row.names=1)

#explore the data, dim function give you the dimension of your data; how 
#many columns(samples) do you have and how many row(genes) do you have
dim(data)

#explore the data distribution using the histogram plot
hist(data, col = "orange", main="Histogram", breaks = 100)

#scaling the data using log2 transformation to better visulization
# we use (+1) to avoid the infinity character when we log zero valus 
hist(log2(data+1), col = "orange", main="Histogram")

#It is absolutely critical that the columns of the "data" and the rows of 
#the "pheno" (information about samples) are in the same order.DESeq2 will
#not make guesses as to which column of the count matrix belongsto which row
#of the column data, these must be provided to DESeq2 already in consistent order
pheno=pheno[colnames(data),]

#The deseq2 package require the count data values to be integers 
#save the gene names in a variable
genes=row.names(data)
#convert the data values to integers
data=apply(data,2,as.integer)
#view the data
head(data)
#rename the rows of the data
row.names(data)=genes
#specify how many conditions do you want to compare according to 
#the phenotypic table
cond1="Primary Tumor" 
cond2="Solid Tissue Normal"

#creat a deseq dataset object

dds= DESeqDataSetFromMatrix( countData = data , colData = pheno, design = ~ Sample.Type)
#run the deseq2 worflow
dds.run = DESeq(dds)
#specifying teh contrast (to make a res object based on two specific conditions)
res=results(dds.run, contrast = c("Sample.Type",cond1 ,cond2))
# remove nulls
res=as.data.frame(res[complete.cases(res), ])

#chose the statstical significant differentaily expressed genes (DEGs) based
#on the p adjusted value less than 0.05 and biological significance  based
#on the fold change more than 2
deseq.deg=res[res$padj < 0.05 & abs(res$log2FoldChange)>2,]

#export the Degs into your current folder for further analysthis
write.csv(as.matrix(deseq.deg),file="deseq.deg.csv", quote=F,row.names=T)

# install and load the pheatmap package
install.packages("pheatmap")
library(pheatmap)

# Filter the data for top 2000 differentially expressed genes
deseq.deg <- read.csv("deseq.deg.csv", row.names=1)
top_genes <- head(rownames(deseq.deg), 2000)
heatmap_data <- deseq.deg[top_genes, ]

# Scale the data
scaled_data <- scale(heatmap_data)

# Cluster the rows and columns
row_clusters <- hclust(dist(scaled_data))
col_clusters <- hclust(dist(t(scaled_data)))

# Create the heatmap
pheatmap(scaled_data, clustering_distance_rows = as.dist(1-cor(t(scaled_data))),
         clustering_distance_cols = as.dist(1-cor(scaled_data)),
         color = colorRampPalette(c("navy", "white", "firebrick"))(100),
         main = "Top 2000 Differentially Expressed Genes(Rana)",
         breaks = seq(-3, 3, length.out = 101)) # update the breaks argument

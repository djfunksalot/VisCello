library(irlba)
library(VisCello)
library(Matrix)
library(Biobase)
#gene count matrices 
raw_df <- read.table("data-raw/gene_by_cell_count_matrix.txt")
norm_df <- read.table("data-raw/gene_by_cell_count_matrix_batch_norm.txt")
#metadata
pmeta <- read.table("data-raw/cell_annotation.txt")
fmeta <- read.table("data-raw/gene_annotation.txt")
#sort and clean column/row names
raw_df<-raw_df[ , order(names(raw_df))]
pmeta <- pmeta[ order(row.names(pmeta)), ]
pmeta$name<-gsub("-",".",rownames(pmeta))
rownames(pmeta)<-gsub("-",".",rownames(pmeta))
pmeta$name<-NULL
fmeta$id<-rownames(fmeta)
#create expression set
eset <- new("ExpressionSet",
             assayData = assayDataNew("environment", exprs=Matrix(as.matrix(raw_df), sparse = T), norm_exprs = Matrix(as.matrix(norm_df), sparse = T)),
             phenoData =  new("AnnotatedDataFrame", data = pmeta),
             featureData = new("AnnotatedDataFrame", data = fmeta))
#eset for PT
#dimension reduction
umap_proj <- read.table("data-raw/umap.txt")
umap_proj <- umap_proj[ order(row.names(umap_proj)), ]
umap_proj$name<-gsub("-",".",rownames(umap_proj))
rownames(umap_proj) <- umap_proj$name
umap_proj$name<-NULL
#
#create cello object
cello <- new("Cello", name = "nonPT", idx = 1:ncol(eset)) 
cello@proj <- list("UMAP [2D]" = umap_proj)
clist <- list()
clist=list("nonPT"=cello)
#save to web directory
saveRDS(clist, "clist.rds",compress=FALSE) 
saveRDS(eset, "eset.rds",compress=FALSE) 
#gene count matrices 
pt_raw_df <- read.table("data-raw/gene_by_cell_count_matrix_rawUMI_withPT.txt")
pt_norm_df <- read.table("data-raw/gene_by_cell_count_matrix_batch_norm_withPT.txt")
#order columns
pt_raw_df <- pt_raw_df[,order(colnames(pt_raw_df))]
pt_norm_df <- pt_norm_df[,order(colnames(pt_norm_df))]
#phenotype metadata
pt_pmeta <- read.table("data-raw/cell_annotation_withPT.txt")
#genotype (feature) metadata
pt_fmeta <- read.table("data-raw/gene_annotation_withPT.txt")
#sort and clean column/row names
#order raw dataframe by barcode
pt_raw_df<-pt_raw_df[ , order(names(pt_raw_df))]
#pt_raw_df<-pt_raw_df[ , order(names(pt_norm_df))]
#order phenotype metadata by barcode
pt_pmeta <- pt_pmeta[ order(row.names(pt_pmeta)), ]
# replace occurance of "-" in barcode with "." to match gene count matrix
pt_pmeta$name<-gsub("-",".",rownames(pt_pmeta))
rownames(pt_pmeta)<-gsub("-",".",rownames(pt_pmeta))
pt_pmeta$name<-NULL
# don't know if we actually need this but the instructions say to do it
pt_fmeta$id<-rownames(pt_fmeta)
#subset pt_raw_df. use barcodes present in pt_norm_df
pt_raw_df<-pt_raw_df[,names(pt_norm_df)]
#create expression set
pt_eset <- new("ExpressionSet",
             assayData = assayDataNew("environment", exprs=Matrix(as.matrix(pt_raw_df), sparse = T), norm_exprs = Matrix(as.matrix(pt_norm_df), sparse = T)),
             phenoData =  new("AnnotatedDataFrame", data = pt_pmeta),
             featureData = new("AnnotatedDataFrame", data = pt_fmeta))
#dimension reduction
pt_umap_proj <- read.table("data-raw/tsne_withPT.txt")
pt_umap_proj <- pt_umap_proj[ order(row.names(pt_umap_proj)), ]
pt_umap_proj$name<-gsub("-",".",rownames(pt_umap_proj))
rownames(pt_umap_proj) <- pt_umap_proj$name
pt_umap_proj$name<-NULL
#create cello object
pt_cello <- new("Cello", name = "PT", idx = 1:ncol(pt_eset)) 
pt_cello@proj <- list("UMAP [2D]" = pt_umap_proj)
pt_clist <- list()
pt_clist=list("PT"=pt_cello)
#save to web directory
saveRDS(pt_clist, "pt_clist.rds",compress=FALSE) 
saveRDS(pt_eset, "pt_eset.rds",compress=FALSE) 

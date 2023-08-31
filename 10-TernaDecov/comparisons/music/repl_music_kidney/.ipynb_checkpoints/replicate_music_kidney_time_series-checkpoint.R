

library(data.table)

#########################################
## Load data
#########################################

# load reference (Park) data
park_data_path <- '/mnt/disk1/nbarkas/deconvolution_method/datasets/Park_etal/GSE107585_Mouse_kidney_single_cell_datamatrix.txt'
park_data <- fread(park_data_path)
park_data <- as.data.frame(park_data)


# Get a cell cluster factor
park_cluster <- as.character(t(park_data[1,])[,1])
names(park_cluster) <- colnames(park_data)
park_cluster <- park_cluster[-1]
park_cluster <- as.factor(park_cluster)

# remove the cluster number row
park_tmp <- park_data[-1,]

# Copy V1 into the rownames
rownames(park_tmp) <- as.character(park_tmp$V1)
head(rownames(park_tmp))

# Remove V1 columns
park_tmp <- park_tmp[,-1]
park_tmp[1:3,1:3]
park_counts <- as.matrix(park_tmp)
park_counts[1:3,1:3]

# load data to be deconvolved (craciun)
path.to.data <- '/mnt/disk1/nbarkas/deconvolution_method/datasets/craciun_etal/raw/'
files <- list.files(path.to.data)
names(files) <- nbHelpers::strpart(files,'.',1,fixed=TRUE)
full.files <- paste0(path.to.data, files)
names(full.files) <- names(files)
full.files
raw_files_read <- lapply(full.files, function(f) {read.table(f, header=TRUE)})
str(raw_files_read,1)
head(raw_files_read[[1]])

# check that all genes are in the same order
gene_names <- lapply(raw_files_read, function(x){x$ensemblid})
all(unlist(lapply(seq(1,length(gene_names)), function(i) { all(gene_names[[1]] == gene_names[[i]]) })))

# make data.frame
craciun_tmp <- do.call(cbind,lapply(raw_files_read, function(x) {x[2]}))
rownames(craciun_tmp) <- raw_files_read[[1]]$ensemblid
colnames(craciun_tmp) <- names(raw_files_read)
craciun_counts <- as.matrix(craciun_tmp)

n0 <- colnames(craciun_counts)
n1 <- nbHelpers::strpart(n0, '_',2)
replicate <- nbHelpers::strpart(n1, 'rep',2)
timemap <- c("Normal" =0,  "FA1day"=1,  "FA2day"=2,  "FA3day" =3, "FA7day"=7,  "FA14day"=14)
craniun_sample_time <- timemap[nbHelpers::strpart(n1, 'rep',1)]
names( craniun_sample_time ) <- colnames(craciun_counts)


#########################################
## Make Expression Set
#########################################

# Make eset for park
library(Biobase)

library(org.Mm.eg.db)

columns(org.Mm.eg.db)

ens.ids <- mapIds(x = org.Mm.eg.db, 
                  keys = rownames(park_counts), 
                  keytype = 'SYMBOL', 
                  column = 'ENSEMBL', multiVals = 'first')

length(ens.ids)
dim(park_counts)

rownames(park_counts) == names(ens.ids)

rownames(park_counts) <- ens.ids



head(rownames(park_counts))

park.es <- ExpressionSet(assayData = park_counts,
                         phenoData = AnnotatedDataFrame(data.frame(sampleID=names(park_cluster),
                                                                   cluster=park_cluster)))

# Filter and subsample park data
park.filter.es <- park.es[,colSums(exprs(park.es)) > 1000]
park.filter.sampled.es <- park.filter.es[,sample(ncol(exprs(park.filter.es)),10000)]



craciun.es <- ExpressionSet(assayData = craciun_counts,
                            phenoData = AnnotatedDataFrame(data.frame(time=craniun_sample_time,
                            sampleID=names(craniun_sample_time))))


library(MuSiC)
library(xbioc)

select.ct <- c("3", "5", "6", "10", "13", "14", "4", "1", "7", "12", "11", "15", "9", "8", "16", "2")

assayData(park.es)




# Estimate cell type proportions
Est.prop = music_prop(bulk.eset = craciun.es, 
                      sc.eset = park.filter.sampled.es, 
                      clusters = 'cluster',
                      samples = 'sampleID',
                      select.ct = select.ct, verbose = F)

rownames(craciun.es)

rownames(park.es)

names(Est.prop.GSE50244)

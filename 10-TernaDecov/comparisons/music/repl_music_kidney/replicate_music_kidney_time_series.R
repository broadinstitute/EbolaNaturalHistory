



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

# Remove V1 columnsb
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

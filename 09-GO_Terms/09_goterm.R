library(topGO)
library(ALL)
data(ALL)
data(geneList)
library(biomaRt)
library(Rgraphviz)
library(ggplot2)
library(ggpubr)
library(dplyr)

colMap <- function(x) {
  .col <- rep(rev(heat.colors(length(unique(x)))), time = table(x))
  return(.col[match(1:length(x), order(x))])
}

# get top genes by filtering pval < 0.01, get top N of them
topN = 100
select_genes <- function(genelist) {
  # expects a named vector of pvals
  gene_ordered = all_genes_org[order(all_genes_org, decreasing = FALSE)]
  sig = (genelist < gene_ordered[topN + 1]) # to get top N genes
  return(sig)
}

# this time we want to look at the time DE genes
res <- read.csv('~/Desktop/Harvard/SabetiLab/Ebola_DE/EbolaFiles/DE_Compiled.csv') # adjust to path

# examine p value data
head(res$Infection.pvalue) 
head(res$impulseTOsigmoid_p) 
head(res$sigmoidTOconst_p) 
head(res$p) 

# need to map go terms to genes
bm <- useMart("ensembl")
bm <- useDataset("mmulatta_gene_ensembl", mart=bm)
EG2GO <- getBM(mart=bm, attributes=c('ensembl_gene_id','external_gene_name','go_id'))
head(EG2GO,15)
EG2GO <- EG2GO[EG2GO$go_id != '',]
geneID2GO <- by(EG2GO$go_id,
                EG2GO$ensembl_gene_id,
                function(x) as.character(x))
GO2geneID <- by(EG2GO$ensembl_gene_id,
                EG2GO$go_id,
                function(x) as.character(x))
head(geneID2GO)
head(GO2geneID)

organs = c("Adrenal Gland", "Skin:Rash", "Sex Organ:Sex-Organ", "Sex Organ:Ovary", 
           "Kidney", "Skin:Non-Rash", "Brain:Brain-Wh", "Brain:Brain-Gr",
           "Whole blood", "Spleen", "Lymph node:LN-AX-R", "Lymph node:LN-ING-L",
           "Lymph node:LN-MES", "Liver", "Lung", "Sex Organ:Testis") #ordered
myList = list()

for(x in 1:length(organs)) {
  
  # set up genes
  curr_org = organs[x]
  print(curr_org) 
  res_org = subset(res, tissue == curr_org & !is.na(p))
  all_genes_org = res_org$p
  names(all_genes_org) = sub("\\..*", "", as.vector(res_org$X)) 
  
  # build the GO object
  sampleGOdata_org <- new("topGOdata",
                          description = "Simple session", ontology = "BP",
                          allGenes = all_genes_org, 
                          geneSel = select_genes,
                          nodeSize = 10, annot = annFUN.gene2GO,
                          gene2GO = geneID2GO)
  
  # tests
  resultFisher_org <- runTest(sampleGOdata_org, algorithm = "classic", statistic = "fisher")
  resultKS_org <- runTest(sampleGOdata_org, algorithm = "classic", statistic = "ks") 
  resultKS.elim_org <- runTest(sampleGOdata_org, algorithm = "elim", statistic = "ks") 
  
  # organize results
  ntop = 3 # reduce number of terms considered
  allRes_org <- GenTable(sampleGOdata_org, classicFisher = resultFisher_org,
                         classicKS = resultKS_org, elimKS = resultKS.elim_org,
                         orderBy = "elimKS", ranksOf = "classicFisher", 
                         topNodes = ntop, numChar = 99)
  allRes_org$KS <- as.numeric(allRes_org$classicKS)
  allRes_org <- allRes_org[allRes_org$KS < 0.05,] # filter terms for KS p < 0.05
  allRes_org <- allRes_org[,c("GO.ID","Term","KS")]
  allRes_org$Organ = curr_org
  
  myList[[curr_org]] = allRes_org
  
}

for (j in 1:length(myList)) {
  myList[[j]] = subset(myList[[j]], !is.na(GO.ID))
}

ggdata <- bind_rows(myList, .id = "column_label")
ggdata$enrichment = -log10(ggdata$KS)
ggdata$expected = 0
for(y in 1:nrow(ggdata)) {
  ggdata$expected[y] = as.numeric(termStat(sampleGOdata_org, ggdata$GO.ID[y])[3])
}
# write this data to table so we don't need to run through all the computations again!!
write.csv(ggdata, '~/Desktop/GO_balloon_input_20220504_N100.csv', sep = '\t')
ggdata <- read.csv('~/Desktop/GO_balloon_input_20220504_N100.csv', sep = '\t')

ggdata$Organ = factor(ggdata$Organ, levels = c("Adrenal Gland", "Skin:Rash", "Sex Organ:Sex-Organ", "Sex Organ:Ovary", 
                                               "Kidney", "Skin:Non-Rash", "Brain:Brain-Wh", "Brain:Brain-Gr",
                                               "Whole blood", "Spleen", "Lymph node:LN-AX-R", "Lymph node:LN-ING-L",
                                               "Lymph node:LN-MES", "Liver", "Lung", "Sex Organ:Testis"))
ggdata$Term = factor(ggdata$Term, levels = unique(ggdata[order(ggdata$enrichment, decreasing = TRUE), ]$Term))


ggballoonplot(ggdata, x = "Organ", y = "Term", size = "expected",
              fill = "enrichment",
              size.range = c(1,3.5),
              ggtheme = theme_bw()) +
  #scale_fill_viridis_c(option = "C")
  scale_fill_gradientn(colors = colorRampPalette(c("navy", "white","firebrick3"))(100))
ggsave("~/Desktop/panel2d_N100_w6.pdf",width=6,height=3) # try 8, 6, 4






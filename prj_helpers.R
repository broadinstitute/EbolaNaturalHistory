#' Get the specified substring of string delimited by a character
#' @param x the strings to breakdown now (can be vector)
#' @param split the character to breakdown by
#' @param integer specifying which substing to get
#' @return extracted substrings
#' @export strpart
strpart <- function(x, split, n, fixed=FALSE) {
    sapply(strsplit(as.character(x),split,fixed=fixed),'[',n)
}

#' Convert NA values in x to the value specified by val
#' @export NA2VALUE
NA2VALUE <- function(x, val) {
  x[is.na(x)] <- c(val); x
}

#' Convert NA values in x to FALSE
#' @export NA2FALSE
NA2FALSE <- function(x) {
  NA2VALUE(x,FALSE)
}

# a function to consistently fix column name
fix_colnames <- function(column_names) {
  make.names(gsub(".bam","",gsub("bams_renamed/","",column_names)))
}


# getMethod(plotPCA, "DESeqTransform")


customPlotPCA <- function (object, intgroup = "condition", ntop = 500, returnData = FALSE,fixCoord=TRUE) {
    rv <- rowVars(assay(object))
    select <- order(rv, decreasing = TRUE)[seq_len(min(ntop,
                                                       length(rv)))]
    pca <- prcomp(t(assay(object)[select,]))
    percentVar <- pca$sdev ^ 2 / sum(pca$sdev ^ 2)
    if (!all(intgroup %in% names(colData(object)))) {
      stop("the argument 'intgroup' should specify columns of colData(dds)")
    }
    intgroup.df <- as.data.frame(colData(object)[, intgroup,
                                                 drop = FALSE])
    group <- if (length(intgroup) > 1) {
      factor(apply(intgroup.df, 1, paste, collapse = ":"))
    }
    else {
      colData(object)[[intgroup]]
    }
    d <-
      data.frame(
        PC1 = pca$x[, 1],
        PC2 = pca$x[, 2],
        group = group,
        intgroup.df,
        name = colnames(object)
      )
    if (returnData) {
      attr(d, "percentVar") <- percentVar[1:2]
      return(d)
    }
    
    p <- ggplot2::ggplot(data = d, ggplot2::aes_string(x = "PC1", y = "PC2", color = "group")) +
      ggplot2::geom_point(size = 3) + 
      ggplot2::xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) + 
      ggplot2::ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance"))
    
    if (fixCoord) {
        p <- p + ggplot2::coord_fixed()
    }
    
    p
}


getPCA <- function (object, ntop = 500) {
  rv <- rowVars(assay(object))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca <- prcomp(t(assay(object)[select,]))
  percentVar <- pca$sdev ^ 2 / sum(pca$sdev ^ 2)
  data.frame(
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    PC3 = pca$x[, 3],
    PC4 = pca$x[, 4],
    PC5 = pca$x[, 5]
  )
}

###
plotGene_custom01 <- function(gene,dds) {
    gene_symbol <- rowData(dds)[gene,]$external_gene_name
    
    plotCounts(dds,gene,intgroup = c('tissue','infected'), returnData = TRUE) %>% 
      ggplot(aes(x=tissue,fill=infected,y=log10(count+1) )) +
        geom_boxplot() + theme_bw() + ggtitle(paste0(gene, ' (', gene_symbol, ')'))
}

plotTopNGenes <- function(res, dds, n=10, returnData = FALSE, pseudocount = 1, facet_by = 'gene'){
    top_genes <- rownames(head(res,n=n))
    d <- do.call(rbind,lapply(top_genes, function(g){
        d1 <- plotCounts(dds,g,intgroup = c('tissue','infected'), returnData = TRUE)
        d1$gene <- c(g)
        d1
    }))
    d$symbol <- rowData(dds)$external_gene_name[match(d$gene, rowData(dds)$Geneid)]
    d$label <- paste0(d$symbol,' (',d$gene,')')
    if (returnData) {
        return(d)
    } else {
        if (facet_by == 'gene') {
            d %>% ggplot(aes(x=tissue,fill=infected,y=log2(count + pseudocount) )) +
                 geom_boxplot() + 
                theme_bw() + 
                theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
            scale_y_continuous(name="log2(Normalized Counts)") + facet_wrap(~label) -> p
        } else {
            d %>% ggplot(aes(x=label,fill=infected,y=log2(count + pseudocount) )) +
                 geom_boxplot() + 
                theme_bw() + 
                theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
            scale_y_continuous(name="log2(Normalized Counts)") + facet_wrap(~tissue) -> p
        }
        p
    }
}

post_process_res <- function(res,se) {
  res.ordered <- as.data.frame(res[order(res$padj),])
  ordered.row.data <- as.data.frame(rowData(se)[match(rownames(res.ordered),rownames(rowData(se))),])
  data.frame(cbind(res.ordered,ordered.row.data))
}

plotVolcano <- function(res,cutoff=0.05,xlim=NULL,n.label=10) {
  # Volcano Plot
  res %>% 
    ggplot(aes(x=log2FoldChange,y=-log10(padj),color=padj < cutoff)) + 
    geom_point()  + 
    scale_x_continuous(lim=xlim) + 
    theme_bw() + 
    scale_color_manual(values=c('black','red')) + 
    ggrepel::geom_label_repel(data=head(res,n=n.label), aes(label=external_gene_name)) +
    theme(legend.position = 'none')
}


runMacaqueGO <- function(res, test.cats=c("GO:CC", "GO:BP", "GO:MF") ) {
  # get the genes
  genes <- as.integer(res$padj[!is.na(res$log2FoldChange)] < .05)
  names(genes) <- rownames(res)[!is.na(res$log2FoldChange)]
  genes <- genes[!is.na(genes)]
  # Get Gene Lengths
  gene.lengths <- res$Length
  names(gene.lengths) <- rownames(res)
  # Manually provide the gene lengths to nullp as there is no db for rheMac10
  pwf <- nullp(genes, "rheMac10", "ensGene", bias.data = gene.lengths[names(genes)])
  # Remove ensembl transcript versions
  rownames(pwf) <- nbHelpers::strpart(rownames(pwf), '.', 1, fixed = T)
  # An example of the gene2cat data
  gene2cat <-getgo(genes = nbHelpers::strpart(names(genes),split = '.',fixed = T,n = 1),genome = 'rheMac10',id = 'ensGene',fetch.cats=test.cats)
  # Run the go term analysis
  goseq_res1 <- goseq(pwf = pwf,genome = 'rheMac10',gene2cat = gene2cat, test.cats = test.cats )
  invisible(goseq_res1)
}

#' Convert NA values in x to the value specified by val
#' @export NA2VALUE
NA2VALUE <- function(x, val) {
  x[is.na(x)] <- c(val); x
}

#' Convert NA values in x to FALSE
#' @export NA2FALSE
NA2FALSE <- function(x) {
  NA2VALUE(x,FALSE)
}

#' Get a vector of the names of an object named by the names themselves
#' @description Get a named vector of the names of an object. This is useful
#' with lapply when passing names of objects as it ensures that the output list
#' is also named
#' @param g an objects on which we can call names()
#' @export namedNames
namedNames <- function(g) {
  n <- names(g)
  names(n) <- n;
  n
}


draw.double.venn.fromVectors <- function(v1,v2,...) {
    v1 <- unique(v1)
    v2 <- unique(v2)
    
    n1 <- length(v1)
    n2 <- length(v2)
    
    n12 <- length(intersect(v1,v2))
    
    VennDiagram::draw.pairwise.venn(area1 = n1,area2=n2,cross.area = n12, ...)
}


get_genes_in_go_cat <- function(gocat) {
    getGenesInGO <- function(goid) {select(org.Mmu.eg.db, keys=goid, columns=c("ENSEMBL"), keytype='GOALL')}
    g_tmp <- intersect(getGenesInGO(gocat)$ENSEMBL, strpart(names(genes[genes]),'.',1,fixed=T))
    brain_de$brain_gray[brain_de$brain_gray$Geneid_noversion %in% g_tmp,]$external_gene_name -> genes_of_interest_tmp
    genes_of_interest_tmp
}

get_genes_in_go_cat_ens <- function (gocat) 
{
    getGenesInGO <- function(goid) {
        select(org.Mmu.eg.db, keys = goid, columns = c("ENSEMBL"), 
            keytype = "GOALL")
    }
    getGenesInGO(gocat)$ENSEMBL
}
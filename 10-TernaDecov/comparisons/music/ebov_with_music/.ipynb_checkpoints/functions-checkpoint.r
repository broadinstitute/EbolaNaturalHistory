library(tidyverse)

plot_ebov_deconvolution <- function(music_results, se.obj, return_data = FALSE) {
    # prepare results into long format
    res1 <- as.data.frame(music_results$Est.prop.weighted)
    res1$file_identifiers <- rownames(res1)
    rownames(res1) <- NULL
    res1 %>% pivot_longer(cols=!file_identifiers,names_to="celltype") -> res1_long
    
    # prepare metadata
    meta <- colData(se.obj)
    meta <- meta[,c('file_identifiers','id.individual','dpi_time_factor')]
    res1_merged <- merge(meta, res1_long)
    res1_merged <- as.data.frame(res1_merged)
    
    if (return_data) {
        res1_merged
    } else {
        # generate plot
        p <- ggplot(res1_merged, aes(x=dpi_time_factor,y=value)) + geom_point() + facet_wrap(~celltype)
        p
    }

}
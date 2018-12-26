#!/usr/bin/env Rscript

echo <- function(...) {
    cat(..., "\n")
}


print.usage <- function() {
    echo("Perform the DESeq2 analyze.")
    echo()
    echo("Usage:")
    echo("   $ ./deseq2_analyze.R EXP_TABLE CONFIG.json")
    echo()
}


read.exptable <- function(path) {
    exptable <- read.table(path, header = TRUE, quote = '\t')
    rownames(exptable) <-exptable$gene
    exptable <- exptable[,-1]
    return(exptable)
}


read.config <- function(path) {
    library("rjson")
    json.data <- fromJSON(file=path)
    if (length(json.data$conditions) != length(json.data$samples)) {
        stop("Config file is not valid, conditions and samples must have same size.")
    }
    return(json.data)
}


main <- function() {

    # Parse arguments
    args <- commandArgs(trailingOnly=TRUE)
    if (length(args) < 2) {
        print.usage()
        stop()
    }
    exptable.path <- args[1]
    config.path <- args[2]

    # read expression table
    exptable <- read.exptable(exptable.path)
    echo("Shape of full expression table:")
    echo(dim(exptable))
    echo()

    # read config json file
    config <- read.config(config.path)
    echo("samples:")
    echo(config$samples, sep="\t")
    echo()
    echo("conditions:")
    echo(config$conditions, sep="\t")
    echo()
    echo("comparisons:")
    for (comp in config$comparisons) {
        cat(comp, sep="_VS_")
        cat("\t")
    }
    echo("\n")

    # take subset of the expression table
    exptable <- exptable[config$samples]
    echo("Shape of expression table subset:")
    echo(dim(exptable))

}

main()


#library(DESeq2)

#condition <- c(rep("bend3 day0 uninfect", 2), rep("bend3 day0 infect", 2), 
#               rep("bend3 day3 uninfect", 2), rep("bend3 day3 infect", 2),
#               rep("raw day0 uninfect", 2), rep("raw day0 infect", 2),
#               rep("raw day3 uninfect", 2), rep("raw day3 infect", 2))
#table.condition <- data.frame(name = colnames(merged), condition = condition)
#rownames(table.condition) <- table.condition$name
#
## translate to deseq2 data format.
#dds <- DESeqDataSetFromMatrix(merged, colData=table.condition, design= ~ condition)
#dds <- dds[ rowSums(counts(dds)) > 1, ] # filter out the zero count genes
#
## do PCA analyze
#library(ggplot2)
#rld <- rlog(dds)
## plot 
#data.pca <- plotPCA(rld, intgroup=c("condition", "name"), returnData=TRUE)
#percentVar <- round(100 * attr(data.pca, "percentVar"))
#fig.pca <- ggplot(data.pca, aes(PC1, PC2, color=condition, shape=name)) +
#    scale_shape_manual(values=seq(0,15)) +
#    geom_point(size=5) +
#    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
#    ylab(paste0("PC2: ",percentVar[2],"% variance")) 
#ggsave("pca.pdf")
#
## DE analyze
#library(DESeq)
#dds <- DESeq(dds)
#infect_bend3d3 <- results(dds, contrast=c("condition", "bend3 day3 uninfect", "bend3 day3 infect"))
#infect_rawd0   <- results(dds, contrast=c("condition", "raw day0 uninfect", "raw day0 infect"))
#infect_bend3d0 <- results(dds, contrast=c("condition", "bend3 day0 uninfect", "bend3 day0 infect"))
#infect_rawd3   <- results(dds, contrast=c("condition", "raw day3 uninfect", "raw day3 infect"))
#time_bend3  <- results(dds, contrast=c("condition", "bend3 day0 infect", "bend3 day3 infect"))
#time_bend3u <- results(dds, contrast=c("condition", "bend3 day0 uninfect", "bend3 day3 uninfect"))
#time_raw    <- results(dds, contrast=c("condition", "raw day0 infect", "raw day3 infect"))
#time_rawu   <- results(dds, contrast=c("condition", "raw day0 uninfect", "raw day3 uninfect"))
#list.res <- list(infect_bend3d3, infect_rawd0, infect_bend3d0, infect_rawd3, 
#                 time_bend3, time_bend3u, time_raw, time_rawu)
#names <- c("infect_bend3d3", "infect_rawd0", "infect_bend3d0", "infect_rawd3",
#           "time_bend3", "time_bend3u", "time_raw", "time_rawu")
#
#lapply(seq_along(list.res),
#       function(x) write.table(list.res[[x]], paste(names[x], ".txt", sep=""),
#                                         col.names=TRUE, row.names=TRUE, sep="\t",
#                                         quote=FALSE))
#
## plot MA
#library(geneplotter)
#lapply(seq_along(list.res),
#        function(x) {
#            pdf(paste("MA_", names[x], ".pdf", sep=""))
#            plotMA(list.res[[x]], main="DESeq2")
#            dev.off()
#        })
#
## plot Heatmap
#library("pheatmap")
#select <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)[1:1000]
#nt <- normTransform(dds) # defaults to log2(x+1)
#log2.norm.counts <- assay(nt)[select,]
#df <- as.data.frame(colData(dds)[,c("name","condition")])
#pdf('heatmap1000.pdf',width = 6, height = 7)
#pheatmap(log2.norm.counts, cluster_rows=TRUE, show_rownames=FALSE,
#cluster_cols=TRUE, annotation_col=df)
#dev.off()
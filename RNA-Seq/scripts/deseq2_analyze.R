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


plot.pca <- function(dds, fig.path) {
    library(ggplot2)
    rld <- rlog(dds)
    data.pca <- plotPCA(rld, intgroup=c("condition", "name"), returnData=TRUE)
    percentVar <- round(100 * attr(data.pca, "percentVar"))
    fig.pca <- ggplot(data.pca, aes(PC1, PC2, color=condition, shape=name)) +
        scale_shape_manual(values=seq(0,15)) +
        geom_point(size=5) +
        xlab(paste0("PC1: ",percentVar[1],"% variance")) +
        ylab(paste0("PC2: ",percentVar[2],"% variance")) 
    ggsave(fig.path)
}


plot.heatmap <- function(dds, fig.path, row.num) {
    if (missing(row.num)) {
        row.num <- 1000
    }
    library("pheatmap")
    dds <- estimateSizeFactors(dds)
    select <- order( rowMeans( counts(dds, normalized=TRUE) ), decreasing=TRUE )[1:row.num]
    nt <- normTransform(dds) # defaults to log2(x+1)
    log2.norm.counts <- assay(nt)[select,]
    df <- as.data.frame(colData(dds)[,c("name","condition")])
    pdf(fig.path, width = 6, height = 7)
        pheatmap(log2.norm.counts, cluster_rows=TRUE, show_rownames=FALSE,
        cluster_cols=TRUE, annotation_col=df)
    dev.off()
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

    # convert to DESeq DataSet
    library(DESeq2)
    table.condition <- data.frame(name = colnames(exptable), condition=config$conditions)
    rownames(table.condition) <- table.condition$name
    dds <- DESeqDataSetFromMatrix(exptable, colData=table.condition, design= ~ condition)
    dds <- dds[ rowSums(counts(dds)) > 1, ] # filter out the zero count genes

    # output normalized counts
    dds <- estimateSizeFactors(dds)
    counts.normalized <- counts(dds, normalized=TRUE)
    write.table(counts.normalized, "./normalized_counts.tsv", col.names=NA, row.names=TRUE, sep="\t", quote=FALSE)

    # PCA plot
    plot.pca(dds, "pca.pdf")

    # plot heatmap
    plot.heatmap(dds, "heatmap.pdf")

    # DESeq analyze
    dds <- DESeq(dds)
    deseq.analyze <- function(compare) {
        results(dds, contrast=c("condition", compare[2], compare[1]))
    }
    deseq.results <- lapply(config$comparisons, deseq.analyze)
    deseq.save <- function(idx) {
        compare.name <- paste(config$comparisons[[idx]], collapse="-")
        out.filename <- paste(compare.name, ".tsv", sep="")
        write.table(deseq.results[idx], out.filename,
                                        col.names=NA, row.names=TRUE, sep="\t",
                                        quote=FALSE)
    }
    save_ <- lapply(seq_along(deseq.results), deseq.save)

    # plot MA
    library(geneplotter)
    plot.MA <- function(idx) {
        compare.name <- paste(config$comparisons[[idx]], collapse="-")
        fig.path <- paste("MA_", compare.name, ".pdf", sep="")
        pdf(fig.path)
            plotMA(deseq.results[[idx]], main="DESeq2 MA plot")
        dev.off()
    }
    plot_ <- lapply(seq_along(deseq.results), plot.MA)
}

main()

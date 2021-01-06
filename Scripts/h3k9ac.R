## ----style, echo=FALSE, results='hide', message=FALSE-------------------------
library(BiocStyle)
library(knitr)
opts_chunk$set(error=FALSE, message=FALSE, warning=FALSE, cache=TRUE)

## -----------------------------------------------------------------------------
library(chipseqDBData)
acdata <- H3K9acData()
acdata

## ----mapstat------------------------------------------------------------------
library(Rsamtools)
diagnostics <- list()
for (b in seq_along(acdata$Path)) {
    bam <- acdata$Path[[b]]
    total <- countBam(bam)$records
    mapped <- countBam(bam, param=ScanBamParam(
        flag=scanBamFlag(isUnmapped=FALSE)))$records
    marked <- countBam(bam, param=ScanBamParam(
        flag=scanBamFlag(isUnmapped=FALSE, isDuplicate=TRUE)))$records
    diagnostics[[b]] <- c(Total=total, Mapped=mapped, Marked=marked)
}

diag.stats <- data.frame(do.call(rbind, diagnostics))
rownames(diag.stats) <- acdata$Name
diag.stats$Prop.mapped <- diag.stats$Mapped/diag.stats$Total*100
diag.stats$Prop.marked <- diag.stats$Marked/diag.stats$Mapped*100
diag.stats

## -----------------------------------------------------------------------------
library(BiocFileCache)
bfc <- BiocFileCache("local", ask=FALSE)
black.path <- bfcrpath(bfc, file.path("https://www.encodeproject.org",
    "files/ENCFF547MET/@@download/ENCFF547MET.bed.gz"))

## -----------------------------------------------------------------------------
library(rtracklayer)
blacklist <- import(black.path)
blacklist

## ---- message=FALSE-----------------------------------------------------------
library(csaw)
standard.chr <- paste0("chr", c(1:19, "X", "Y"))
param <- readParam(minq=20, discard=blacklist, restrict=standard.chr)

## -----------------------------------------------------------------------------
x <- correlateReads(acdata$Path, param=reform(param, dedup=TRUE))
frag.len <- maximizeCcf(x)
frag.len

## ----ccfplot, fig.cap="Cross-correlation function (CCF) against delay distance for the H3K9ac data set. The delay with the maximum correlation is shown as the red line."----
plot(1:length(x)-1, x, xlab="Delay (bp)", ylab="CCF", type="l")
abline(v=frag.len, col="red")
text(x=frag.len, y=min(x), paste(frag.len, "bp"), pos=4, col="red")

## ----directional, echo=FALSE, fig.asp=0.4, fig.wide=TRUE, fig.cap="Directional extension of reads by the average fragment length `ext` in single-end ChIP-seq data. Each extended read represents an imputed fragment, and the number of fragments overlapping a window of a given `width` is counted."----
par(mar=c(0,0,0,0), cex=1.2)
plot(0,0,xlim=c(-2, 12), ylim=c(-1, 3.5), type="n", axes=FALSE, xlab="", ylab="", xpd=TRUE)

# Filling the left region.
rect(0, 0, 10, 1, col="grey")
lines(c(0, 0, 10, 10), 
    c(-.1, -.5, -.5, -.1))
lines(c(5, 5),
    c(-.5, -.7))
text(5, -.7, "width",pos=1)

# Adding the left read.
lines(c(-2, 2.5), 
    c(1.5, 1.5), 
    col=rgb(1, 0.7, 0.7), lwd=5)
lines(c(-2, -0.5, -0.8),
    c(1.5, 1.5, 1.8),
    col="red", lwd=5)
lines(c(-2, -2, 2.5, 2.5),
    c(1.8, 2, 2, 1.8))
text(-2, 1.45, "forward\nread", pos=1)

# Adding the right read.
lines(c(12, 7.5),
    c(1.8, 1.8),
    col=rgb(0.7, 0.7, 1), lwd=5)
lines(c(12, 10.5, 10.8),
    c(1.8, 1.8, 2.1),
    col="blue", lwd=5)
lines(c(12, 12, 7.5, 7.5),
    c(2.1, 2.3, 2.3, 2.1))
text(12, 1.75, "reverse\nread", pos=1)

# Describing the extension.
lines(c(0, 0, 3.1),
    c(2, 2.8, 2.8))
lines(c(10, 10, 6.9),
    c(2.3, 2.8, 2.8))
text(5, 2.8, "fragment length (ext)")

## -----------------------------------------------------------------------------
win.data <- windowCounts(acdata$Path, param=param, width=150, ext=frag.len)
win.data

## -----------------------------------------------------------------------------
bins <- windowCounts(acdata$Path, bin=TRUE, width=2000, param=param)
filter.stat <- filterWindowsGlobal(win.data, bins)
min.fc <- 3
keep <- filter.stat$filter > log2(min.fc)
summary(keep)

## ----bghistplot, fig.cap="Histogram of average abundances across all 2 kbp genomic bins. The filter threshold is shown as the red line."----
hist(filter.stat$back.abundances, main="", breaks=50,
    xlab="Background abundance (log2-CPM)")
threshold <- filter.stat$abundances[1] - filter.stat$filter[1] + log2(min.fc)
abline(v=threshold, col="red")

## -----------------------------------------------------------------------------
filtered.data <- win.data[keep,]

## ----trendplot, fig.cap="Abundance-dependent trend in the log-fold change between two H3K9ac samples (mature B over pro-B), across all windows retained after filtering. A smoothed spline fitted to the log-fold change against the average abundance is also shown in red."----
win.ab <- scaledAverage(filtered.data)
adjc <- calculateCPM(filtered.data, use.offsets=FALSE)
logfc <- adjc[,4] - adjc[,1]
smoothScatter(win.ab, logfc, ylim=c(-6, 6), xlim=c(0, 5),
    xlab="Average abundance", ylab="Log-fold change")

lfit <- smooth.spline(logfc~win.ab, df=5)
o <- order(win.ab)
lines(win.ab[o], fitted(lfit)[o], col="red", lty=2)

## -----------------------------------------------------------------------------
filtered.data <- normOffsets(filtered.data)
offsets <- assay(filtered.data, "offset")
head(offsets)

## ----normplot, fig.cap="Effect of non-linear normalization on the trended bias between two H3K9ac samples. Normalized log-fold changes are shown for all windows retained after filtering. A smoothed spline fitted to the log-fold change against the average abundance is also shown in red."----
norm.adjc <- calculateCPM(filtered.data, use.offsets=TRUE)
norm.fc <- norm.adjc[,4]-norm.adjc[,1]
smoothScatter(win.ab, norm.fc, ylim=c(-6, 6), xlim=c(0, 5),
    xlab="Average abundance", ylab="Log-fold change")

lfit <- smooth.spline(norm.fc~win.ab, df=5)
lines(win.ab[o], fitted(lfit)[o], col="red", lty=2)

## -----------------------------------------------------------------------------
celltype <- acdata$Description
celltype[grep("pro", celltype)] <- "proB"
celltype[grep("mature", celltype)] <- "matureB"

celltype <- factor(celltype)
design <- model.matrix(~0+celltype)
colnames(design) <- levels(celltype)
design

## -----------------------------------------------------------------------------
library(edgeR)
y <- asDGEList(filtered.data)
str(y)
y <- estimateDisp(y, design)
summary(y$trended.dispersion)

## ----bcvplot, fig.cap="Abundance-dependent trend in the BCV for each window, represented by the blue line. Common (red) and tagwise estimates (black) are also shown."----
plotBCV(y)

## ----qlplot, fig.cap="Effect of EB shrinkage on the raw QL dispersion estimate for each window (black) towards the abundance-dependent trend (blue) to obtain squeezed estimates (red)."----
fit <- glmQLFit(y, design, robust=TRUE)
plotQLDisp(fit)

## -----------------------------------------------------------------------------
summary(fit$df.prior)

## ----mdsplot, fig.cap="MDS plot with two dimensions for all samples in the H3K9ac data set. Samples are labelled and coloured according to the cell type."----
plotMDS(norm.adjc, labels=celltype,
    col=c("red", "blue")[as.integer(celltype)])

## -----------------------------------------------------------------------------
contrast <- makeContrasts(proB-matureB, levels=design)
res <- glmQLFTest(fit, contrast=contrast)
head(res$table)

## -----------------------------------------------------------------------------
merged <- mergeResults(filtered.data, res$table, tol=100, 
    merge.args=list(max.width=5000))
merged$regions

## -----------------------------------------------------------------------------
tabcom <- merged$combined
tabcom

## -----------------------------------------------------------------------------
is.sig <- tabcom$FDR <= 0.05
summary(is.sig)

## -----------------------------------------------------------------------------
table(tabcom$direction[is.sig])

## -----------------------------------------------------------------------------
tabbest <- merged$best
tabbest

## -----------------------------------------------------------------------------
is.sig.pos <- (tabbest$rep.logFC > 0)[is.sig]
summary(is.sig.pos)

## -----------------------------------------------------------------------------
out.ranges <- merged$regions
mcols(out.ranges) <- DataFrame(tabcom,
    best.pos=mid(ranges(rowRanges(filtered.data[tabbest$rep.test]))),
    best.logFC=tabbest$rep.logFC)
saveRDS(file="h3k9ac_results.rds", out.ranges)

## -----------------------------------------------------------------------------
simplified <- out.ranges[is.sig]
simplified$score <- -10*log10(simplified$FDR)
export(con="h3k9ac_results.bed", object=simplified)

## ---- message=FALSE-----------------------------------------------------------
library(org.Mm.eg.db)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
anno <- detailRanges(out.ranges, orgdb=org.Mm.eg.db,
    txdb=TxDb.Mmusculus.UCSC.mm10.knownGene)
head(anno$overlap)

## -----------------------------------------------------------------------------
head(anno$left)
head(anno$right)

## -----------------------------------------------------------------------------
meta <- mcols(out.ranges)
mcols(out.ranges) <- data.frame(meta, anno)

## ---- message=FALSE-----------------------------------------------------------
library(ChIPpeakAnno)
data(TSS.mouse.GRCm38)
minimal <- out.ranges
elementMetadata(minimal) <- NULL
anno.regions <- annotatePeakInBatch(minimal, AnnotationData=TSS.mouse.GRCm38)
colnames(elementMetadata(anno.regions))

## -----------------------------------------------------------------------------
prom <- suppressWarnings(promoters(TxDb.Mmusculus.UCSC.mm10.knownGene,
    upstream=3000, downstream=1000, columns=c("tx_name", "gene_id")))
entrez.ids <- sapply(prom$gene_id, FUN=function(x) x[1]) # Using the first Entrez ID.
gene.name <- select(org.Mm.eg.db, keys=entrez.ids, keytype="ENTREZID", column="SYMBOL")
prom$gene_name <- gene.name$SYMBOL[match(entrez.ids, gene.name$ENTREZID)]
head(prom)

## -----------------------------------------------------------------------------
olap.out <- overlapResults(filtered.data, regions=prom, res$table)
olap.out
simple <- DataFrame(ID=prom$tx_name, Gene=prom$gene_name, olap.out$combined)
simple[!is.na(simple$PValue),]

## ---- message=FALSE-----------------------------------------------------------
library(Gviz)
gax <- GenomeAxisTrack(col="black", fontsize=15, size=2)
greg <- GeneRegionTrack(TxDb.Mmusculus.UCSC.mm10.knownGene, showId=TRUE,
    geneSymbol=TRUE, name="", background.title="transparent")
symbols <- unlist(mapIds(org.Mm.eg.db, gene(greg), "SYMBOL",
    "ENTREZID", multiVals = "first"))
symbol(greg) <- symbols[gene(greg)]

## -----------------------------------------------------------------------------
o <- order(out.ranges$PValue)
sorted.ranges <- out.ranges[o]
sorted.ranges

## -----------------------------------------------------------------------------
cur.region <- sorted.ranges[1]
cur.region

## ---- echo=FALSE, results="hide"----------------------------------------------
if (cur.region!=GRanges("chr17", IRanges(34285101, 34290050))) {
    stop("first region does not match expectations")
}

## ----simplebroadplot, fig.width=8, fig.asp=0.75, fig.cap="Coverage tracks for a simple DB event between pro-B and mature B cells, across a broad region in the H3K9ac data set. Read coverage for each sample is shown as a per-million value at each base."----
collected <- list()
lib.sizes <- filtered.data$totals/1e6
for (i in seq_along(acdata$Path)) {
    reads <- extractReads(bam.file=acdata$Path[[i]], cur.region, param=param)
    cov <- as(coverage(reads)/lib.sizes[i], "GRanges")
    collected[[i]] <- DataTrack(cov, type="histogram", lwd=0, ylim=c(0,10),
        name=acdata$Description[i], col.axis="black", col.title="black",
        fill="darkgray", col.histogram=NA)
}
plotTracks(c(gax, collected, greg), chromosome=as.character(seqnames(cur.region)),
    from=start(cur.region), to=end(cur.region))

## -----------------------------------------------------------------------------
complex <- sorted.ranges$num.up.logFC > 0 & sorted.ranges$num.down.logFC > 0
cur.region <- sorted.ranges[complex][2]
cur.region

## ---- echo=FALSE, results="hide"----------------------------------------------
if (cur.region!=GRanges("chr5", IRanges(122987201, 122991450))) {
    stop("second region does not match expectations")
}

## ----complexplot, fig.width=8, fig.asp=0.75, fig.cap="Coverage tracks for a complex DB event in the H3K9ac data set, shown as per-million values."----
collected <- list()
for (i in seq_along(acdata$Path)) {
    reads <- extractReads(bam.file=acdata$Path[[i]], cur.region, param=param)
    cov <- as(coverage(reads)/lib.sizes[i], "GRanges")
    collected[[i]] <- DataTrack(cov, type="histogram", lwd=0, ylim=c(0,3),
        name=acdata$Description[i], col.axis="black", col.title="black",
        fill="darkgray", col.histogram=NA)
}
plotTracks(c(gax, collected, greg), chromosome=as.character(seqnames(cur.region)),
    from=start(cur.region), to=end(cur.region))

## -----------------------------------------------------------------------------
sharp <- sorted.ranges$num.tests < 20
cur.region <- sorted.ranges[sharp][1]
cur.region

## ---- echo=FALSE, results="hide"----------------------------------------------
if (cur.region!=GRanges("chr16", IRanges(36665551, 36666200))) {
    stop("second region does not match expectations")
}

## ----simplesharpplot, fig.width=8, fig.asp=0.75, fig.cap="Coverage tracks for a sharp and simple DB event in the H3K9ac data set, shown as per-million values."----
collected <- list()
for (i in seq_along(acdata$Path)) {
    reads <- extractReads(bam.file=acdata$Path[[i]], cur.region, param=param)
    cov <- as(coverage(reads)/lib.sizes[i], "GRanges")
    collected[[i]] <- DataTrack(cov, type="histogram", lwd=0, ylim=c(0,3),
        name=acdata$Description[i], col.axis="black", col.title="black",
        fill="darkgray", col.histogram=NA)
}
plotTracks(c(gax, collected, greg), chromosome=as.character(seqnames(cur.region)),
    from=start(cur.region), to=end(cur.region))

## -----------------------------------------------------------------------------
sessionInfo()


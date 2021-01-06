## ----style, echo=FALSE, results='hide', message=FALSE-------------------------
library(BiocStyle)
library(knitr)
opts_chunk$set(error=FALSE, message=FALSE, warning=FALSE, cache=TRUE)

## -----------------------------------------------------------------------------
library(chipseqDBData)
h3k27me3data <- H3K27me3Data()
h3k27me3data

## -----------------------------------------------------------------------------
library(Rsamtools)
diagnostics <- list()
for (b in seq_along(h3k27me3data$Path)) {
    bam <- h3k27me3data$Path[[b]]
    total <- countBam(bam)$records
    mapped <- countBam(bam, param=ScanBamParam(
        flag=scanBamFlag(isUnmapped=FALSE)))$records
    marked <- countBam(bam, param=ScanBamParam(
        flag=scanBamFlag(isUnmapped=FALSE, isDuplicate=TRUE)))$records
    diagnostics[[b]] <- c(Total=total, Mapped=mapped, Marked=marked)
}

diag.stats <- data.frame(do.call(rbind, diagnostics))
rownames(diag.stats) <- h3k27me3data$Name
diag.stats$Prop.mapped <- diag.stats$Mapped/diag.stats$Total*100
diag.stats$Prop.marked <- diag.stats$Marked/diag.stats$Mapped*100
diag.stats

## -----------------------------------------------------------------------------
library(BiocFileCache)
bfc <- BiocFileCache("local", ask=FALSE)
black.path <- bfcrpath(bfc, file.path("http://hgdownload.cse.ucsc.edu",
    "goldenPath/mm10/bigZips/chromOut.tar.gz"))
tmpdir <- tempfile()
dir.create(tmpdir)
untar(black.path, exdir=tmpdir)

# Iterate through all chromosomes.
collected <- list()
for (x in list.files(tmpdir, full=TRUE)) {
    f <- list.files(x, full=TRUE, pattern=".fa.out")
    to.get <- vector("list", 15)
    to.get[[5]] <- "character"
    to.get[6:7] <- "integer"
    collected[[length(collected)+1]] <- read.table(f, skip=3, 
        stringsAsFactors=FALSE, colClasses=to.get)
}

collected <- do.call(rbind, collected)
blacklist <- GRanges(collected[,1], IRanges(collected[,2], collected[,3]))
blacklist

## -----------------------------------------------------------------------------
library(csaw)
param <- readParam(minq=10, discard=blacklist)
param

## -----------------------------------------------------------------------------
win.data <- windowCounts(h3k27me3data$Path, param=param, width=2000,
    spacing=500, ext=200)
win.data

## -----------------------------------------------------------------------------
bins <- windowCounts(h3k27me3data$Path, bin=TRUE, width=10000, param=param)
win.data <- normFactors(bins, se.out=win.data)
(normfacs <- win.data$norm.factors)

## ----compoplot, fig.width=12, fig.asp=0.5, fig.cap="Mean-difference plots for the bin counts, comparing sample 1 to all other samples. The red line represents the log-ratio of the normalization factors between samples."----
bin.ab <- scaledAverage(bins)
adjc <- calculateCPM(bins, use.norm.factors=FALSE)

par(cex.lab=1.5, mfrow=c(1,3))
smoothScatter(bin.ab, adjc[,1]-adjc[,2], ylim=c(-6, 6),
    xlab="Average abundance", ylab="Log-ratio (1 vs 2)")
abline(h=log2(normfacs[1]/normfacs[4]), col="red")

smoothScatter(bin.ab, adjc[,1]-adjc[,3], ylim=c(-6, 6),
    xlab="Average abundance", ylab="Log-ratio (1 vs 3)")
abline(h=log2(normfacs[2]/normfacs[4]), col="red")

smoothScatter(bin.ab, adjc[,1]-adjc[,4], ylim=c(-6, 6),
    xlab="Average abundance", ylab="Log-ratio (1 vs 4)")
abline(h=log2(normfacs[3]/normfacs[4]), col="red")

## -----------------------------------------------------------------------------
filter.stat <- filterWindowsGlobal(win.data, bins)
min.fc <- 2 

## ----bghistplot, fig.cap="Histogram of average abundances across all 10 kbp genomic bins. The filter threshold is shown as the red line."----
hist(filter.stat$back.abundances, main="", breaks=50,
    xlab="Background abundance (log2-CPM)")
threshold <- filter.stat$abundances[1] - filter.stat$filter[1] + log2(min.fc)
abline(v=threshold, col="red")

## -----------------------------------------------------------------------------
keep <- filter.stat$filter > log2(min.fc)
summary(keep)
filtered.data <- win.data[keep,]

## -----------------------------------------------------------------------------
library(edgeR)
y <- asDGEList(filtered.data)
str(y)

## -----------------------------------------------------------------------------
genotype <- h3k27me3data$Description
genotype[grep("control", genotype)] <- "wt"
genotype[grep("knock-out", genotype)] <- "ko"

genotype <- factor(genotype)
design <- model.matrix(~0+genotype)
colnames(design) <- levels(genotype)
design

## ----bcvplot, fig.cap="Abundance-dependent trend in the BCV for each window, represented by the blue line. Common (red) and tagwise estimates (black) are also shown."----
y <- estimateDisp(y, design)
summary(y$trended.dispersion)
plotBCV(y)

## ----qlplot, fig.cap="Effect of EB shrinkage on the raw QL dispersion estimate for each window (black) towards the abundance-dependent trend (blue) to obtain squeezed estimates (red)."----
fit <- glmQLFit(y, design, robust=TRUE)
summary(fit$df.prior)
plotQLDisp(fit)

## ----mdsplot, fig.cap="MDS plot with two dimensions for all samples in the H3K27me3 data set. Samples are labelled and coloured according to the genotype. A larger top set of windows was used to improve the visualization of the genome-wide differences between the WT samples."----
plotMDS(cpm(y, log=TRUE), top=10000, labels=genotype,
    col=c("red", "blue")[as.integer(genotype)])

## -----------------------------------------------------------------------------
contrast <- makeContrasts(wt-ko, levels=design)
res <- glmQLFTest(fit, contrast=contrast)

## -----------------------------------------------------------------------------
# Counting into 500 bp windows.
win.data2 <- windowCounts(h3k27me3data$Path, param=param, width=500,
    spacing=100, ext=200)

# Re-using the same normalization factors.
win.data2$norm.factors <- win.data$norm.factors

# Filtering on abundance.
filter.stat2 <- filterWindowsGlobal(win.data2, bins)
keep2 <- filter.stat2$filter > log2(min.fc)
filtered.data2 <- win.data2[keep2,]

# Performing the statistical analysis.
y2 <- asDGEList(filtered.data2)
y2 <- estimateDisp(y2, design)
fit2 <- glmQLFit(y2, design, robust=TRUE)
res2 <- glmQLFTest(fit2, contrast=contrast)

## -----------------------------------------------------------------------------
merged <- mergeResultsList(list(filtered.data, filtered.data2), 
    tab.list=list(res$table, res2$table),
    equiweight=TRUE, tol=100, merge.args=list(max.width=30000))
merged$regions

## -----------------------------------------------------------------------------
tabcom <- merged$combined
is.sig <- tabcom$FDR <= 0.05
summary(is.sig)

## -----------------------------------------------------------------------------
table(tabcom$direction[is.sig])

## -----------------------------------------------------------------------------
tabbest <- merged$best
is.sig.pos <- (tabbest$rep.logFC > 0)[is.sig]
summary(is.sig.pos)

## -----------------------------------------------------------------------------
out.ranges <- merged$regions
mcols(out.ranges) <- data.frame(tabcom, best.logFC=tabbest$rep.logFC)
saveRDS(file="h3k27me3_results.rds", out.ranges)

## -----------------------------------------------------------------------------
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
anno <- detailRanges(out.ranges, orgdb=org.Mm.eg.db,
    txdb=TxDb.Mmusculus.UCSC.mm10.knownGene)
mcols(out.ranges) <- cbind(mcols(out.ranges), anno)

## -----------------------------------------------------------------------------
cdx2 <- genes(TxDb.Mmusculus.UCSC.mm10.knownGene)["12591"] # Cdx2 Entrez ID
cur.region <- subsetByOverlaps(out.ranges, cdx2)[1]
cur.region

## ---- results="hide", echo=FALSE----------------------------------------------
if (cur.region$FDR > 0.05 || cur.region$best.logFC < 0) {
    stop("Cdx2 should be significantly upregulated in WT")
}

## -----------------------------------------------------------------------------
library(Gviz)
gax <- GenomeAxisTrack(col="black", fontsize=15, size=2)
greg <- GeneRegionTrack(TxDb.Mmusculus.UCSC.mm10.knownGene, showId=TRUE,
    geneSymbol=TRUE, name="", background.title="transparent")
symbols <- unlist(mapIds(org.Mm.eg.db, gene(greg), "SYMBOL",
    "ENTREZID", multiVals = "first"))
symbol(greg) <- symbols[gene(greg)]

## ----tfplot, fig.width=8, fig.asp=0.75, fig.cap="Coverage tracks for a region with H3K27me3 enrichment in KO (top two tracks) against the WT (last two tracks)."----
collected <- list()
lib.sizes <- filtered.data$totals/1e6
for (i in seq_along(h3k27me3data$Path)) {
    reads <- extractReads(bam.file=h3k27me3data$Path[[i]], cur.region, param=param)
    cov <- as(coverage(reads)/lib.sizes[i], "GRanges")
    collected[[i]] <- DataTrack(cov, type="histogram", lwd=0, ylim=c(0,1),
        name=h3k27me3data$Description[i], col.axis="black", col.title="black",
        fill="darkgray", col.histogram=NA)
}

plotTracks(c(gax, collected, greg), chromosome=as.character(seqnames(cur.region)),
    from=start(cur.region), to=end(cur.region))

## -----------------------------------------------------------------------------
col1a2 <- genes(TxDb.Mmusculus.UCSC.mm10.knownGene)["12843"] # Col1a2 Entrez ID
cur.region <- subsetByOverlaps(out.ranges, col1a2)
cur.region

## ---- results="hide", echo=FALSE----------------------------------------------
if (cur.region$FDR < 0.05) {
    stop("Col1a2 should not be significantly DB")
}

## -----------------------------------------------------------------------------
sessionInfo()


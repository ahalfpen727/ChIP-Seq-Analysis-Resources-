## ----style, echo=FALSE, results='hide', message=FALSE-------------------------
library(BiocStyle)
library(knitr)
opts_chunk$set(error=FALSE, message=FALSE, warning=FALSE, cache=TRUE)

## -----------------------------------------------------------------------------
library(chipseqDBData)
cbpdata <- CBPData()
cbpdata

## -----------------------------------------------------------------------------
library(Rsamtools)
diagnostics <- list()
for (b in seq_along(cbpdata$Path)) {
    bam <- cbpdata$Path[[b]]
    total <- countBam(bam)$records
    mapped <- countBam(bam, param=ScanBamParam(
        flag=scanBamFlag(isUnmapped=FALSE)))$records
    marked <- countBam(bam, param=ScanBamParam(
        flag=scanBamFlag(isUnmapped=FALSE, isDuplicate=TRUE)))$records
    diagnostics[[b]] <- c(Total=total, Mapped=mapped, Marked=marked)
}

diag.stats <- data.frame(do.call(rbind, diagnostics))
rownames(diag.stats) <- cbpdata$Name
diag.stats$Prop.mapped <- diag.stats$Mapped/diag.stats$Total*100
diag.stats$Prop.marked <- diag.stats$Marked/diag.stats$Mapped*100
diag.stats

## -----------------------------------------------------------------------------
library(BiocFileCache)
bfc <- BiocFileCache("local", ask=FALSE)
black.path <- bfcrpath(bfc, file.path("https://www.encodeproject.org",
    "files/ENCFF547MET/@@download/ENCFF547MET.bed.gz"))

library(rtracklayer)
blacklist <- import(black.path)

## -----------------------------------------------------------------------------
library(csaw)
param <- readParam(minq=10, discard=blacklist)
param

## -----------------------------------------------------------------------------
x <- correlateReads(cbpdata$Path, param=reform(param, dedup=TRUE))
frag.len <- maximizeCcf(x)
frag.len

## ----ccfplot, fig.cap="Cross-correlation function (CCF) against delay distance for the CBP dataset. The delay with the maximum correlation is shown as the red line."----
plot(1:length(x)-1, x, xlab="Delay (bp)", ylab="CCF", type="l")
abline(v=frag.len, col="red")
text(x=frag.len, y=min(x), paste(frag.len, "bp"), pos=4, col="red")

## -----------------------------------------------------------------------------
win.data <- windowCounts(cbpdata$Path, param=param, width=10, ext=frag.len)
win.data

## -----------------------------------------------------------------------------
bins <- windowCounts(cbpdata$Path, bin=TRUE, width=10000, param=param)
win.data <- normFactors(bins, se.out=win.data)
(normfacs <- win.data$norm.factors)

## ----compoplot, fig.width=12, fig.asp=0.5, fig.cap="Mean-difference plots for the bin counts, comparing sample 4 to all other samples. The red line represents the log-ratio of the normalization factors between samples."----
bin.ab <- scaledAverage(bins)
adjc <- calculateCPM(bins, use.norm.factors=FALSE)

par(cex.lab=1.5, mfrow=c(1,3))
smoothScatter(bin.ab, adjc[,1]-adjc[,4], ylim=c(-6, 6),
    xlab="Average abundance", ylab="Log-ratio (1 vs 4)")
abline(h=log2(normfacs[1]/normfacs[4]), col="red")

smoothScatter(bin.ab, adjc[,2]-adjc[,4], ylim=c(-6, 6),
    xlab="Average abundance", ylab="Log-ratio (2 vs 4)")
abline(h=log2(normfacs[2]/normfacs[4]), col="red")

smoothScatter(bin.ab, adjc[,3]-adjc[,4], ylim=c(-6, 6),
    xlab="Average abundance", ylab="Log-ratio (3 vs 4)")
abline(h=log2(normfacs[3]/normfacs[4]), col="red")

## -----------------------------------------------------------------------------
filter.stat <- filterWindowsGlobal(win.data, bins)
min.fc <- 3
keep <- filter.stat$filter > log2(min.fc)
summary(keep)
filtered.data <- win.data[keep,]

## -----------------------------------------------------------------------------
library(edgeR)
y <- asDGEList(filtered.data)
summary(y)

## -----------------------------------------------------------------------------
genotype <- cbpdata$Description
genotype[grep("wild-type", genotype)] <- "wt"
genotype[grep("knock-out", genotype)] <- "ko"

genotype <- factor(genotype)
design <- model.matrix(~0+genotype)
colnames(design) <- levels(genotype)
design

## ----bcvplot, fig.cap="Abundance-dependent trend in the biological coefficient of variation (i.e., the root-NB dispersion) for each window, represented by the blue line. Common (red) and tagwise estimates (black) are also shown."----
y <- estimateDisp(y, design)
summary(y$trended.dispersion)
plotBCV(y)

## ----qlplot, fig.cap="Effect of EB shrinkage on the raw QL dispersion estimate for each window (black) towards the abundance-dependent trend (blue) to obtain squeezed estimates (red). Quarter-root estimates are shown for greater dynamic range."----
fit <- glmQLFit(y, design, robust=TRUE)
summary(fit$df.prior)
plotQLDisp(fit)

## ----mdsplot, fig.cap="MDS plot with two dimensions for all samples in the CBP dataset. Samples are labelled and coloured according to the genotype. A larger top set of windows was used to improve the visualization of the genome-wide differences between the WT samples."----
plotMDS(cpm(y, log=TRUE), top=10000, labels=genotype,
    col=c("red", "blue")[as.integer(genotype)])

## -----------------------------------------------------------------------------
contrast <- makeContrasts(wt-ko, levels=design)
res <- glmQLFTest(fit, contrast=contrast)

## -----------------------------------------------------------------------------
merged <- mergeResults(filtered.data, res$table, tol=100, 
    merge.args=list(max.width=5000))
merged$regions
tabcom <- merged$combined
is.sig <- tabcom$FDR <= 0.05
summary(is.sig)

## -----------------------------------------------------------------------------
table(tabcom$direction[is.sig])

# Direction according the best window in each cluster.
tabbest <- merged$best
is.sig.pos <- (tabbest$rep.logFC > 0)[is.sig]
summary(is.sig.pos)

## -----------------------------------------------------------------------------
out.ranges <- merged$regions
mcols(out.ranges) <- DataFrame(tabcom,
    best.pos=mid(ranges(rowRanges(filtered.data[tabbest$rep.test]))),
    best.logFC=tabbest$rep.logFC)
saveRDS(file="cbp_results.rds", out.ranges)

## -----------------------------------------------------------------------------
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
anno <- detailRanges(out.ranges, orgdb=org.Mm.eg.db,
    txdb=TxDb.Mmusculus.UCSC.mm10.knownGene)
mcols(out.ranges) <- cbind(mcols(out.ranges), anno)

## -----------------------------------------------------------------------------
o <- order(out.ranges$PValue)    
cur.region <- out.ranges[o[2]]
cur.region

## ---- results="hide", echo=FALSE----------------------------------------------
if (!overlapsAny(cur.region, GRanges("chr16", IRanges(70313851, 70314860)), type="equal")) {
        warning("first region does not match expectations")
}

## -----------------------------------------------------------------------------
library(Gviz)
gax <- GenomeAxisTrack(col="black", fontsize=15, size=2)
greg <- GeneRegionTrack(TxDb.Mmusculus.UCSC.mm10.knownGene, showId=TRUE,
    geneSymbol=TRUE, name="", background.title="transparent")
symbols <- unlist(mapIds(org.Mm.eg.db, gene(greg), "SYMBOL",
    "ENTREZID", multiVals = "first"))
symbol(greg) <- symbols[gene(greg)]

## ----tfplot, fig.width=8, fig.asp=0.75, fig.cap="Coverage tracks for TF binding sites that are differentially bound in the WT (top two tracks) against the KO (last two tracks). Blue and red tracks represent forward- and reverse-strand coverage, respectively, on a per-million scale (capped at 5 in SRR1145788, for visibility)."----
collected <- list()
lib.sizes <- filtered.data$totals/1e6

for (i in seq_along(cbpdata$Path)) {
    reads <- extractReads(bam.file=cbpdata$Path[[i]], cur.region, param=param)
    pcov <- as(coverage(reads[strand(reads)=="+"])/lib.sizes[i], "GRanges")
    ncov <- as(coverage(reads[strand(reads)=="-"])/-lib.sizes[i], "GRanges")
    ptrack <- DataTrack(pcov, type="histogram", lwd=0, ylim=c(-5, 5),
        name=cbpdata$Description[i], col.axis="black", col.title="black",
        fill="blue", col.histogram=NA)
    ntrack <- DataTrack(ncov, type="histogram", lwd=0, ylim=c(-5, 5),
        fill="red", col.histogram=NA)
    collected[[i]] <- OverlayTrack(trackList=list(ptrack, ntrack))
}

plotTracks(c(gax, collected, greg), chromosome=as.character(seqnames(cur.region)),
    from=start(cur.region), to=end(cur.region))

## -----------------------------------------------------------------------------
sessionInfo()


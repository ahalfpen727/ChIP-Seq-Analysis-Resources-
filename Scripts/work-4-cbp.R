## ----style, echo=FALSE, results='hide', message=FALSE----------------------
library(BiocStyle)
library(knitr)
opts_chunk$set(error=FALSE, message=FALSE, warning=FALSE)
opts_chunk$set(fig.asp=1)

## --------------------------------------------------------------------------
sra.numbers <- c("SRR1145787", "SRR1145788", "SRR1145789", "SRR1145790")
genotype <- c("wt", "wt", "ko", "ko")
all.sra <- paste0(sra.numbers, ".sra")
data.frame(SRA=all.sra, Condition=genotype)

## ---- echo=FALSE, results="hide"-------------------------------------------
remap <- FALSE
redownload <- any(!file.exists(paste0(sra.numbers, ".bam")))

## ---- eval=remap-----------------------------------------------------------
#  for (sra in all.sra) {
#      code <- system(paste("fastq-dump", sra))
#      stopifnot(code==0L)
#  }
#  all.fastq <- paste0(sra.numbers, ".fastq")

## ---- eval=remap, message=FALSE, results="hide"----------------------------
#  bam.files <- paste0(sra.numbers, ".bam")
#  align(index="index/mm10", readfile1=all.fastq, type=1, phredOffset=64,
#      input_format="FASTQ", output_file=bam.files)

## ---- eval=remap, results="hide"-------------------------------------------
#  temp.bam <- "cbp_temp.bam"
#  temp.file <- "cbp_metric.txt"
#  temp.dir <- "cbp_working"
#  dir.create(temp.dir)
#  for (bam in bam.files) {
#      out <- suppressWarnings(sortBam(bam, "cbp_temp"))
#      file.rename(out, bam)
#      code <- system(sprintf("MarkDuplicates I=%s O=%s M=%s \\
#          TMP_DIR=%s AS=true REMOVE_DUPLICATES=false \\
#          VALIDATION_STRINGENCY=SILENT",
#          bam, temp.bam, temp.file, temp.dir))
#      stopifnot(code==0L)
#      file.rename(temp.bam, bam)
#  }
#  indexBam(bam.files)

## ---- eval=remap, echo = FALSE, results = 'hide'---------------------------
#  # Cleaning up
#  unlink(all.fastq)
#  unlink(temp.dir, recursive=TRUE)
#  unlink(temp.file)

## ---- echo=FALSE, results="hide", eval=!remap, message=FALSE---------------
bam.files <- paste0(sra.numbers, ".bam")

## ---- echo=FALSE, results="hide", eval=redownload, message=FALSE-----------
core.loc <- "http://s3.amazonaws.com/chipseqdb-bamfiles/"
for (bam in bam.files) {
    bam.url <- paste0(core.loc, bam)
    download.file(bam.url, bam)
    download.file(paste0(bam.url, ".bai"), paste0(bam, ".bai"))
}

## --------------------------------------------------------------------------
library(Rsamtools)
diagnostics <- list()
for (bam in bam.files) {
    total <- countBam(bam)$records
    mapped <- countBam(bam, param=ScanBamParam(
        flag=scanBamFlag(isUnmapped=FALSE)))$records
    marked <- countBam(bam, param=ScanBamParam(
        flag=scanBamFlag(isUnmapped=FALSE, isDuplicate=TRUE)))$records
    diagnostics[[bam]] <- c(Total=total, Mapped=mapped, Marked=marked)
}
diag.stats <- data.frame(do.call(rbind, diagnostics))
diag.stats$Prop.mapped <- diag.stats$Mapped/diag.stats$Total*100
diag.stats$Prop.marked <- diag.stats$Marked/diag.stats$Mapped*100
diag.stats

## --------------------------------------------------------------------------
library(rtracklayer)
blacklist <- import("mm10.blacklist.bed.gz")
library(csaw)
param <- readParam(minq=50, discard=blacklist)

## --------------------------------------------------------------------------
x <- correlateReads(bam.files, param=reform(param, dedup=TRUE))
frag.len <- maximizeCcf(x)
frag.len

## --------------------------------------------------------------------------
win.data <- windowCounts(bam.files, param=param, width=10, ext=frag.len)
win.data

## --------------------------------------------------------------------------
bins <- windowCounts(bam.files, bin=TRUE, width=10000, param=param)
win.data <- normFactors(bins, se.out=win.data)
normfacs <- win.data$norm.factors
normfacs

## ----compoplot, fig.width=12, fig.asp=0.5, fig.cap="Mean-difference plots for the bin counts, comparing library 4 to all other libraries. The red line represents the log-ratio of the normalization factors between libraries."----
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

## --------------------------------------------------------------------------
filter.stat <- filterWindows(win.data, bins, type="global")
min.fc <- 3
keep <- filter.stat$filter > log2(min.fc)
summary(keep)
filtered.data <- win.data[keep,]

## --------------------------------------------------------------------------
genotype <- factor(genotype)
design <- model.matrix(~0+genotype)
colnames(design) <- levels(genotype)
design

## --------------------------------------------------------------------------
library(edgeR)
y <- asDGEList(filtered.data, norm.factors=normfacs)
y <- estimateDisp(y, design)
summary(y$trended.dispersion)
fit <- glmQLFit(y, design, robust=TRUE)
summary(fit$df.prior)

## ----mdsplot, fig.cap="MDS plot with two dimensions for all libraries in the CBP data set. Libraries are labelled and coloured according to the genotype. A larger top set of windows was used to improve the visualization of the genome-wide differences between the WT libraries."----
plotMDS(cpm(y, log=TRUE), top=10000, labels=genotype,
    col=c("red", "blue")[as.integer(genotype)])

## --------------------------------------------------------------------------
contrast <- makeContrasts(wt-ko, levels=design)
res <- glmQLFTest(fit, contrast=contrast)
merged <- mergeWindows(rowRanges(filtered.data), tol=100, max.width=5000)
tabcom <- combineTests(merged$id, res$table)
tabbest <- getBestTest(merged$id, res$table)
is.sig <- tabcom$FDR <= 0.05
summary(is.sig)
table(tabcom$direction[is.sig])
is.sig.pos <- (tabbest$logFC > 0)[is.sig]
summary(is.sig.pos)

## --------------------------------------------------------------------------
out.ranges <- merged$region
elementMetadata(out.ranges) <- data.frame(tabcom,
    best.pos=mid(ranges(rowRanges(filtered.data[tabbest$best]))),
    best.logFC=tabbest$logFC)
saveRDS(file="cbp_results.rds", out.ranges)
save(file="cbp_objects.Rda", win.data, bins, y)

## --------------------------------------------------------------------------
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
anno <- detailRanges(out.ranges, orgdb=org.Mm.eg.db,
            txdb=TxDb.Mmusculus.UCSC.mm10.knownGene)
meta <- elementMetadata(out.ranges)
elementMetadata(out.ranges) <- data.frame(meta, anno)

## --------------------------------------------------------------------------
o <- order(out.ranges$PValue)    
cur.region <- out.ranges[o[2]]
cur.region

## ---- results="hide", echo=FALSE-------------------------------------------
if (!overlapsAny(cur.region, GRanges("chr16", IRanges(70313851, 70314860)), type="equal")) {
        warning("first region does not match expectations")
}

## ----tfplot, fig.width=8, fig.asp=0.75, fig.cap="Coverage tracks for TF binding sites that are differentially bound in the WT (top two tracks) against the KO (last two tracks). Blue and red tracks represent forward- and reverse-strand coverage, respectively, on a per-million scale (capped at 5 in SRR1145788, for visibility)."----
library(Gviz)
collected <- list()
lib.sizes <- filtered.data$totals/1e6

for (i in seq_along(bam.files)) {
    reads <- extractReads(bam.file=bam.files[i], cur.region, param=param)
    pcov <- as(coverage(reads[strand(reads)=="+"])/lib.sizes[i], "GRanges")
    ncov <- as(coverage(reads[strand(reads)=="-"])/-lib.sizes[i], "GRanges")
    ptrack <- DataTrack(pcov, type="histogram", lwd=0, ylim=c(-5, 5),
        name=bam.files[i], col.axis="black", col.title="black",
        fill="blue", col.histogram=NA)
    ntrack <- DataTrack(ncov, type="histogram", lwd=0, ylim=c(-5, 5),
        fill="red", col.histogram=NA)
    collected[[i]] <- OverlayTrack(trackList=list(ptrack, ntrack))
}

gax <- GenomeAxisTrack(col="black", fontsize=15, size=2)
greg <- GeneRegionTrack(TxDb.Mmusculus.UCSC.mm10.knownGene, showId=TRUE,
    geneSymbol=TRUE, name="", background.title="transparent")
plotTracks(c(gax, collected, greg), chromosome=as.character(seqnames(cur.region)),
    from=start(cur.region), to=end(cur.region))

## --------------------------------------------------------------------------
sessionInfo()


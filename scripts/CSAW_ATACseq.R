# library(chipseqDBData)
library(csaw)
library(edgeR)
# tf.data <- NFYAData()
# tf.data <- head(tf.data, -1) # skip the input.
# bam.files <- tf.data$Path

# cell.type <- sub("NF-YA ([^ ]+) .*", "\\1", tf.data$Description)
# design <- model.matrix(~factor(cell.type))
# colnames(design) <- c("intercept", "cell.type")

makeQCplots_chip_PE <- function(bam.file, outplot, pe.param){

    ## Histogram to check frag size cutoff
    message("Checking fragment sizes")
    fragsize <- csaw::getPESizes(bam.file)

    ## Checking cross-correlation
    message("Checking strand cross-correlation")
    max.delay <- 500
    dedup.on <- csaw::readParam(dedup = TRUE, minq = 20)
    CCF <- csaw::correlateReads(bam.file, max.delay, param = dedup.on)

    ## Choosing appropriate window size
    message("Checking read distribution around putative peaks")
    plotwc <- function(curbam){
        windowed <- csaw::windowCounts(curbam, spacing = 50, param = pe.param, filter = 20)
        rwsms <- rowSums(SummarizedExperiment::assay(windowed))
        maxed <- csaw::findMaxima(SummarizedExperiment::rowRanges(windowed), range = 1000, metric = rwsms)
        curbam.out <- csaw::profileSites(curbam, SummarizedExperiment::rowRanges(windowed)[maxed], param = pe.param)
        return(curbam.out)
    }
    collected <- plotwc(bam.file)
    xranged <- as.integer(names(collected))

    ## plot
    message("Plotting")
    pdf(outplot)
    # fragment sizes
    hist(log10(fragsize$sizes),
         breaks = 50,
         xlab= " log10(Fragment sizes)",
         ylab = "Frequency",
         main = "fragment sizes",
         col = "steelblue")
    abline(v = 400, col = "red",lwd = 3)

    # cross correlation
    plot(0:max.delay, CCF, type = "l", ylab = "CCF", xlab = "Delay (bp)", main = "PE-Cross-correlation")

    # coverage in windows
    plot(xranged, collected, type = "l", col = "blue", xlim = c(-1000, 1000), lwd = 2,
         xlab = "Distance (bp)", ylab = "Relative coverage per base")
    abline(v = c(-150,200), col = "dodgerblue", lty = 2)
    legend("topright", col = "dodgerblue", legend = "specified window size")

    dev.off()

}

tmmNormalize_chip <- function(chipCountObject, binsize, plotfile){


    bam.files <- SummarizedExperiment::colData(chipCountObject$windowCounts)$bam.files
    # Get norm factors
    wider <- csaw::windowCounts(bam.files, bin = TRUE, width = binsize, param = chipCountObject$pe.param)
    normfacs <- csaw::normFactors(wider, se.out=FALSE)

    chipCountObject$normFactors <- normfacs

    # get norm counts
    adj.counts <- edgeR::cpm(csaw::asDGEList(wider), log = TRUE)
    chipCountObject$background_logcpm <- adj.counts

    # plot normalized counts
    pdf(plotfile, width=20, height=20)
    par(mfrow = c(2, 2),mar = c(5, 4, 2, 1.5))
    #par(mfrow = c(3, 3), mar = c(5, 4, 2, 1.5))
    #for (i in 1:(length(bam.files) - 1)) {
    #    cur.x <- adj.counts[, 1]
    #    cur.y <- adj.counts[, 1 + i]
    #    smoothScatter(x = (cur.x + cur.y)/2 + 6*log2(10),
    #                y = cur.x-cur.y, xlab = "A",
    #                ylab = "M",
    #                main = paste("1 vs", i + 1))

    #    all.dist <- diff(log2(normfacs[c(i + 1, 1)]))
    #    abline(h = all.dist, col = "red")
    #}
    ## MDS plot to check for replicate variability
    for (top in c(100, 500, 1000, 5000)) {
        limma::plotMDS(adj.counts, main = top,
                   col = as.numeric(chipCountObject$sampleSheet$condition),
                   labels = chipCountObject$sampleSheet$name, top = top)
    }
    dev.off()

    ## Return normfactors
    return(chipCountObject)
}

getDBregions_chip <- function(chipCountObject, plotfile = NULL){

    # Make DGElist
    y <- csaw::asDGEList(chipCountObject$windowCounts, norm.factors = chipCountObject$normFactors)
    if(chipCountObject$designType=="condition"){
    colnames(y)<-chipCountObject$sampleSheet$name}else{
    colnames(y)<-paste0(rep(chipCountObject$sampleSheet$name,each=2),".genome",c(1,2))}
    design <- chipCountObject$design
    # Estimate dispersions
    y <- edgeR::estimateDisp(y, design)
    o <- order(y$AveLogCPM)
    fit <- edgeR::glmQLFit(y, design, robust=TRUE)

    # and plot dispersions
    if(!(is.null(plotfile))){
        pdf(plotfile)
        par(mfrow = c(1,2))
        plot(y$AveLogCPM[o], sqrt(y$trended.dispersion[o]), type = "l", lwd = 2,
             ylim = c(0, 1), xlab = expression("Ave."~Log[2]~"CPM"),
             ylab = ("Biological coefficient of variation"))
        edgeR::plotQLDisp(fit)
        dev.off()
    }

    ### TEST for DB windows
    # check design type
    if(chipCountObject$designType != "condition") {
        results <- edgeR::glmQLFTest(fit, coef = paste0("allelegenome2"))
    } else {
        results <- edgeR::glmQLFTest(fit, coef = paste0("condition",unique(chipCountObject$sampleSheet$condition)[2]))
    }

    # Merge DB windows into regions: Using quick and dirty method
    merged <- csaw::mergeWindows(SummarizedExperiment::rowRanges(chipCountObject$windowCounts), tol = 100L)
    # get combined test p-value for merged windows
    tabcom <- csaw::combineTests(merged$id, results$table, pval.col = 4, fc.col = 1)
    # get fold change of the best window within each combined cluster
    tab.best <- csaw::getBestTest(merged$id, results$table,pval.col=4,cpm.col=1)
    tabcom$best.logFC <- tab.best$rep.logFC
    tabcom$best.test <- tab.best$rep.test
    tabcom$best.start <- GenomicRanges::start(SummarizedExperiment::rowRanges(chipCountObject$windowCounts))[tab.best$rep.test]

    # Return all results
    chipResultObject <- list(fit = fit, results = results, mergedRegions = merged, combinedPvalues = tabcom)
    return(chipResultObject)
}

writeOutput_chip <- function(chipResultObject, outfile_prefix, fdrcutoff,lfccutoff){
    # get merged regions
    merged <- chipResultObject$mergedRegions
    tabcom <- chipResultObject$combinedPvalues
    merged$region$score <- -10*log10(tabcom$FDR)
    names(merged$region) <- paste0("region", 1:length(merged$region))
    tabcom$name <- names(merged$region)
    ## export merged data
    rtracklayer::export.bed(merged$region, paste0(outfile_prefix, "_allregions.bed"))
    write.table(tabcom, file = paste0(outfile_prefix,"_scores.txt"),
            row.names = FALSE, quote = FALSE, sep = "\t")

    ## export FDR significant data
    is.sig <- tabcom$FDR <= fdrcutoff
    test <- merged$region[is.sig]

    if(length(test) > 0){
        rtracklayer::export.bed(test, paste0(outfile_prefix,"_significant.bed"))
    } else {
        warning("output empty! please lower the fdr threshold.")
    }
    ##merge regions with stats
    print(head(as.data.frame(merged$region)))
    print(head(tabcom))
    tabx<-as.data.frame(merged$region,stringsAsFactors=FALSE)
    tabx$name<-rownames(tabx)
    full_res<-as.data.frame(merge(x=tabx,y=tabcom,by.x="name",by.y="name"),stringsAsFactors=FALSE) 
    full_res<-full_res[,c(2:ncol(full_res),1)]
    print(sprintf("Colnames of result file are %s",colnames(full_res)))
    full_res[,2]<-full_res[,2]-1
    full_res[,2]<-format(full_res[,2], scientific = FALSE,trim=TRUE)
    full_res[,3]<-format(full_res[,3], scientific = FALSE,trim=TRUE)
    #write full results
    write.table(full_res,file="Full.results.bed",row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE)
    ##filter full result for FDR and LFC and write to output
    full_res.filt<-subset(full_res,(FDR<=fdrcutoff)&(abs(best.logFC)>=lfccutoff))
    if(nrow(full_res.filt)>0){
    write.table(full_res.filt,file="Filtered.results.bed",row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE)}else{system("touch Filtered.results.bed")}
    res.filt.up<-subset(full_res.filt,direction %in% "up")
    if(nrow(res.filt.up)>0){
    write.table(res.filt.up,file="Filtered.results.UP.bed",row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE)}else{system("touch Filtered.results.UP.bed")}
    res.filt.down<-subset(full_res.filt,direction %in% "down")
    if(nrow(res.filt.down)>0){
    write.table(res.filt.down,file="Filtered.results.DOWN.bed",row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE)}else{system("touch Filtered.results.DOWN.bed")}
    res.filt.mixed<-subset(full_res.filt,direction %in% "mixed")
    if(nrow(res.filt.mixed)>0){
    write.table(res.filt.mixed,file="Filtered.results.MIXED.bed",row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE)}else{system("touch Filtered.results.MIXED.bed")}
}

argv<-commandArgs(T)
sampleInfoFilePath<-argv[1]
mapped_dir<-argv[2]
peaks_dir<-argv[3]
fragmentSize.metric<-argv[4]
sampleSheet <- read.table(sampleInfoFilePath, header = TRUE, colClasses = c("character", "character"))
# name   condition
# SRR13171471    control
# SRR13171472    control
# SRR13171473    test
# SRR13171474    test
rownames(sampleSheet)<-sampleSheet$name
# param <- readParam(minq=20)
param <- readParam(max.frag = 500, pe = 'both')

cnames.sub<-unique(colnames(sampleSheet)[2:which(colnames(sampleSheet) %in% "condition")])
d<-as.formula(noquote(paste0("~",paste(cnames.sub,collapse="+"))))
sampleSheet$condition = factor(sampleSheet$condition )
sampleSheet$condition <- relevel(sampleSheet$condition, ref = as.character(sampleSheet$condition[1]))
designm <- model.matrix(d, data = sampleSheet)
rownames(designm)<-sampleSheet$name
designType <- "condition"

# bam.files<-c('/mnt/liuq/ATAC_seq/rawdata/PRJNA681776/mapped/SRR13171471.bam',
# '/mnt/liuq/ATAC_seq/rawdata/PRJNA681776/mapped/SRR13171472.bam',
# '/mnt/liuq/ATAC_seq/rawdata/PRJNA681776/mapped/SRR13171473.bam',
# '/mnt/liuq/ATAC_seq/rawdata/PRJNA681776/mapped/SRR13171474.bam')

# bam.files<-c('/mnt/liuq/ATAC_seq/rawdata/PRJNA681776/mapped/SRR13171471.bam',
# '/mnt/liuq/ATAC_seq/rawdata/PRJNA681776/mapped/SRR13171472.bam',
# '/mnt/liuq/ATAC_seq/rawdata/PRJNA681776/mapped/SRR13171475.bam',
# '/mnt/liuq/ATAC_seq/rawdata/PRJNA681776/mapped/SRR13171476.bam')

bam.files<-paste0(mapped_dir,sampleSheet$name,'.sorted.bam')

# https://github.com/maxplanck-ie/snakepipes/blob/master/snakePipes/shared/rscripts/DB_functions.R

# https://github.com/maxplanck-ie/snakepipes/blob/master/snakePipes/shared/rscripts/CSAW.R

pe = "both"
insert_size_metrics<-fragmentSize.metric
# insert_size_metrics<-'/mnt/liuq/ATAC_seq/rawdata/PRJNA681776/deepTools_qc/bamPEFragmentSize/HET1A_KYSE150/fragmentSize.metric.tsv'
d = read.delim(insert_size_metrics)
fraglength = median(d[,6])

data <- windowCounts(bam.files = bam.files, ext=fraglength, spacing=20, param=param, filter = 20)

# colnames(data) <- c('SRR13171471','SRR13171472','SRR13171473','SRR13171474')
# colnames(data) <- c('SRR13171471','SRR13171472','SRR13171475','SRR13171476')
colnames(data) <- sampleSheet$name
c
data<-data[,sampleSheet$name]
chip_object <- list(windowCounts = data, sampleSheet = sampleSheet,
                        design = designm, designType = designType, pe.param = param)

first_bam <- head(SummarizedExperiment::colData(chip_object$windowCounts)$bam.files, n = 1)
last_bam <- tail(SummarizedExperiment::colData(chip_object$windowCounts)$bam.files, n = 1)
# makeQCplots_chip_PE(bam.file = first_bam, outplot = "QCplots_first_sample.pdf", pe.param = param)
# makeQCplots_chip_PE(bam.file = last_bam, outplot = "QCplots_last_sample.pdf", pe.param = param)
fnames<-sampleSheet$name

# peaks <- c('/mnt/liuq/ATAC_seq/rawdata/PRJNA681776/peaks/genrich/SRR13171471/SRR13171471_peaks.narrowPeak','/mnt/liuq/ATAC_seq/rawdata/PRJNA681776/peaks/genrich/SRR13171472/SRR13171472_peaks.narrowPeak','/mnt/liuq/ATAC_seq/rawdata/PRJNA681776/peaks/genrich/SRR13171473/SRR13171473_peaks.narrowPeak','/mnt/liuq/ATAC_seq/rawdata/PRJNA681776/peaks/genrich/SRR13171474/SRR13171474_peaks.narrowPeak')
# peaks <- c('/mnt/liuq/ATAC_seq/rawdata/PRJNA681776/peaks/genrich/SRR13171471/SRR13171471_peaks.narrowPeak','/mnt/liuq/ATAC_seq/rawdata/PRJNA681776/peaks/genrich/SRR13171472/SRR13171472_peaks.narrowPeak','/mnt/liuq/ATAC_seq/rawdata/PRJNA681776/peaks/genrich/SRR13171475/SRR13171475_peaks.narrowPeak','/mnt/liuq/ATAC_seq/rawdata/PRJNA681776/peaks/genrich/SRR13171476/SRR13171476_peaks.narrowPeak')
peaks <- paste0(peaks_dir,sampleSheet$name,'/',sampleSheet$name,'_peaks.narrowPeak')

# allpeaks = lapply(snakemake@input[['peaks']], function(x) {
allpeaks = lapply(peaks, function(x) {
        peakfile<-x
        if(file.exists(peakfile) & file.info(peakfile)$size > 0){
            bed = read.delim(peakfile, header=FALSE)
            bed.gr = GRanges(seqnames = bed$V1, ranges = IRanges(start = bed$V2, end = bed$V3), name = bed$V4)
            }else{message(paste0("Skipping peakfile ",peakfile))
                    bed.gr=GRanges(c(seqnames=NULL,ranges=NULL,strand=NULL,name=NULL))}
            return(bed.gr)
    })
allpeaks <- Reduce(function(x,y) GenomicRanges::union(x,y), allpeaks)
keep <- overlapsAny(SummarizedExperiment::rowRanges(chip_object$windowCounts), allpeaks)
chip_object$windowCounts <- chip_object$windowCounts[keep,]
chip_object <- tmmNormalize_chip(chip_object, binsize = 10000, plotfile = "TMM_normalizedCounts.pdf")
chip_results <- getDBregions_chip(chip_object, plotfile = "DiffBinding_modelfit.pdf")
fdr<-0.05
writeOutput_chip(chip_results, outfile_prefix = "DiffBinding", fdrcutoff = fdr,lfccutoff=1)
## save data
message("Saving data")
save(chip_object, chip_results, file = "DiffBinding_analysis.Rdata")

#### SESSION INFO

sink("CSAW.session_info.txt")
sessionInfo()
sink()

# binned <- windowCounts(bam.files, bin=TRUE, width=10000, param=param)
# keep <- filterWindowsGlobal(data, binned)$filter > log2(5)
# data <- data[keep,]
# # 计算归一化因子。
# data <- normFactors(binned, se.out=data)

# y <- asDGEList(data)
# y <- estimateDisp(y, design)
# fit <- glmQLFit(y, design, robust=TRUE)
# results <- glmQLFTest(fit)
# # 多重检验校正。
# merged <- mergeResults(data, results$table, tol=1000L)
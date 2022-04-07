library(karyoploteR)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(data.table)
args = commandArgs(trailingOnly=TRUE)

# usage: Rscript karyoploter_methylation.R [methylation_file1] [methylation_file2] [dmr_sites] [ctcf_sites] [merged_sites] [sample1_name] [sample2_name] [extend_length] [output_folder]

print("Loading datasets...")
methylation_normal <- fread(args[1], select = c("chromosome", "start", "methylated_frequency"))
methylation_diagnosis <- fread(args[2], select = c("chromosome", "start", "methylated_frequency"))
dmr_sites <- read.table(args[3], sep="\t", strip.white=TRUE)
ctcf_sites <- read.table(args[4], sep="\t", strip.white=TRUE)
merged_sites <- read.table(args[5], sep="\t", strip.white=TRUE)
dir.create(as.character(args[9]))

for(i in 1:nrow(merged_sites)) {
    methylation_normal_subset <- subset(methylation_normal, chromosome==as.character(merged_sites[i,1]) & start >= as.integer(merged_sites[i,2]) - as.integer(args[8]) & start <= as.integer(merged_sites[i,3]) + as.integer(args[8]))
    methylation_diagnosis_subset <- subset(methylation_diagnosis, chromosome==as.character(merged_sites[i,1]) & start >= as.integer(merged_sites[i,2]) - as.integer(args[8]) & start <= as.integer(merged_sites[i,3]) + as.integer(args[8]))

    #put coverage data into vectors
    chr_meth_normal <- unlist(methylation_normal_subset[, "chromosome"], use.names = FALSE)
    x_meth_normal <- unlist(methylation_normal_subset[, "start"], use.names = FALSE)
    y_meth_normal <- unlist(methylation_normal_subset[, "methylated_frequency"], use.names = FALSE)

    chr_meth_diagnosis <- unlist(methylation_diagnosis_subset[, "chromosome"], use.names = FALSE)
    x_meth_diagnosis <- unlist(methylation_diagnosis_subset[, "start"], use.names = FALSE)
    y_meth_diagnosis <- unlist(methylation_diagnosis_subset[, "methylated_frequency"], use.names = FALSE)

    output_filename <- paste(as.character(args[9]), "/", as.character(args[7]), ".", as.character(merged_sites[i,1]), "_", as.integer(merged_sites[i,2]), "_", as.integer(merged_sites[i,3]), ".pdf", sep ='')
    pdf(output_filename, width=20, height=5)

    #plot the karyogram
    pp <- getDefaultPlotParams(plot.type=1)
    pp$leftmargin <- 0.12
    pp$topmargin <- 10
    pp$bottommargin <- 15
    pp$ideogramheight <- 10
    pp$data1inmargin <- 2

    zoom.region <- toGRanges(data.frame(as.character(merged_sites[i,1]), as.integer(merged_sites[i,2])-as.integer(args[8]), as.integer(merged_sites[i,3])+as.integer(args[8])))

    kp <- plotKaryotype(genome = "hg38", zoom=zoom.region, plot.params = pp)

    #add legend to plot
    par(xpd=TRUE)
    legend(0.4, 1.25, fill = c("blue", "red"), legend = c(args[6], args[7]), bty="n", ncol=2, cex=1.3)

    print("Plotting methylation...")

    kpAddLabels(kp, labels = "Methylation %", r0=0.23, r1=0.98, label.margin = 0.035, cex=1.2)
    kpAxis(kp, r0=0.23, r1=0.98, ymin=0, ymax=100, side=1, cex=0.8)
    kpAxis(kp, r0=0.23, r1=0.98, ymin=0, ymax=100, side=2, cex=0.8)
    kpPlotLoess(kp, chr=chr_meth_normal, x=x_meth_normal, y=y_meth_normal, r0=0.23, r1=0.98, conf.interval=NULL, span=0.05, ymin=0.0, ymax=1.0, col="blue")
    kpPlotLoess(kp, chr=chr_meth_diagnosis, x=x_meth_diagnosis, y=y_meth_diagnosis, r0=0.23, r1=0.98, conf.interval=NULL, span=0.05, ymin=0.0, ymax=1.0, col="red")

    kpPoints(kp, chr=chr_meth_normal, x=x_meth_normal, y=y_meth_normal, r0=0.23, r1=0.98, cex=0.4, col="blue")
    kpPoints(kp, chr=chr_meth_diagnosis, x=x_meth_diagnosis, y=y_meth_diagnosis, r0=0.23, r1=0.98, cex=0.4, col="red")

    kpPlotRegions(kp, data=dmr_sites, r0=0.23, r1=0.98, col=rgb(red = 1, green = 0, blue = 0, alpha = 0.2), border=rgb(red = 1, green = 0, blue = 0, alpha = 0.2))

    kpAddLabels(kp, labels = "CTCF Bindind Sites", r0=0.15, r1=0.18, label.margin = 0.035, cex=1.1)
    kpPlotRegions(kp, data=ctcf_sites, r0=0.15, r1=0.18, col="darkred")

    print("Plotting genes...")
    #plot gene locations
    genes.data <- makeGenesDataFromTxDb(TxDb.Hsapiens.UCSC.hg38.knownGene, karyoplot=kp, plot.transcripts = TRUE, plot.transcripts.structure = TRUE)
    genes.data <- addGeneNames(genes.data)
    genes.data <- mergeTranscripts(genes.data)
    kpAddLabels(kp, labels = "Genes", r0=0, r1=0.08, label.margin = 0.035, cex=1.2)
    kpPlotGenes(kp, data=genes.data, r0=0, r1=0.08)

    #plot base positions
    kpAddBaseNumbers(kp, tick.dist = 5000, minor.tick.dist = 2500, add.units = TRUE, cex=0.9, digits = 6)

    dev.off()
}

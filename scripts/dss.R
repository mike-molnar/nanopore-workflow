library(DSS)
require(bsseq)
require(tidyr)
require(dplyr)

tidy_df <- function(df, sampleName) {
  Cov <- rlang::sym(paste0("Cov.",sampleName))
  M <- rlang::sym(paste0("M.",sampleName))
  df %>%
    select(chr, start, end, 
           !!Cov:=called_sites, 
           !!M:=called_sites_methylated)
}

make_bs <- function(df) {
  M <- df %>% 
    select(starts_with("M.")) %>%
    rename_all(function(x) gsub("M.", "", x)) %>%
    mutate_all(function(x) ifelse(is.na(x), 0, x)) %>%
    as.matrix()
  Cov <- df %>% 
    select(starts_with("Cov.")) %>%
    rename_all(function(x) gsub("Cov.", "", x)) %>%
    mutate_all(function(x) ifelse(is.na(x), 0, x)) %>%
    as.matrix()
  BSseq(chr = df$chr, 
        pos = df$start,
        Cov = Cov,
        M = M,
        sampleNames = colnames(M))
}

args<-commandArgs(TRUE)

normal1_df <- read.table(args[1], header=TRUE, sep="\t", stringsAsFactors = FALSE) %>%
  mutate(start=start+1, end=end+2, sampleName="Normal") %>%
  dplyr::rename(chr=chromosome) %>%
  arrange(chr, start)

tumor1_df <- read.table(args[2], header=TRUE, sep="\t", stringsAsFactors = FALSE) %>%
  mutate(start=start+1, end=end+2, sampleName="Diagnosis") %>%
  dplyr::rename(chr=chromosome) %>%
  arrange(chr, start)

bs <- tidy_df(normal1_df, "Normal") %>%
  full_join(tidy_df(tumor1_df, "Diagnosis"),
            by=c("chr", "start", "end")) %>%
  group_by(chr, start) %>%
  summarise_if(is.numeric, sum) %>%
  ungroup() %>%
  make_bs()  

pData(bs)$phenoData <- c("Normal", "Diagnosis")
loci.idx <- which(rowSums(getCoverage(bs, type="Cov")==0) == 0)

mParam = MulticoreParam(workers=as.numeric(args[4]), progressbar=TRUE)

dml_test <- DMLtest(BSobj = bs[loci.idx,], 
           group1=c("Normal"), 
           group2=c("Diagnosis"), 
           BPPARAM=mParam, 
           smoothing=TRUE,
           smoothing.span=as.numeric(args[9]))

dmrs <- callDMR(dml_test, p.threshold=as.numeric(args[5]), dis.merge=as.numeric(args[6]), minCG=as.numeric(args[7]), minlen=as.numeric(args[8]))
write.csv(dmrs, file=args[3])

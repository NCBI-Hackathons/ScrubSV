#!/usr/local/bin/R
# check phased SNPs first and and if there is a phase swap flag
# then check unphased SNPs for homo/hetero swap

# install required packages
installRequiredPackages <- function(myRequired){
  myPacks <- rownames(installed.packages())
  toInstall <- myRequired[!myRequired %in% myPacks]
  install.packages(toInstall, repos = "https://cran.cnr.berkeley.edu/")
}

# strsplit2 borrowed from edgeR
strsplit2 <- function (x, split, ...) {
  x <- as.character(x)
  n <- length(x)
  s <- strsplit(x, split = split, ...)
  nc <- unlist(lapply(s, length))
  out <- matrix("", n, max(nc))
  for (i in 1:n) {
    if (nc[i]) 
      out[i, 1:nc[i]] <- s[[i]]
  }
  out
}

# filter for non-phased SNPs
filterNonPhased <- function(snpCharVec){
  (any(snpCharVec == "0/1") | any(snpCharVec == "1/0")) & 
    (any(snpCharVec == "1/1") | any(snpCharVec == "1|1"))
}

# filter for phased SNPs
filterPhased <- function(snpCharVec){
  any(snpCharVec == "0|1") & any(snpCharVec == "1|0") | 
    ((any(snpCharVec == "0|1") | any(snpCharVec == "1|0")) & 
       (any(snpCharVec == "1|1") | any(snpCharVec == "1/1")))
}

# function for reading in VCF containing SNPs in Deleted SVs
scrub_del <- function(delVCF){
  suppressMessages(library(magrittr))
  suppressMessages(library(plyr))
  # is input going to be gzipped?
  input <- 
    read.table(gzfile(delVCF), header = F, sep = "\t",
               colClasses = c(rep("NULL", 9), 
                              "character", "NULL", "integer", "NULL", 
                              rep("factor", 2), 
                              rep("NULL", 1))) %>% 
    `colnames<-` (., c("SNP", "POS", "ID", "FLAG")) %>% 
    mutate(SNP = strsplit2(.$SNP, ":")[,1], FLAG = NA)
  
  noPhaseFilterFLAG <-
    split(input, input$ID) %>% 
    lapply(., function(x){filterNonPhased(x[,1])}) %>% 
    unlist()
  
  phaseFilterFLAG <-
    split(input, input$ID) %>% 
    lapply(., function(x){filterPhased(x[,1])}) %>% 
    unlist()
  
  out <- 
    input[input$ID %in% c(names(which(noPhaseFilterFLAG)), names(which(phaseFilterFLAG))),] %>% 
    mutate(SNP = NULL, FLAG = "scrub_del") %>% 
    unique() %>% 
    `rownames<-` (., c())
  
  #saving a file could be a problem for parallel execution
  write.table(x = out, file = "scrub_del.tsv", quote = F, row.names = FALSE, sep = "\t")
}

installRequiredPackages(myRequired = c("magrittr", "plyr"))

args <- commandArgs(TRUE)

scrub_del(delVCF = as.character(args[1]))


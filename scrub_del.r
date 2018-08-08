#!/usr/local/bin/R
# check phased SNPs first and and if there is a phase swap flag
# then check unphased SNPs for homo/hetero swap

# install required packages
installRequiredPackages <- function(myRequired){
  myPacks <- rownames(installed.packages())
  toInstall <- myRequired[!myRequired %in% myPacks]
  if(length(toInstall != 0)){
    install.packages(toInstall, repos = "https://cran.cnr.berkeley.edu/")
  }
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
# minSNPs is the minimum number of erroneous SNPs to call a FLAG
filterErroneousSNPs <- function(snpCharVec, minSNPs){
  if(nrow(snpCharVec) == 1){
    # print("length is 1")
    return(FALSE)
  }
  tabSNPs <- table(snpCharVec)
  phased <- tabSNPs[grepl("[|]", names(tabSNPs))]
  notPhased <- tabSNPs[grepl("[/]", names(tabSNPs))]
  # if the Deletion doesn't contain enough SNPs don't flag it
  if(max(tabSNPs) <= minSNPs){
    # print("max SNPs is less than threshold")
    return(FALSE)
  }
  # if the deletion is multallelic flag it!
  if(any(grepl("2|3", names(tabSNPs)))){
    # print("this is multiallelic")
    return(TRUE)
  }
  # are these phased SNPs?
  if(any(phased) & length(phased) > 1) {
    # are there incongruous phased SNP calls?
    if(any(snpCharVec == "0|1") & any(snpCharVec == "1|0") | 
       ((any(snpCharVec == "0|1") | any(snpCharVec == "1|0")) & 
        (any(snpCharVec == "1|1") | any(snpCharVec == "1/1")))) {
      # are there enough to meet minSNPs requirement
      if(max(phased[-which.max(phased)]) > minSNPs){
        return(TRUE)
        # if there is no additional evidence then return false
      } else if(!any(notPhased)){
        return(FALSE)
      }
    } 
  }
  # if we made it this far and there are non-phased SNP calls
  if(any(notPhased) & length(notPhased) > 1){
    # if there are conflicting non-phased SNPs
    if((any(snpCharVec == "0/1") | any(snpCharVec == "1/0")) & 
       (any(snpCharVec == "1/1") | any(snpCharVec == "1|1"))) {
      # if the conflicting SNP call is greater than theshold call it TRUE
      if(max(notPhased[-which.max(notPhased)]) > minSNPs){
        return(TRUE)
        # else return false
      } else { return(FALSE) }
    }
  } 
  # If we made it this far, let's just say the SV deserves another shot
  # print("we made it to the end")
  return(FALSE)
}

# function for reading in VCF containing SNPs in Deleted SVs
# nSNPs = the number of SNPs required to determine flagged DEL SVs
scrub_del <- function(delVCF, minSNPs){
  suppressMessages(library(magrittr))
  suppressMessages(library(plyr))
  suppressMessages(library(data.table))
  # is input going to be gzipped?
  input <- 
    fread(delVCF, header = F, sep = "\t",
          select = c(10,12,14,15)) %>% 
    `colnames<-` (., c("SNP", "POS", "ID", "FLAG")) %>% 
    mutate(SNP = strsplit2(.$SNP, ":")[,1], FLAG = NA)
  
  filterFLAG <-
    split(input, input$ID) %>% 
    lapply(., function(x){filterErroneousSNPs(x[,1], minSNPs)}) %>% 
    unlist()
  
  out <- 
    input[input$ID %in% names(which(filterFLAG)),] %>% 
    mutate(SNP = NULL, FLAG = "scrub_del") %>% 
    unique() %>% 
    `rownames<-` (., c())
  
  #saving a file could be a problem for parallel execution
  # write.table(x = out, file = "scrub_del.tsv", quote = F, row.names = FALSE, sep = "\t")
  return(out)
}

installRequiredPackages(myRequired = c("magrittr", "plyr", "data.table"))

args <- commandArgs(TRUE)

scrub_del(delVCF = as.character(args[1]), minSNPs = as.numeric(args[2]))


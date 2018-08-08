#This script reads in a file from the output of intersecting inversions with overlapping SNVs, and then performs the following: 
#Identify duplications that overlap a user-defined minimum number of SNPs (num_snps_sd, at least 4), and flag duplications with highly fluctuating SNP coverages (read depth SDs for a given duplication larger than 2 standard deviations from the mean)
#Usage: Rscript scrub_inv.R vcffile num_snps_sd

#load required libraries, install them if not available
packages <- c("dplyr","tidyr","data.table")
package.check <- lapply(packages, FUN = function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dependencies = TRUE)
    library(x, character.only = TRUE)
  }
})

#store user arguments
args = commandArgs(trailingOnly=TRUE)
vcffile = args[1]
num_snps_sd = args[2]

#read file containing SNPs overlapping with SVs
dups <- fread(vcffile,sep="\t",header=F)
#remove SNPs with DP 0
dups <- dups %>% filter(!(V8==0))

#identifying SVs with SD of read depths greater than SD threshold (mean + 2SD)
sd_min_numsnp_threshold = ifelse(num_snps_sd > 4, num_snps_sd, 4)
dups_min_snps <- dups %>% 
  group_by(V6) %>% 
  summarize(nsnp=n()) %>% 
  filter(nsnp > sd_min_numsnp_threshold)
DP_SD <- dups %>% 
  filter(V6 %in% dups_min_snps$V6) %>% 
  group_by(V6) %>% 
  summarize(sdDP=sd(as.numeric(V8)))
#get SD cut-off value
DPSD_threshold <- mean(DP_SD$sdDP) + (2 * sd(DP_SD$sdDP)) 
DP_SD_thresh <- filter(DP_SD, sdDP>DPSD_threshold)
dups_DP_SD <- filter(dups, V6 %in% DP_SD_thresh$V6) %>% 
  mutate(FLAG = "scrubsv_highSD") %>% 
  select(V6,V5,FLAG) %>% 
  distinct()
colnames(dups_DP_SD) <- c("ID","POS","FLAG")

#write to table
write.table(dups_DP_SD, "flagged_inversions.txt",sep="\t",col.names = T,row.names = F,quote=F)
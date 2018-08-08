#This script reads in a file from the output of intersecting duplications with overlapping SNVs, and then performs the following two parts: 
#(i) Identify duplications that overlap a user-defined minimum number of SNPs (num_snps_sd, at least 4), and flag duplications with highly fluctuating SNP coverages (read depth SDs for a given duplication larger than 2 standard deviations from the mean)
#(ii) Flag duplications with a minimum number (user specified, num_multiallelic_snps) of multiallelic SNPs.
#Usage: Rscript scrub_dup.R vcffile num_snps_sd num_multiallelic_snps 

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
num_multiallelic_snps = args[3]

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
#write.table(dups_DP_SD, "dups_DP_highSD.txt",sep="\t",col.names = T,row.names = F,quote=F)

#group multiallelic SNVs
#to do: edit to allow for other multiallelic genotype combinations
multiallelic <- c("1/2","2/1","2|1","1|2","1/3","2/3","3/1","3/2")

#identify SVs overlapping multiallelic SNVs
dups_multiallelic <- dups %>% 
  separate(V3,sep=":",into=c("GT")) %>% 
  filter(GT %in% multiallelic) %>% 
  group_by(V6) %>% 
  summarize(num_multiallelic = n()) %>% 
  filter(num_multiallelic >= num_multiallelic_snps)
#flag only those SVs with multiallelic SNVs that have not been flagged for high DP SD before
dups_MA <- dups %>% 
  filter(V6 %in% dups_multiallelic$V6 && !(V6 %in% DP_SD_thresh$V6)) %>% 
  mutate(FLAG = "scrubdup_multiallelic") %>% 
  select(V6,V5,FLAG) %>% 
  distinct()
colnames(dups_MA) <- c("ID","POS","FLAG")

#write to table
#write.table(dups_MA, "dups_multiallelic_SNPs.txt", sep = "\t", col.names = T, row.names = F, quote = F)

flagged_duplications <- rbind(dups_DP_SD, dups_MA)
write.table(flagged_duplications, "flagged_duplications.txt", sep = "\t", col.names = T, row.names = F, quote = F)
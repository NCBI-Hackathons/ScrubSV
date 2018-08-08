##################################################################
##################################################################
###                                                            ###
### Rscript overlap_sv.R --arg1=input file --arg2=SVType       ###
### input file: output of vcf-query as a tab delimited table   ###
### with columns of ID, chromosome, POS, END, SVTYPE from vcf  ###
### file.                                                      ###
### SVType: one of DEL, DUP, INS, INV, TRA, BND                ###
### output files:                                              ###
###       1) *_flag.txt: a tab delmited table with additional  ###
###          flag for SVs overlapping with smaller ones.       ###
###       2) *_overlap.txt: a tab delimited table with info    ###
###          of ID, POS, and FLAG only.                        ###
###                                                            ###
##################################################################
##################################################################
 
args<-commandArgs(TRUE)

if(length(args)<1){
        args<-c("--help")
}

if("--help"%in%args){
        cat("The R script
                Arguments: --arg1=tab delimited text file from vcf-to-query output
                           --arg2=type of SV want to filter (DEL or DUP)
                Example: Rscript duplicated_sv.R --arg1=HG002_delly.tab --arg2=DEL\n")
        q(save="no")
}

first<-args[grep("arg1",args)]
filename<-sub("--arg1=","",first)
if(length(args)==2){
	second<-args[grep("arg2",args)]
	type=sub("--arg2=","",second)
}

library("GenomicRanges")

tab<-read.table(filename,header=F,sep="\t")
if(length(args)==2){
	tab<-tab[which(tab$V5==type),]
}

gchr<-paste("chr",tab[,2],sep="")
gstt<-ifelse(as.numeric(as.character(tab$V3))<as.numeric(as.character(tab$V4)),as.numeric(as.character(tab$V3)),as.numeric(as.character(tab$V4)))
gend<-ifelse(as.numeric(as.character(tab$V3))<as.numeric(as.character(tab$V4)),as.numeric(as.character(tab$V4)),as.numeric(as.character(tab$V3)))
glen<-gend-gstt
grns<-GRanges(gchr,IRanges(gstt,gend),symbol=tab$V1)

flag<-rep("",dim(tab)[[1]])

for(i in 1:(length(flag)-1)){
	sel_rns<-GRanges(gchr[i],IRanges(gstt[i],gend[i]))
	sel_len<-glen[i]
	new_ref<-grns[(i+1):length(flag),]
	hits<-findOverlaps(new_ref,sel_rns)
	if(length(hits)>0){
		oid<-as.character(mcols(new_ref)[queryHits(hits),1])
		overlaps<-pintersect(new_ref[queryHits(hits)],sel_rns[subjectHits(hits)])
		overlen<-glen[match(oid,tab$V1)]
		flag[i]<-ifelse(length(which(overlen<sel_len))>0,"OVLP","")
	}
}

new_tab<-cbind(tab,flag)
new_out<-paste(strsplit(filename,"[.]")[[1]][1],"_",type,"_flag.txt",sep="")
colnames(new_tab)<-c("ID","Chromosome","POS","End","SVTYPE","FLAG")
write.table(new_tab,new_out,row.names=F,col.names=T,quote=F,sep="\t")
new_tab<-new_tab[,c(1,3,6)]
new_tab<-new_tab[which(new_tab[,3]=="OVLP"),]
new_out<-gsub("flag","overlap",new_out)
write.table(new_tab,new_out,row.names=F,col.names=T,quote=F,sep="\t")


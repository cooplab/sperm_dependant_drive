#https://github.com/fmichonneau/rotl

#library(devtools)
#install_github("fmichonneau/rncl")
#install_github("fmichonneau/rotl")


library("rotl")
library("rncl")
library("ape")



write.csv(file="~/Dropbox/Ideas/Om/scripts/redone_meiotic_table.csv",cbind(meiotic.table[reordered.meiotic.table,],resolved_names.no.NAs))

new.meiotic.table<-read.csv(file="~/Dropbox/Ideas/Om/scripts/redone_meiotic_table.csv",as.is=TRUE)[,-1]

new.meiotic.table[grep("Lymnea",new.meiotic.table$unique_name),colnames(resolved_names.no.NAs)]<-tnrs_match_names("Lymnaeidae")
new.meiotic.table[grep("Pectinaria",new.meiotic.table$unique_name),colnames(resolved_names.no.NAs)]<-tnrs_match_names("Pectinariidae")

tr <- tol_induced_subtree(ott_ids=new.meiotic.table$ott_id)

old.names<- tr$tip.label
 new.names<-sapply(tr$tip.label,function(x){strsplit(x,"_")[[1]][1]})
 phylogeny.ott<-sapply(tr$tip.label,function(x){strsplit(x,"ott")[[1]][2]})
 
order.meiotic.table.phy<- match(phylogeny.ott,new.meiotic.table$ott_id)
 tr$tip.label<-paste(new.meiotic.table[order.meiotic.table.phy,1],new.names)


sperm.entry.stage<-new.meiotic.table[order.meiotic.table.phy,]$V3
sperm.entry.stage[sperm.entry.stage ==" GV (-MI)"]<-" GV-MI"
sperm.entry.stage[sperm.entry.stage ==" Pron. "]<- " Pron."
 sperm.entry.stage[sperm.entry.stage =="GV"]<- " GV"
stage.cols<-rainbow(length(unique( sperm.entry.stage)))


unfertilized.stage<-new.meiotic.table[order.meiotic.table.phy,]$V4
unfertilized.stage[unfertilized.stage =="Pron. "]<-" Pron."
unfertilized.stage[unfertilized.stage ==" ?"]<-"?"

stage.cols<-c(rainbow(length(unique( sperm.entry.stage))),"grey")
names(stage.cols)<- c(unique(sperm.entry.stage),"?")
fertilization.site<-new.meiotic.table[order.meiotic.table.phy,]$V5
  fertilization.site<-new.meiotic.table[order.meiotic.table.phy,]$V5
  fertilization.site[grep("Internal",fertilization.site)]<-"Internal"
  fertilization.site[grep("External",fertilization.site)]<-"External"
    fertilization.pch<-c(24,25,23,NA)
names(fertilization.pch)<- c( "External"   ,  "Internal" ,    "Intraovarian", " ? " )
# Abbreviations: GV germinal vesicle stage; Pron. pronucleus; MI metaphase I; Al anaphase I; MII metaphase II; AII anaphase II. 

pdf(file=paste(directory,"../meiotic_phylogeny.pdf",sep=""))

 par(mar=c(0,0,0,0))
 plot.phylo(tr,cex=0.75,label.offset=7) #,type="fan")

legend("topleft",pch=c(rep(19,length(stage.cols)),fertilization.pch) ,legend=c(names(stage.cols),names(fertilization.pch)),col=c(stage.cols,rep("black",length(fertilization.pch))))
 points(rep(78, length(tr$tip.label)), 1:length(tr$tip.label), pch=21, bg=stage.cols[sperm.entry.stage], cex=1) 
 points(rep(80, length(tr$tip.label)), 1:length(tr$tip.label), pch=22, bg=stage.cols[unfertilized.stage], cex=1) 
 points(rep(82, length(tr$tip.label)), 1:length(tr$tip.label), pch= fertilization.pch, col="black",bg="black", cex=0.7) 

dev.off()



#	reading in OCR text

if(FALSE){
meiotic.table<-read.csv("~/Dropbox/Ideas/Om/writeup_and_notes/evolution/meiotic_arrest_in_animal_oocytes_Masui_1_table.csv",head=FALSE,as.is=TRUE,skip=3)

genus<-unlist(sapply(meiotic.table[,2],strsplit," "))

#genus<-read.csv("~/Dropbox/Ideas/Om/writeup_and_notes/evolution/genus.txt",head=FALSE,as.is=TRUE)

genus<-unlist(c(genus))
genus<-genus[genus!=" "]
genus<-genus[genus!=""]
genus[genus=="bat"]<-"Chiroptera"

resolved_names <- tnrs_match_names(genus)

resolved_names.no.NAs<-resolved_names[!is.na(resolved_names$ott_id),]

reordered.meiotic.table<-sapply(resolved_names.no.NAs$search_string,function(my.name){
	this.one<-grep(my.name,meiotic.table[,2],ignore.case=TRUE)
	if(length(this.one)>1){ if(my.name=="sycon"){this.one<-1};if(my.name=="mus"){this.one<-35}  }
	if(length(this.one)>1)   recover()
	return(this.one)
})

write.csv(file="~/Dropbox/Ideas/Om/scripts/redone_meiotic_table.csv",cbind(meiotic.table[reordered.meiotic.table,],resolved_names.no.NAs))

}
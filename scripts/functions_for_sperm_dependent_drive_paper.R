

library("RColorBrewer")   #color scheme used in many figures

###BASIC FUNCTIONS TO SETUP AND ITERATE DRIVE MODEL

iterate.1.locus.drive<-function(s.array,num.iterations,female.transmission.probs,my.geno.freqs,initialize.allele.freqs=c(0.999,0.001,0),selfing.rate=0){
	geno.freqs<-rep(NA,6)
	names(geno.freqs)<-geno.names

	my.freqs<-numeric()
	if(missing(my.geno.freqs)){
       #     print("running 2 allele system to stablity")
		allele.freqs<-initialize.allele.freqs	
		names(allele.freqs)<-c("1","2","3")	

	
		geno.freqs["11"]<-allele.freqs["1"]^2
		geno.freqs["22"]<-allele.freqs["2"]^2
		geno.freqs["33"]<-allele.freqs["3"]^2

		geno.freqs["12"]<-2*allele.freqs["1"]*allele.freqs["2"]
		geno.freqs["13"]<-2*allele.freqs["1"]*allele.freqs["3"]
		geno.freqs["23"]<-2*allele.freqs["3"]*allele.freqs["2"]
		stopifnot(sum(geno.freqs)>0.9999999999 & sum(geno.freqs)< 1.0000000001)

	}else{
		print("initiating from my.geno.freqs")
		geno.freqs[c("11","12","22")]<-my.geno.freqs[["geno.freqs"]][nrow(my.geno.freqs[["geno.freqs"]]),c("11","12","22")]*.99
		geno.freqs[c("13","23")]<-my.geno.freqs[["allele.freqs"]][nrow(my.geno.freqs[["geno.freqs"]]),1:2]*.99*.01
		geno.freqs["33"]<-.01^2
		
		}

	my.freqs<-rbind(my.freqs,geno.freqs)

	for(i in 1:num.iterations){
	#selection
		geno.freqs["11"]<-geno.freqs["11"] * (1-s.array["11"])
		geno.freqs["12"]<-geno.freqs["12"] * (1-s.array["12"])
		geno.freqs["13"]<-geno.freqs["13"] * (1-s.array["13"])

		geno.freqs["22"]<-geno.freqs["22"] * (1-s.array["22"])
		geno.freqs["23"]<-geno.freqs["23"] * (1-s.array["23"])

		geno.freqs["33"]<-geno.freqs["33"] * (1-s.array["33"])
	
		geno.freqs<-geno.freqs/sum(geno.freqs)	

		geno.freqs.array<-outer(rep(geno.freqs,each=2),rep(geno.freqs,each=2),"*")  ##one for each allele

	geno.freqs.array<-(1-selfing.rate) * geno.freqs.array  ##selfing  (1-s) * pi *pj
	diag(geno.freqs.array)<-diag(geno.freqs.array) +  selfing.rate * rep(geno.freqs,each=2)  #selfing  (1-s) * pi *pj + delta_{ij} s pi

	##transmission
		tot.trans<-geno.freqs.array*female.transmission.probs*0.5

		new.geno.freqs<-rep(NA,6)	
		names(new.geno.freqs)<-geno.names
		new.geno.freqs["11"]<-sum(tot.trans[transmitted.allele==1,transmitted.allele==1])
		new.geno.freqs["22"]<-sum(tot.trans[transmitted.allele==2,transmitted.allele==2])
		new.geno.freqs["33"]<-sum(tot.trans[transmitted.allele==3,transmitted.allele==3])

		new.geno.freqs["12"]<-sum(tot.trans[transmitted.allele==1,transmitted.allele==2])+sum(tot.trans[transmitted.allele==2,transmitted.allele==1])
		new.geno.freqs["13"]<-sum(tot.trans[transmitted.allele==1,transmitted.allele==3])+sum(tot.trans[transmitted.allele==3,transmitted.allele==1])
		new.geno.freqs["23"]<-sum(tot.trans[transmitted.allele==2,transmitted.allele==3])+sum(tot.trans[transmitted.allele==3,transmitted.allele==2])

		geno.freqs<-new.geno.freqs
		my.freqs<-rbind(my.freqs,geno.freqs)

	}

	freq.1<-my.freqs[,"11"]+0.5*rowSums(my.freqs[,c("12","13")])
	freq.2<-my.freqs[,"22"]+0.5*rowSums(my.freqs[,c("12","23")])
	freq.3<-my.freqs[,"33"]+0.5*rowSums(my.freqs[,c("13","23")])


	geno.and.allele.freqs<-list()
	geno.and.allele.freqs[["geno.freqs"]]<-my.freqs
	geno.and.allele.freqs[["allele.freqs"]]<-cbind(freq.1,freq.2,freq.3)

	return(geno.and.allele.freqs) 
 }
 
 
  make.sperm.dep.female.transmission.prob<-function(D){
 	geno.names<-c("11", "12","13", "22", "23", "33")
	transmission.probs<-matrix(0.5,nrow=12,ncol=12)
	rownames(transmission.probs)<-rep(geno.names,each=2)
	colnames(transmission.probs)<-rep(geno.names,each=2)
	transmitted.allele<-c(1,1,1,2,1,3,2,2,2,3,3,3)
 	
	##sperm alleles 1 and 2 allow: 2 and 3 female hets 12 to drive at D1
	transmission.probs[transmitted.allele %in% 1:2,colnames(transmission.probs) =="12" & transmitted.allele == 1] <- 1-D["1 or 2","12"]
	transmission.probs[transmitted.allele %in% 1:2,colnames(transmission.probs) =="12" & transmitted.allele == 2] <- D["1 or 2","12"]
	
	##sperm alleles 1 and 2 allow: 2 and 3 female hets 12 to drive at D13
	transmission.probs[transmitted.allele %in% 1:2,colnames(transmission.probs) =="13" & transmitted.allele == 1] <- 1-D["1 or 2","13"]
	transmission.probs[transmitted.allele %in% 1:2,colnames(transmission.probs) =="13" & transmitted.allele == 3] <- D["1 or 2","13"]
	
	##sperm allele 3 allow: 2 and 3 hets 12 and 13 to drive at D2
	transmission.probs[transmitted.allele %in% 3,colnames(transmission.probs) == "13" & transmitted.allele == 1] <- 1-D["3","13"]
	transmission.probs[transmitted.allele %in% 3,colnames(transmission.probs) =="13" & transmitted.allele == 3] <- D["3","13"]
	
	transmission.probs[transmitted.allele %in% 3,colnames(transmission.probs) == "12" & transmitted.allele == 1] <- 1-D["3","12"]
	transmission.probs[transmitted.allele %in% 3,colnames(transmission.probs) =="12" & transmitted.allele == 2] <- D["3","12"]
	
	transmission.probs[transmitted.allele %in% 1:2,colnames(transmission.probs) == "23" & transmitted.allele == 2] <- 1-D["1 or 2","23"]
	transmission.probs[transmitted.allele %in% 1:2,colnames(transmission.probs) =="23" & transmitted.allele == 3] <- D["1 or 2","23"]
	
	transmission.probs[transmitted.allele %in% 3,colnames(transmission.probs) == "23" & transmitted.allele == 2] <- 1-D["3","23"]
	transmission.probs[transmitted.allele %in% 3,colnames(transmission.probs) =="23" & transmitted.allele == 3] <- D["3","23"]
	
	
	
transmission.probs
}
 
make.male.geno.dep.female.transmission.prob<-function(D){
 	geno.names<-c("11", "12","13", "22", "23", "33")
	transmission.probs<-matrix(0.5,nrow=12,ncol=12)
	rownames(transmission.probs)<-rep(geno.names,each=2)
	colnames(transmission.probs)<-rep(geno.names,each=2)
	transmitted.allele<-c(1,1,1,2,1,3,2,2,2,3,3,3)
 	
 	
	transmission.probs[,colnames(transmission.probs) =="12" & transmitted.allele == 1] <- 1-D["1 or 2","12"]
	transmission.probs[,colnames(transmission.probs) =="12" & transmitted.allele == 2] <- D["1 or 2","12"]

	transmission.probs[,colnames(transmission.probs) =="13" & transmitted.allele == 1] <- 1-D["1 or 2","13"]
	transmission.probs[,colnames(transmission.probs) =="13" & transmitted.allele == 3] <- D["1 or 2","13"]

	transmission.probs[,colnames(transmission.probs) == "23" & transmitted.allele == 2] <- 1-D["1 or 2","23"]
	transmission.probs[,colnames(transmission.probs) =="23" & transmitted.allele == 3] <- D["1 or 2","23"]


	##males carrying allele 3 alter rate of female drive in dominant fashion. 
	
	male.geno.3<-grep("3",rownames(transmission.probs))
	
	transmission.probs[male.geno.3,colnames(transmission.probs) == "13" & transmitted.allele == 1] <- 1-D["3","13"]
	transmission.probs[male.geno.3,colnames(transmission.probs) =="13" & transmitted.allele == 3] <- D["3","13"]

	transmission.probs[male.geno.3,colnames(transmission.probs) == "12" & transmitted.allele == 1] <- 1-D["3","12"]
	transmission.probs[male.geno.3,colnames(transmission.probs) =="12" & transmitted.allele == 2] <- D["3","12"]

	transmission.probs[male.geno.3,colnames(transmission.probs) == "23" & transmitted.allele == 2] <- 1-D["3","23"]
	transmission.probs[male.geno.3,colnames(transmission.probs) =="23" & transmitted.allele == 3] <- D["3","23"]
	recover()
transmission.probs
}


make.male.geno.dep.female.transmission.prob.w.dom<-function(D){
 	geno.names<-c("11", "12","13", "22", "23", "33")
	transmission.probs<-matrix(0.5,nrow=12,ncol=12)
	rownames(transmission.probs)<-rep(geno.names,each=2)
	colnames(transmission.probs)<-rep(geno.names,each=2)
	transmitted.allele<-c(1,1,1,2,1,3,2,2,2,3,3,3)
 	
 	
	transmission.probs[,colnames(transmission.probs) =="12" & transmitted.allele == 1] <- 1-D["1 or 2","12"]
	transmission.probs[,colnames(transmission.probs) =="12" & transmitted.allele == 2] <- D["1 or 2","12"]

	transmission.probs[,colnames(transmission.probs) =="13" & transmitted.allele == 1] <- 1-D["1 or 2","13"]
	transmission.probs[,colnames(transmission.probs) =="13" & transmitted.allele == 3] <- D["1 or 2","13"]

	transmission.probs[,colnames(transmission.probs) == "23" & transmitted.allele == 2] <- 1-D["1 or 2","23"]
	transmission.probs[,colnames(transmission.probs) =="23" & transmitted.allele == 3] <- D["1 or 2","23"]


	##males carrying allele 3 alter rate of female drive in dominant fashion. 
	
	male.heterozy<-rownames(transmission.probs) %in%  c("13","23")
	male.homozy<-rownames(transmission.probs) =="33"
	
	
	transmission.probs[male.heterozy,colnames(transmission.probs) == "13" & transmitted.allele == 1] <- 1-D["3.het","13"]
	transmission.probs[male.heterozy,colnames(transmission.probs) =="13" & transmitted.allele == 3] <- D["3.het","13"]

	transmission.probs[male.heterozy,colnames(transmission.probs) == "12" & transmitted.allele == 1] <- 1-D["3.het","12"]
	transmission.probs[male.heterozy,colnames(transmission.probs) =="12" & transmitted.allele == 2] <- D["3.het","12"]

	transmission.probs[male.heterozy,colnames(transmission.probs) == "23" & transmitted.allele == 2] <- 1-D["3.het","23"]
	transmission.probs[male.heterozy,colnames(transmission.probs) =="23" & transmitted.allele == 3] <- D["3.het","23"]

##homozy. male effect
	transmission.probs[male.homozy,colnames(transmission.probs) == "13" & transmitted.allele == 1] <- 1-D["3","13"]
	transmission.probs[male.homozy,colnames(transmission.probs) =="13" & transmitted.allele == 3] <- D["3","13"]

	transmission.probs[male.homozy,colnames(transmission.probs) == "12" & transmitted.allele == 1] <- 1-D["3","12"]
	transmission.probs[male.homozy,colnames(transmission.probs) =="12" & transmitted.allele == 2] <- D["3","12"]

	transmission.probs[male.homozy,colnames(transmission.probs) == "23" & transmitted.allele == 2] <- 1-D["3","23"]
	transmission.probs[male.homozy,colnames(transmission.probs) =="23" & transmitted.allele == 3] <- D["3","23"]

#	recover()
transmission.probs
}


run.iterations<-function(D,s.array,sperm.or.geno.dependent,num.iter=4000,plot.initial.rise=FALSE){
	if(sperm.or.geno.dependent=="sperm"){ female.transmission.probs<-make.sperm.dep.female.transmission.prob(D)}
	if(sperm.or.geno.dependent=="geno"){ female.transmission.probs<-make.male.geno.dep.female.transmission.prob(D)}
	
	old.geno.freqs<-iterate.1.locus.drive(s.array,num.iterations=num.iter,female.transmission.probs=female.transmission.probs)
	new.geno.freqs<-iterate.1.locus.drive(s.array=s.array,num.iterations=num.iter,female.transmission.probs=female.transmission.probs,my.geno.freqs=old.geno.freqs)
	
	if(plot.initial.rise) new.geno.freqs[["allele.freqs"]]<-rbind(old.geno.freqs[["allele.freqs"]],new.geno.freqs[["allele.freqs"]])
	plot.freqs(new.geno.freqs)
	
}


directory<-"~/Dropbox/Ideas/Om/scripts/"
geno.names<-c("11", "12","13", "22", "23", "33")
transmitted.allele<-c(1,1,1,2,1,3,2,2,2,3,3,3)


####FUNCTIONS TO TEST WHETHER DRIVER SPREADS
#### test.invasion.self.prop used to construct invasion diagrams, allele intro'd at freq. x
#### checks whether it increases in frequency
test.invasion.self.prop<-function(d,s.het,s.hom,x=.001,steps=10,sperm.or.male="sperm",selfing.rate=0){
	s.array<-rep(0,6)
	names(s.array)<-c("11","12","13","23","22","33")
	s.array["13"]<- s.het
	s.array["33"]<- s.hom
	
	D<-matrix(1/2,nrow=2,ncol=3,dimnames=list(c("1 or 2","3"),c("12","13","23"))) ##entries are drive coeff of the 2nd 
	D["1 or 2","12"]<-0.5
	D["1 or 2","13"]<-0.5
	D["3","12"]<-.5
	D["3","13"]<-d
	D["1 or 2","23"]<- 1-0.5
	D["3","23"]<-0.5
	stopifnot(sperm.or.male=="sperm" | sperm.or.male=="male")
	if(sperm.or.male=="sperm"){ female.transmission.probs<-make.sperm.dep.female.transmission.prob(D)}
	if(sperm.or.male=="male"){ female.transmission.probs<-make.male.geno.dep.female.transmission.prob(D) }
	old.geno.freqs<-iterate.1.locus.drive(s.array,num.iterations=steps,female.transmission.probs=female.transmission.probs,initialize.allele.freqs =c(1-x, 0,x),selfing.rate=selfing.rate)
#recover()
	invading<-old.geno.freqs[["allele.freqs"]][steps,3]>old.geno.freqs[["allele.freqs"]][1,3]
	return(invading)
}

test.fixation.of.simple.driver<-function(d,s.het,s.hom,x=.999,steps=10,selfing.rate=0){
	s.array<-rep(0,6)
	names(s.array)<-c("11","12","13","23","22","33")
	s.array["12"]<- s.het
	s.array["22"]<- s.hom
	
	D<-matrix(1/2,nrow=2,ncol=3,dimnames=list(c("1 or 2","3"),c("12","13","23"))) ##entries are drive coeff of the 2nd 
	D["1 or 2","12"]<- d
	D["1 or 2","13"]<-0.5
	D["3","12"]<- 0.5
	D["3","13"]<- 0.5
	D["1 or 2","23"]<- 1-0.5
	D["3","23"]<-0.5
	female.transmission.probs<-make.sperm.dep.female.transmission.prob(D)
	old.geno.freqs<-iterate.1.locus.drive(s.array=s.array,num.iterations=steps,female.transmission.probs=female.transmission.probs,initialize.allele.freqs =c(1-x, x,0),selfing.rate=selfing.rate)
	fixing<-old.geno.freqs[["allele.freqs"]][steps,2]>old.geno.freqs[["allele.freqs"]][1,2]
	return(fixing)
}

####FUNCTIONS TO MAKE FIGS FOR PAPER

##############Phase diagram figure for Figure 1 with inset

make.invasion.grid.fig.1<-function(selfing.rate=0){
	D<-matrix(1/2,nrow=2,ncol=3,dimnames=list(c("1 or 2","3"),c("12","13","23"))) ##entries are drive coeff of the 2nd allele listed against the 1st.
	
	d.range<-seq(0.5,1,length=200)
	s3.range<-seq(0,1,length=200)
	
	invasion.grid<-sapply(s3.range,function(s3){
		sapply(d.range,function(d){
			test.invasion.self.prop(d=d,s.het=0,s.hom=s3,selfing.rate=selfing.rate)
		})
	})
	
	invasion.lines.s3.cutoff<-s3.range[apply(invasion.grid,1,function(x){max(which(x))})]
	
	fixation.grid<-sapply(s3.range,function(s3){
		sapply(d.range,function(d){
			test.fixation.of.simple.driver(d=d,s.het=0,s.hom=s3,selfing.rate=selfing.rate)
		})
	})
	fixation.lines.s3<-s3.range[apply(fixation.grid,1,function(x){max(which(x))})]
save(file=paste(directory,"invasion_grid_homozy_cost",format(selfing.rate,dig=3),".Robj",sep=""),fixation.lines.s3,fixation.grid,invasion.lines.s3.cutoff,invasion.grid,selfing.rate)
}




make.Figure.1<-function(grid.file="invasion_grid_homozy_cost.Robj"){
	show(load(paste(directory,grid.file,sep="")))

my.cols<-c("red","blue")  #brewer.pal(2,name="Dark2")
plot(c(0.5,1),c(0,1),xlab="",ylab="",type="n")
mtext("Drive coefficient, d",side=1,line=2)
mtext("selection against homozygotes, s",side=2,line=2)

polygon(x=c(d.range,rev(d.range)),y=c(fixation.lines.s3,rep(1,length(d.range))),col="white")
polygon(x=c(d.range,rev(d.range)),y=c(fixation.lines.s3,rev(invasion.lines.s3.cutoff)),col=my.cols[1],border=NA)
polygon(x=c(d.range,rev(d.range)),y=c(invasion.lines.s3.cutoff,rep(0,length(d.range))),col=my.cols[2],border=NA)
text(0.8,0.35,"Simple driver invades but can't fix",srt=25)
text(0.85,0.27,"Simple driver can invade & fix",col="white",srt=15)
text(0.8,0.05,"Simple driver & Self promoter driver can invade & fix",col="white")
lines(d.range,(2*d.range- 1)/2,col="black",lwd=3,lty=1) ##simple driver fixes below this line
lines(d.range,(2*d.range- 1)/(4*d.range),col="black",lwd=3,lty=2) ##self prom. fixes below this line

text(x=0.95,y=0.95,"B.",cex=2)

par(fig=c(grconvertX(0.51, from = "user", to = "ndc"), grconvertX(0.75, from = "user", to = "ndc"), grconvertY(0.4, from = "user", to = "ndc"), grconvertY(1, from = "user", to = "ndc")), new = TRUE); 
par(mar=c(2,2,1,0))

drive.coeffs<-c(0.6,.7,.8)
my.point.cols<-c("red","blue","orange")

	s.array<-rep(0,6)
	names(s.array)<-c("11","12","13","23","22","33")
	s.array["33"]<- 0.1


plot(c(0,2000),c(0,1),type="n",xlab="",ylab="",axes=FALSE)
box();
axis(1,tick=FALSE,line=-0.5);axis(1,label=FALSE);
axis(2,tick=FALSE,line=-0.5);axis(2,label=FALSE);
mtext("generations",side=1,line=1.5)
mtext("frequency",side=2,line=1.5)

for(i in 1:3){
		 D["1 or 2","12"]<-0.5
	 D["1 or 2","13"]<-0.5
	 D["3","12"]<-.5
	 D["3","13"]<-drive.coeffs[i]
	 D["1 or 2","23"]<- 1-0.5
	 D["3","23"]<-0.5
		 female.transmission.probs<-make.sperm.dep.female.transmission.prob(D)

	old.geno.freqs<-iterate.1.locus.drive(s.array=s.array,num.iterations=2000,female.transmission.probs=female.transmission.probs,initialize.allele.freqs =c(0.99, 0,0.01 ))
	lines(old.geno.freqs[["allele.freqs"]][,3],col=my.point.cols[i],lwd=2)
	}
text(x=50,y=0.95,"C.",cex=1)
#dev.copy2eps(file=paste(directory,"invasion_space_recessive_driver.eps",sep=""))
}


plot.freqs<-function(new.geno.freqs,plot.HWE.dev=FALSE){
	plot(new.geno.freqs[["allele.freqs"]][,1],ylim=c(0,1),type="l",ylab="frequency",xlab="generations",lwd=2)
 	lines(new.geno.freqs[["allele.freqs"]][,2],col="red",lwd=2)
 	lines(new.geno.freqs[["allele.freqs"]][,3],col="blue",lwd=2)

	if(plot.HWE.dev){
#####deviations away from HWe 
 		plot(new.geno.freqs[["geno.freqs"]][,"11"] / new.geno.freqs[["allele.freqs"]][,1]^2,ylim=c(.85,1.15),type="l",col="red",ylab="frequency",xlab="generations")
		lines(new.geno.freqs[["geno.freqs"]][,"22"] / new.geno.freqs[["allele.freqs"]][,2]^2,col="blue")
		lines(new.geno.freqs[["geno.freqs"]][,"33"] /new.geno.freqs[["allele.freqs"]][,3]^2,col="green")

		lines(new.geno.freqs[["geno.freqs"]][,"12"] / (2*new.geno.freqs[["allele.freqs"]][,1]*new.geno.freqs[["allele.freqs"]][,2]),col="red",lty=2)
		lines(new.geno.freqs[["geno.freqs"]][,"13"] / (2*new.geno.freqs[["allele.freqs"]][,1]*new.geno.freqs[["allele.freqs"]][,3]),col="blue",lty=2)
		lines(new.geno.freqs[["geno.freqs"]][,"23"] / (2*new.geno.freqs[["allele.freqs"]][,3]*new.geno.freqs[["allele.freqs"]][,2]),col="green",lty=2)
	}

}

figures.of.traj.for.figure.2<-function(){
	layout(t(1:2)); par(mar=c(3,3,1,1))


	##driver that avoids itself 
	D<-matrix(1/2,nrow=2,ncol=3,dimnames=list(c("1 or 2","3"),c("12","13","23"))) ##entries are drive coeff of the 2nd allele listed against the 1st.
	D1 = 1
	D2 = 0.5
	s.array<-rep(0,6);names(s.array)<-c("11","12","13","23","22","33")
	s.array["23"]<- 1; s.array["33"]<- 1; s.array["22"]<- 1
	
	 D["1 or 2","12"]<-D1
	 D["1 or 2","13"]<-D1
	 D["3","12"]<-D2
	 D["3","13"]<-D2
	 
	run.iterations(D=D,s.array=s.array,sperm.or.geno.dependent="sperm",num.iter=200,plot.initial.rise=TRUE)
	text(x=125,y=0.9,"Simple Driver\n invades")
	
	text(100,0.3,"A")
	text(100,0.7,"B")	
	text(270,0.65,"B-")	

	text(x=310,y=0.8,"Sperm-acting supressor\n in coupling phase \n invades")
	mtext("generations",side=1,line=2)
	mtext("Frequency",side=2,line=2)
	
		##female driver followed by sperm supressor
	
	s.array<-rep(0,6);names(s.array)<-c("11","12","13","23","22","33")
	s.array["22"]<- 1

	D1 = 1
	D2 = 0.5
	D<-matrix(1/2,nrow=2,ncol=3,dimnames=list(c("1 or 2","3"),c("12","13","23")))
	 D["1 or 2","12"]<-D1
	 D["1 or 2","13"]<-D1
	 D["3","12"]<-D2
	 D["3","13"]<-D2
	 D["1 or 2","23"]<- 1-D1
	 D["3","23"]<-0.5
	
	run.iterations(D=D,s.array=s.array,sperm.or.geno.dependent="sperm",num.iter=200,plot.initial.rise=TRUE)
	text(x=125,y=0.9,"Simple Driver\n invades")
	text(100,0.3,"A")
	text(100,0.7,"B")	
	text(250,0.9,"A-")	
	
	text(x=325,y=0.5,"Sperm-acting \n drive supressor\n in repulsion phase \n invades")
	mtext("generations",side=1,line=2)
	mtext("Frequency",side=2,line=2)
	dev.copy2eps(file=paste(directory,"trajectories_of_sperm_based_supressors.eps",sep=""))
}




#####  Supplementary figure of bistability of self promoting driver with recessive fitness cost.

make.grid.for.bistability.of.self.promoting.driver.w.recess.cost<-function(){
	d.range<-seq(0.5,1,length=200)
	s3.range<-seq(0,1,length=200)

	x.range<-seq(0.00001,1,length=100);

	my.x.bistable<-numeric()

	for(i in 1:length(s3.range)){
		
		invasion.grid.bistable<-sapply(d.range,function(d){
			sapply(x.range,function(x){
				test.invasion.self.prop(d=d,s.het=0,s.hom=s3.range[i],x=x,step=10)
			})
		})
		cat(i," ")
		#image(d.range,x.range,invasion.grid.bistable,log="y")
		invasion.x.lines<-x.range[apply(invasion.grid.bistable,2,function(x){min(which(x))})]
	#	text(d.range[round(length(d.range)/2)],invasion.x.lines[round(length(d.range)/2)],s,col="red")
		my.x.bistable<-rbind(my.x.bistable,invasion.x.lines)
		save(file="homozyg_bistability.Robj",my.x.bistable)
	}
}

make.supp.fig.bistability.of.self.promoting.driver.w.recess.cost<-function(){
	d.range<-seq(0.5,1,length=200)
	s3.range<-seq(0,1,length=200)

	show(load("homozyg_bistability.Robj"))
	layout(t(1:2))
	par(mar=c(3,3,1,1))
	my.x.bistable<-t(my.x.bistable)
	my.x.bistable<-cbind(my.x.bistable,matrix(nrow=200,ncol=67)) ##fill out rest of array 
	image(d.range,s3.range,my.x.bistable,xlab="",ylab="")
	lines(d.range,(2*d.range- 1)/2,col="black",lwd=3,lty=1) ##simple driver fixes below this line
	lines(d.range,(2*d.range- 1)/(4*d.range),col="black",lwd=3,lty=2) ##self prom. fixes below this line
	mtext("Bistability of invasion w. self promoter w. homozyg. cost",side=3,line=0)
	mtext("Drive coefficient, d",side=1,line=2)
	mtext("Selection against homozygotes, s",side=2,line=2)
	contour(d.range,s3.range,my.x.bistable,add=TRUE,lty=2)
	text(0.8,0.6,"Self promoter cannot invade\n no matter what freq.")
	text(0.85,0.1,"Self promoter driver can\n invade & fix when rare",col="white")
	
	##HWE
	approx.bistable.freq<-outer(d.range,s3.range,function(d,s){(1-2*d+4*d*s)/(2*s*(2*d-1))})
	approx.bistable.freq<-apply(approx.bistable.freq,1,function(x){x[x>=1]<-NA;x[x<=0]<-0;x})
	image(d.range,s3.range,t(approx.bistable.freq),xlab="",ylab="")
	lines(d.range,(2*d.range- 1)/2,col="black",lwd=3,lty=1) ##simple driver fixes below this line
	lines(d.range,(2*d.range- 1)/(4*d.range),col="black",lwd=3,lty=2) ##self prom. fixes below this line
	mtext("Approximation to region of bistability",side=3,line=0)
	mtext("Drive coefficient, d",side=1,line=2)
	mtext("Selection against homozygotes, s",side=2,line=2)
	contour(d.range,s3.range,t(approx.bistable.freq),add=TRUE,lty=2)
	text(0.8,0.6,"Self promoter cannot invade\n no matter what freq.")
	text(0.85,0.1,"Self promoter driver can\n invade & fix when rare",col="white")
	
	dev.copy2eps(file="homozyg_bistability.eps")
}

#####  Supplementary figure of bistability of self promoting driver with additive fitness cost.
make.grid.for.bistability.of.self.promoting.driver.w.add.cost<- function(){
	my.x.bistable<-numeric()
	
	for(i in 1:length(my.s)){
		
		invasion.grid.bistable<-sapply(d.range,function(d){
			sapply(x.range,function(x){
				test.invasion.self.prop(d=d,s.het=my.s[i],s.hom=my.s[i]*2.0,x=x,step=10)
			})
		})
		#image(d.range,x.range,invasion.grid.bistable,log="y")
		invasion.x.lines<-x.range[apply(invasion.grid.bistable,1,function(x){min(which(x))})]
	#	text(d.range[round(length(d.range)/2)],invasion.x.lines[round(length(d.range)/2)],s,col="red")
		my.x.bistable<-rbind(my.x.bistable,invasion.x.lines)
	}
	
	save(file=paste(directory,"bistable_invasion_grid.Robj",sep=""),my.x.bistable,my.s)
}

make.supp.fig.bistability.of.self.promoting.driver.w.add.cost<- function(){
	show(load(file=paste(directory,"bistable_invasion_grid.Robj",sep="")))
	par(mar=c(3,3,1,1))
	plot(range(d.range),c(10^(-6),1),type="n",log="y",xlab="",ylab="")
	mtext(side=1,line=2,"Drive coefficient, d")
	mtext(side=2,line=2,"Unstable equilibrium frequency")
		my.cols<-brewer.pal(length(my.s),name="Dark2")
	
	for(i in 1:length(my.s)){
		lines(d.range,my.x.bistable[i,],col=my.cols[i],lwd=2)
		
		sh<-my.s[i]
		s<-my.s[i]*2
		a1<-1-2*d.range*sh+4*d.range*s-2*d.range-3*sh;
		a2<-2*(s-sh)*(2*d.range-1)
	
		f.HWE.bistab<- (a1+sqrt(a1^2+8*sh*a2))/(2*a2)
		points(d.range[c(1:4,seq(5,200,by=5))],f.HWE.bistab[c(1:4,seq(5,200,by=5))],col=my.cols[i],lwd=2,lty=2,pch=19,cex=1)
	}
	
	
	legend(x="bottomleft",legend=paste("s=",2*my.s),col=my.cols,lty=1,lwd=2)
	
	dev.copy2eps(file=paste(directory,"bistable_x_vs_d_additive_s.eps",sep=""))
}


###Supplementary figure of male genotype control of female drive, and the effect of dominance.

make.gird.of.effect.male.geno.dominance<-function(){

	invasion.w.dom<-numeric()
	fixation.w.dom<-numeric()
	
	for(d.dom in seq(0,1,length=5)){
		print(d.dom)
		d.range<-seq(0.5,1,length=200)
		s3.range<-seq(0,1,length=200)
		
		invasion.grid.male.control<-sapply(s3.range,function(s3){
			sapply(d.range,function(d){
				test.invasion.self.prop.w.dom(d,d.dom=d.dom,s.het=0,s.hom=s3,x=.001,steps=10)
			})
		})
		
		fixation.grid.male.control<-sapply(s3.range,function(s3){
			sapply(d.range,function(d){
				test.invasion.self.prop.w.dom(d,d.dom=d.dom,s.het=0,s.hom=s3,x=.999,steps=10)
			})
		})
		
		invasion.lines.s3.cutoff.male<-s3.range[apply(invasion.grid.male.control,1,function(x){max(which(x))})]
		fixation.lines.s3.cutoff.male<-s3.range[apply(fixation.grid.male.control,1,function(x){max(which(x))})]
		
		invasion.w.dom<-rbind(invasion.w.dom,invasion.lines.s3.cutoff.male)
		fixation.w.dom<-rbind(fixation.w.dom,fixation.lines.s3.cutoff.male)
		
	}
	save(file=paste(directory,"effect_of_dominance_male_control.Robj",sep=""),invasion.w.dom,fixation.w.dom)
}

make.supp.fig.effect.male.geno.dominance<-function(){
	show(load(file=paste(directory,"effect_of_dominance_male_control.Robj",sep="")))
	my.cols<-brewer.pal(nrow(invasion.w.dom),name="Dark2")
	
	layout(matrix(1:4,nrow=2,byrow=TRUE))
	par(mar=c(3,3,1.5,0.1))
	for(i in c(1,3,4,5)){
		plot(c(0.5,1),c(0,0.5),type="n",ylab="",xlab="")
		lines(d.range,invasion.w.dom[i,],lty=1,col=my.cols[i])
		lines(d.range,fixation.w.dom[i,],lty=2,col=my.cols[i])
		if(i==3){ 
			legend("bottomright",lty=c(1,2),legend=c("invasion line","fixation line")); 
			mtext(side=3,line=0,"h=0.5")
			}	
		if(i==1){	
			mtext(side=2,line=2,"selection coeff against homozygotes, s",cex=0.8)
			mtext(side=3,line=0,"h=0")
			}
		if(i==4){
			mtext(side=2,line=2,"selection coeff against homozygotes, s",cex=0.8)
			text(0.67,0.185,"polymorphic",srt=50,col="black")
			mtext(side=3,line=0,"h=0.75")
			mtext(side=1,line=2,"Drive coefficient,d",cex=0.8)
	
			}
		if(i==5){
			text(0.78,0.3,"polymorphic",srt=50,col="black",cex=1)
			mtext(side=3,line=0,"h=1")		
			mtext(side=1,line=2,"Drive coefficient,d",cex=0.8)
		}
	}
	lines(d.range, (2*d.range-1)/(2*d.range),col=my.cols[i],lty=3)
	legend("bottomright",lty=3,legend=c("HWE invasion line"),col=my.cols[i]); 
	dev.copy2eps(file=paste(directory,"effect_of_dominance_on_invasion_space_four_graph.eps",sep=""))
}

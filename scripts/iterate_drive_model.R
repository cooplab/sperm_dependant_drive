iterate.1.locus.drive<-function(s.array,num.iterations,female.transmission.probs,my.geno.freqs,initialize.allele.freqs=c(0.999,0.001,0)){
	geno.freqs<-rep(NA,6)
	names(geno.freqs)<-geno.names

	my.freqs<-numeric()
	if(missing(my.geno.freqs)){
		print("running 2 allele system to stablity")
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

		geno.freqs.array<-outer(rep(geno.freqs,each=2),rep(geno.freqs,each=2),"*")

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


figures.of.traj.for.paper(){
	layout(t(1:2)); par(mar=c(3,3,1,1))
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
	text(x=125,y=0.9,"Driver invades")
	text(x=325,y=0.5,"Sperm-acting \n drive supressor \n invades")
	mtext("generations",side=1,line=2)
	mtext("Frequency",side=2,line=2)

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
	text(x=125,y=0.9,"Driver invades")
	text(x=310,y=0.8,"Sperm-acting \n driver & supressor \n invades")
	mtext("generations",side=1,line=2)
	mtext("Frequency",side=2,line=2)
	dev.copy2eps(file=paste(directory,"trajectories_of_sperm_based_supressors.eps",sep=""))
}

 		
 		
figs.of.sperm.vs.geno.dep<-function(){ 		
	pdf(file=paste(directory,"tmp.pdf",sep=""))
	#layout(1:2)
	###sperm allele that stops drive, but drives in females
	D<-matrix(1/2,nrow=2,ncol=3,dimnames=list(c("1 or 2","3"),c("12","13","23"))) ##entries are drive coeff of the 2nd allele listed against the 1st.
	D1 = 1
	D2 = 0.5
	s=1
	
	 D["1 or 2","12"]<-D1
	 D["1 or 2","13"]<-D1
	 D["3","12"]<-D2
	 D["3","13"]<-D2
	 
	run.iterations(D=D,s2=s,s3=s,sperm.or.geno.dependent="sperm",num.iter=100) #,plot.initial.rise=TRUE)
	#run.iterations(D=D,s2=s,s3=s,sperm.or.geno.dependent="geno",num.iter=100)
	
	###sperm allele that stops drive, and doesn't drive in females and is driven against by 3.
	D<-matrix(1/2,nrow=2,ncol=3,dimnames=list(c("1 or 2","3"),c("12","13","23")))
	 D["1 or 2","12"]<-D1
	 D["1 or 2","13"]<-D1
	 D["3","12"]<-D2
	 D["3","13"]<-D2
	 D["1 or 2","23"]<- 1-D1
	 D["3","23"]<-0.5
	
	run.iterations(D=D,s2=s,s3=0,sperm.or.geno.dependent="sperm",num.iter=100)
	run.iterations(D=D,s2=s,s3=0,sperm.or.geno.dependent="geno",num.iter=100)
	
	##female allele that blocks drive, and has no effect in sperm.
	D<-matrix(1/2,nrow=2,ncol=3,dimnames=list(c("1 or 2","3"),c("12","13","23")))
	 D["1 or 2","12"]<-D1
	 D["1 or 2","13"]<-0.5
	 D["3","12"]<-D1
	 D["3","13"]<-0.5
	 D["1 or 2","23"]<- 1-0.5
	 D["3","23"]<-0.5
	
	run.iterations(D=D,s2=s,s3=0,sperm.or.geno.dependent="sperm",num.iter=100,plot.initial.rise=TRUE)
	run.iterations(D=D,s2=s,s3=0,sperm.or.geno.dependent="geno",num.iter=100)
	
	dev.off()
}

examples.for.talk <-function(){
	D<-matrix(1/2,nrow=2,ncol=3,dimnames=list(c("1 or 2","3"),c("12","13","23"))) ##entries are drive coeff of the 2nd allele listed against the 1st.
	D1 = 1
	D2 = 0.5
	s=1
	
	 D["1 or 2","12"]<-D1
	 D["1 or 2","13"]<-D1
	 D["3","12"]<-D2
	 D["3","13"]<-D2
	 female.transmission.probs<-make.sperm.dep.female.transmission.prob(D)
	 pdf(file=paste(directory,"driver_to_initial_equilbrium.pdf",sep=""))
	par(mar=c(4,4,1,1))
	layout(1:2)
	 	old.geno.freqs<-iterate.1.locus.drive(s2=0,s3=0.2,num.iterations=100,female.transmission.probs=female.transmission.probs)
	plot.freqs(old.geno.freqs)
	old.geno.freqs<-iterate.1.locus.drive(s2=1,s3=1,num.iterations=100,female.transmission.probs=female.transmission.probs)
	plot.freqs(old.geno.freqs)
	dev.off()
	
	
	 pdf(file=paste(directory,"supressor_invades.pdf",sep=""))
	##female allele that blocks drive, and has no effect in sperm.
	D<-matrix(1/2,nrow=2,ncol=3,dimnames=list(c("1 or 2","3"),c("12","13","23")))
	 D["1 or 2","12"]<-D1
	 D["1 or 2","13"]<-0.5
	 D["3","12"]<-D1
	 D["3","13"]<-0.5
	 D["1 or 2","23"]<- 1-0.5
	 D["3","23"]<-0.5
	
	run.iterations(D=D,s2=1,s3=0,sperm.or.geno.dependent="sperm",num.iter=200,plot.initial.rise=TRUE)
	dev.off()
	
	 pdf(file=paste(directory,"better_driver_invades.pdf",sep=""))
	##female allele that blocks drive, and has no effect in sperm.
	D<-matrix(1/2,nrow=2,ncol=3,dimnames=list(c("1 or 2","3"),c("12","13","23")))
	 D["1 or 2","12"]<-0.8
	 D["1 or 2","13"]<-0.9
	 D["3","12"]<-0.8
	 D["3","13"]<-0.9
	 D["1 or 2","23"]<- 1-0.5
	 D["3","23"]<-0.5
	
	run.iterations(D=D,s2=1,s3=1,sperm.or.geno.dependent="sperm",num.iter=200,plot.initial.rise=TRUE)
	dev.off()
	
	
	 pdf(file=paste(directory,"green_beard_sperm_driver_invades.pdf",sep="")) ##shows it can't invade
	##female allele that blocks drive, and has no effect in sperm.
	D<-matrix(1/2,nrow=2,ncol=3,dimnames=list(c("1 or 2","3"),c("12","13","23")))
	 D["1 or 2","12"]<-0.7
	 D["1 or 2","13"]<-0.7
	 D["3","12"]<-.7
	 D["3","13"]<-1
	 D["1 or 2","23"]<- 1-0.5
	 D["3","23"]<-0.5
	
	run.iterations(D=D,s2=.3,s3=.3,sperm.or.geno.dependent="sperm",num.iter=200,plot.initial.rise=TRUE)
	dev.off()
	
	
	 pdf(file=paste(directory,"driver_which_avoids_itself_invades.pdf",sep=""))
	D<-matrix(1/2,nrow=2,ncol=3,dimnames=list(c("1 or 2","3"),c("12","13","23"))) ##entries are drive coeff of the 2nd allele listed against the 1st.
	D1 = 1
	D2 = 0.5
	s.array<-rep(0,6)
	names(s.array)<-c("11","12","13","23","22","33")
	s.array["23"]<- 1; s.array["33"]<- 1; s.array["22"]<- 1
	
	 D["1 or 2","12"]<-D1
	 D["1 or 2","13"]<-D1
	 D["3","12"]<-D2
	 D["3","13"]<-D2
	 
	run.iterations(D=D,s.array=s.array,sperm.or.geno.dependent="sperm",num.iter=200,plot.initial.rise=TRUE)
	dev.off()
	
	###sperm allele that stops drive, and doesn't drive in females and is driven against by 3.
	 pdf(file=paste(directory,"driver_then_sperm_supressor.pdf",sep=""))
	D1 = 1
	D2 = 0.5
	
	D<-matrix(1/2,nrow=2,ncol=3,dimnames=list(c("1 or 2","3"),c("12","13","23")))
	 D["1 or 2","12"]<-D1
	 D["1 or 2","13"]<-D1
	 D["3","12"]<-D2
	 D["3","13"]<-D2
	 D["1 or 2","23"]<- 1-D1
	 D["3","23"]<-0.5
	
	run.iterations(D=D,s2=s,s3=0,sperm.or.geno.dependent="sperm",num.iter=200,plot.initial.rise=TRUE)
	dev.off()
	
#	pdf(file=paste(directory,"sperm_promotor_driver_invades.pdf",sep=""))

	 D["1 or 2","12"]<-0.5
	 D["1 or 2","13"]<-0.5
	 D["3","12"]<-.5
	 D["3","13"]<-1
	 D["1 or 2","23"]<- 1-0.5
	 D["3","23"]<-0.5
		 female.transmission.probs<-make.sperm.dep.female.transmission.prob(D)
	par(mar=c(4,4,1,1))
	old.geno.freqs<-iterate.1.locus.drive(s2=0,s3=0,num.iterations=10000,female.transmission.probs=female.transmission.probs,initialize.allele.freqs =c(0.99, 0,0.01 ))
	plot.freqs(old.geno.freqs)



		 D["1 or 2","12"]<-0.5
	 D["1 or 2","13"]<-0.5
	 D["3","12"]<-.5
	 D["3","13"]<-0.6
	 D["1 or 2","23"]<- 1-0.5
	 D["3","23"]<-0.5
		 female.transmission.probs<-make.male.geno.dep.female.transmission.prob(D)
	par(mar=c(4,4,1,1))
	old.geno.freqs<-iterate.1.locus.drive(s2=0,s3=.1,num.iterations=3000,female.transmission.probs=female.transmission.probs,initialize.allele.freqs =c(0.99, 0,0.01 ))
	plot.freqs(old.geno.freqs)

	
}

##
x=.01;my.x2<-numeric()
for(i in 1:4000){
x_new = (1-s)*x^2+2*x*(1-x)*(1/4 + 1/2*D1 )
x_new= x_new/(1-s*x^2)
my.x2<-c(my.x2,x)
x=x_new
}
#lines(my.x2,lty=2,col="red")
abline(h=my.x2[4000],lty=2,col="red")

###poor man's iteration of sperm dep. using HWE
x=.01;my.x3<-numeric()
for(i in 1:4000){
x_new = (1-s)*x^2+2*x*(1-x)*(1/4 + 1/2*(x*D2+(1-x)*D1 ))
x_new= x_new/(1-s*x^2)
my.x3<-c(my.x3,x)
x=x_new
}
#lines(my.x3,lty=2,col="green")
abline(h=my.x3[4000],lty=2,col="green")

a= s
b = -((1-s)+2*(1/4+D1/2))
c=2*(1/4+D1/2)-1
(-b + sqrt(b^2-4*a*c))/(2*a)
(-b - sqrt(b^2-4*a*c))/(2*a)
 
 abline(h=(-b - sqrt(b^2-4*a*c))/(2*a))

  HWEx2<-numeric()
  HWEx3<-numeric()
 x2 = .0 #my.x2[4000]*.999
 x3 = .01
 s2=s
s3=s 
 for(i in 1:1000){

 

 wbar<-(1 - s2*x2^2 - s3*x3^2 - s3*2*x2*x3)
 
 x2_new = (1-s2)*x2^2 +(2/2)*(1-s3)*x2*x3+
 2*x2*(1-x2-x3)*(1/4+(1/2)*(x3*D2+(1-x3)*D1))
 
 x3_new = (1-s2)*x3^2 +(2/2)*(1-s3)*x2*x3+
 2*x3*(1-x2-x3)*(1/4+(1/2)*(x3*D2+(1-x3)*D1))
 
 x2=x2_new/wbar
 x3=x3_new/wbar
 HWEx2<-c(HWEx2,x2)
 HWEx3<-c(HWEx3,x3)
 }
 

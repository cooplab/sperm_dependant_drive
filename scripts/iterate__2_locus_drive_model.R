library(abind)

iterate.2.locus.drive<-function(s2,s3,num.iterations,female.transmission.probs,my.geno.freqs){
 	geno.names<-c("11", "12","13", "22", "23", "33")
 	x<- rgamma(18,1)
	geno.freqs<-array(x/sum(x),dim=c(6,3),dimnames=list(geno.names,c("0","1","2")))
	transmitted.allele<-c(1,1,1,2,1,3,2,2,2,3,3,3)
	my.freqs<-numeric()
#	if(missing(my.geno.freqs)){
#		print("running 2 allele system to stablity")
#		allele.freqs<-c(0.999,0.001,0)
#		names(allele.freqs)<-c("1","2","3")	

	
#		geno.freqs["11"]<-allele.freqs["1"]^2
#		geno.freqs["22"]<-allele.freqs["2"]^2
#		geno.freqs["33"]<-allele.freqs["3"]^2

#		geno.freqs["12"]<-2*allele.freqs["1"]*allele.freqs["2"]
#		geno.freqs["13"]<-2*allele.freqs["1"]*allele.freqs["3"]
#		geno.freqs["23"]<-2*allele.freqs["3"]*allele.freqs["2"]

#	}else{
#		print("initiating from my.geno.freqs")
#		geno.freqs[c("11","12","22")]<-my.geno.freqs[["geno.freqs"]][nrow(my.geno.freqs[["geno.freqs"]]),c("11","12","22")]*.99
#		geno.freqs[c("13","23")]<-my.geno.freqs[["allele.freqs"]][nrow(my.geno.freqs[["geno.freqs"]]),1:2]*.99*.01
#		geno.freqs["33"]<-.01^2
		
#		}

female.geno.freqs.array<-array(NA,dim=c(12,12,3),dimnames=list(rep(geno.names,each=2),rep(geno.names,each=2),c("0","1","2")))
	for(i in 1:num.iterations){
	#selection
		geno.freqs["22",]<-geno.freqs["22",] * (1-s2)
		geno.freqs[c("33","23"),]<-geno.freqs[c("33","23"),] * (1-s3)
		geno.freqs<-geno.freqs/sum(geno.freqs)	
		female.geno.freqs.array<-array(NA,dim=c(12,12,3),dimnames=list(rep(geno.names,each=2),rep(geno.names,each=2),c("0","1","2")))

		female.geno.freqs.array[,,"0"]<-matrix( rep(rep(geno.freqs[,"0"],each=2,2),6),nrow= 12,byrow=TRUE)
		female.geno.freqs.array[,,"1"]<-matrix( rep(rep(geno.freqs[,"1"],each=2,2),6),nrow= 12,byrow=TRUE)
		female.geno.freqs.array[,,"2"]<-matrix( rep(rep(geno.freqs[,"2"],each=2,2),6),nrow= 12,byrow=TRUE)

		##transmission 
		#
		##male geno. multiples through by female trans
	
		one.trans<-female.geno.freqs.array*female.transmission.probs*0.5	
		trans<-array(NA,dim=c(12,12,3,3), dimnames=list(rep(geno.names,each=2),rep(geno.names,each=2),c("0","1","2"),c("0","1","2")) )
		trans[,,,"0"]<-one.trans*abind(t(female.geno.freqs.array[,,"0"]),t(female.geno.freqs.array[,,"0"]),t(female.geno.freqs.array[,,"0"]),along=3)  
		
		trans[,,,"1"]<-one.trans*abind(t(geno.freqs.array[,,"1"]),t(geno.freqs.array[,,"1"]),t(geno.freqs.array[,,"1"]),along=3)
		trans[,,,"2"]<-one.trans*abind(t(geno.freqs.array[,,"2"]),t(geno.freqs.array[,,"2"]),t(geno.freqs.array[,,"2"]),along=3)
	

		new.geno.freqs<-array(NA,dim=c(6,3),dimnames=list(geno.names,c("0","1","2")))	
		
		for(loc1.1 in 1:3){
			for(loc1.2 in loc1.1:3){
				for(loc2 in c("0","1","2")){
				
					if(loc1.1 == loc1.2){
						
						new.geno.freqs[paste(loc1.1,loc1.2,sep=""),loc2]<-trans.prob.2nd.locus(trans,make.geno=loc2,i=loc1.1,j=loc1.2)
					}
													
					if(loc1.1 != loc1.2){ 
						
						new.geno.freqs[paste(loc1.1,loc1.2,sep=""),loc2] <-	  trans.prob.2nd.locus(trans,make.geno=loc2,i=loc1.1,j=loc1.2)
														+ trans.prob.2nd.locus(trans,make.geno=loc2,i=loc1.2,j=loc1.1)
					}
					}
				}
			}
													
												

		recover()
		geno.freqs<-new.geno.freqs/sum(new.geno.freqs)


		#my.freqs<-rbind(my.freqs,geno.freqs)

	}

	freq.1<-my.freqs[,"11"]+0.5*rowSums(my.freqs[,c("12","13")])
	freq.2<-my.freqs[,"22"]+0.5*rowSums(my.freqs[,c("12","23")])
	freq.3<-my.freqs[,"33"]+0.5*rowSums(my.freqs[,c("13","23")])


	geno.and.allele.freqs<-list()
	geno.and.allele.freqs[["geno.freqs"]]<-my.freqs
	geno.and.allele.freqs[["allele.freqs"]]<-cbind(freq.1,freq.2,freq.3)

	return(geno.and.allele.freqs) 
 }
 
 
 trans.prob.2nd.locus<-function(trans,make.geno,i,j){
 		transmitted.allele<-c(1,1,1,2,1,3,2,2,2,3,3,3)
	if(make.geno==0){
		freq<- sum(trans[transmitted.allele==i,transmitted.allele==j,"0","0"]
			+0.5*trans[transmitted.allele==i,transmitted.allele==j,"0","1"]
			+0.5*trans[transmitted.allele==i,transmitted.allele==j,"1","0"]
			+0.25*trans[transmitted.allele==i,transmitted.allele==j,"1","1"])
	 }
	 
	if(make.geno==1){
		freq<- 
		sum(0.5*trans[transmitted.allele==i,transmitted.allele==j,"1","1"]
		 +0.5*trans[transmitted.allele==i,transmitted.allele==j,"2","1"]
		 +0.5*trans[transmitted.allele==i,transmitted.allele==j,"1","2"]
		+0.5*trans[transmitted.allele==i,transmitted.allele==j,"0","1"]
		+0.5*trans[transmitted.allele==i,transmitted.allele==j,"1","0"]
		+trans[transmitted.allele==i,transmitted.allele==j,"2","0"]
		+trans[transmitted.allele==i,transmitted.allele==j,"0","2"])
	 }
	  
	  if(make.geno==2){
		 freq<-sum(trans[transmitted.allele==i,transmitted.allele==j,"2","2"]
		 +0.5*trans[transmitted.allele==i,transmitted.allele==j,"2","1"]
		 +0.5*trans[transmitted.allele==i,transmitted.allele==j,"1","2"]
		 +0.25*trans[transmitted.allele==i,transmitted.allele==j,"1","1"])
	 }
	 
	 return(freq)
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
	
transmission.probs
}


plot.freqs<-function(new.geno.freqs,plot.HWE.dev=FALSE){
	plot(new.geno.freqs[["allele.freqs"]][,1],ylim=c(0,1),type="l")
 	lines(new.geno.freqs[["allele.freqs"]][,2],col="red")
 	lines(new.geno.freqs[["allele.freqs"]][,3],col="green")

	if(plot.HWE.dev){
#####deviations away from HWe 
 		plot(new.geno.freqs[["geno.freqs"]][,"11"] / new.geno.freqs[["allele.freqs"]][,1]^2,ylim=c(.85,1.15),type="l",col="red")
		lines(new.geno.freqs[["geno.freqs"]][,"22"] / new.geno.freqs[["allele.freqs"]][,2]^2,col="blue")
		lines(new.geno.freqs[["geno.freqs"]][,"33"] /new.geno.freqs[["allele.freqs"]][,3]^2,col="green")

		lines(new.geno.freqs[["geno.freqs"]][,"12"] / (2*new.geno.freqs[["allele.freqs"]][,1]*new.geno.freqs[["allele.freqs"]][,2]),col="red",lty=2)
		lines(new.geno.freqs[["geno.freqs"]][,"13"] / (2*new.geno.freqs[["allele.freqs"]][,1]*new.geno.freqs[["allele.freqs"]][,3]),col="blue",lty=2)
		lines(new.geno.freqs[["geno.freqs"]][,"23"] / (2*new.geno.freqs[["allele.freqs"]][,3]*new.geno.freqs[["allele.freqs"]][,2]),col="green",lty=2)
	}

}


run.iterations<-function(D,s2,s3,sperm.or.geno.dependent,num.iter=4000){
	if(sperm.or.geno.dependent=="sperm"){ female.transmission.probs<-make.sperm.dep.female.transmission.prob(D)}
	if(sperm.or.geno.dependent=="geno"){ female.transmission.probs<-make.male.geno.dep.female.transmission.prob(D)}
	
	old.geno.freqs<-iterate.1.locus.drive(s2=s2,s3=s3,num.iterations=num.iter,female.transmission.probs=female.transmission.probs)
	new.geno.freqs<-iterate.1.locus.drive(s2=s2,s3=s3,num.iterations=num.iter,female.transmission.probs=female.transmission.probs,my.geno.freqs=old.geno.freqs)
	plot.freqs(new.geno.freqs)
}

##to get this up and running
female.transmission.probs<-make.male.geno.dep.female.transmission.prob(D)
female.transmission.probs<-abind(female.transmission.probs,female.transmission.probs,female.transmission.probs,along=3)
iterate.2.locus.drive(s2=0,s3=0,num.iterations=200,female.transmission.probs=female.transmission.probs,my.geno.freqs=5)



pdf(file="tmp.pdf")
layout(1:2)
###sperm allele that stops drive, but drives in females
D<-matrix(1/2,nrow=2,ncol=3,dimnames=list(c("1 or 2","3"),c("12","13","23"))) ##entries are drive coeff of the 2nd allele listed against the 1st.
D1 = 1
D2 = 0.5
s=1

 D["1 or 2","12"]<-D1
 D["1 or 2","13"]<-D1
 D["3","12"]<-D2
 D["3","13"]<-D2
 
run.iterations(D=D,s2=s,s3=s,sperm.or.geno.dependent="sperm",num.iter=100)
run.iterations(D=D,s2=s,s3=s,sperm.or.geno.dependent="geno",num.iter=100)

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

run.iterations(D=D,s2=s,s3=0,sperm.or.geno.dependent="sperm",num.iter=100)
run.iterations(D=D,s2=s,s3=0,sperm.or.geno.dependent="geno",num.iter=100)

dev.off()


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
 

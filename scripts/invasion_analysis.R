library("RColorBrewer")

test.invasion.self.prop<-function(d,s.het,s.hom,x=.001,steps=10,sperm.or.male="sperm"){
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
	old.geno.freqs<-iterate.1.locus.drive(s.array,num.iterations=steps,female.transmission.probs=female.transmission.probs,initialize.allele.freqs =c(1-x, 0,x))
#recover()
	invading<-old.geno.freqs[["allele.freqs"]][steps,3]>old.geno.freqs[["allele.freqs"]][1,3]
	return(invading)
}

####invasion test for male dependent drive w. dominance controled by d.dom
test.invasion.self.prop.w.dom<-function(d,d.dom,s.het,s.hom,x=.001,steps=10){
	s.array<-rep(0,6)
	names(s.array)<-c("11","12","13","23","22","33")
	s.array["13"]<- s.het
	s.array["33"]<- s.hom
	
		D<-matrix(1/2,nrow=3,ncol=3,dimnames=list(c("1 or 2","3.het","3"),c("12","13","23"))) ##entries are drive coeff of the 2nd 
 ##entries are drive coeff of the 2nd 
	D["3","13"]<-d
	D["3.het","13"]<-(d-0.5)*d.dom + 0.5  ##d.dom modifies the effect of drive 
	D["1 or 2","23"]<- 1-0.5
	D["3","23"]<-0.5
	
	female.transmission.probs<-make.male.geno.dep.female.transmission.prob.w.dom(D) 
	old.geno.freqs<-iterate.1.locus.drive(s.array,num.iterations=steps,female.transmission.probs=female.transmission.probs,initialize.allele.freqs =c(1-x, 0,x))
#recover()
	invading<-old.geno.freqs[["allele.freqs"]][steps,3]>old.geno.freqs[["allele.freqs"]][1,3]
	return(invading)
}



test.fixation.of.simple.driver<-function(d,s.het,s.hom,x=.999,steps=10){
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
	old.geno.freqs<-iterate.1.locus.drive(s.array=s.array,num.iterations=steps,female.transmission.probs=female.transmission.probs,initialize.allele.freqs =c(1-x, x,0))
	fixing<-old.geno.freqs[["allele.freqs"]][steps,2]>old.geno.freqs[["allele.freqs"]][1,2]
	return(fixing)
}


###invasion of self-driver against homozy. cost.
	D<-matrix(1/2,nrow=2,ncol=3,dimnames=list(c("1 or 2","3"),c("12","13","23"))) ##entries are drive coeff of the 2nd allele listed against the 1st.

d.range<-seq(0.5,1,length=200)
s3.range<-seq(0,1,length=200)

invasion.grid<-sapply(s3.range,function(s3){
	sapply(d.range,function(d){
		test.invasion.self.prop(d=d,s.het=0,s.hom=s3)
	})
})

invasion.lines.s3.cutoff<-s3.range[apply(invasion.grid,1,function(x){max(which(x))})]

image(d.range, s3.range,invasion.grid,xlab="Drive coeff. D",ylab="Homozygous cost, s")
theory.s3.cutoff<-d.range-0.5
lines(d.range,theory.s3.cutoff,col="blue")

invasion.grid.HWE<-sapply(s3.range,function(s3){
	sapply(d.range,function(d){
		invading <- (.001*(0.5+s3-d)+(-0.5-s3+d))>0
		invading
	})
})

invasion.lines.s3.HWE.cutoff<-s3.range[apply(invasion.grid.HWE,1,function(x){max(which(x))})]
lines(d.range,invasion.lines.s3.HWE.cutoff)


drive.coeffs<-c(0.6,.7,.8)
selection.coeffs<-rep(0.1,3)
my.point.cols<-c("red","blue","orange")
layout(t(1:2))
image(d.range, s3.range,invasion.grid,xlab="Drive coeff. D",ylab="Homozygous cost, s")
theory.s3.cutoff<-d.range-0.5
lines(d.range,theory.s3.cutoff)

points(drive.coeffs,selection.coeffs,col=my.point.cols,pch=19)
plot(c(0,10000),c(0,1),type="n",xlab="time, generations",ylab="frequency")
for(i in 1:3){
		 D["1 or 2","12"]<-0.5
	 D["1 or 2","13"]<-0.5
	 D["3","12"]<-.5
	 D["3","13"]<-drive.coeffs[i]
	 D["1 or 2","23"]<- 1-0.5
	 D["3","23"]<-0.5
		 female.transmission.probs<-make.male.geno.dep.female.transmission.prob(D)
	par(mar=c(4,4,1,1))
	old.geno.freqs<-iterate.1.locus.drive(s2=0,s3=selection.coeffs[i],num.iterations=10000,female.transmission.probs=female.transmission.probs,initialize.allele.freqs =c(0.999, 0,0.001 ))
	lines(old.geno.freqs[["allele.freqs"]][,3],col=my.point.cols[i],lwd=2)
	}


d.range<-seq(0.5,1,length=100)
s3.range<-seq(0,1,length=200)

invasion.grid<-sapply(s3.range,function(s3){
	sapply(d.range,function(d){
		test.invasion(d=d,s2=0,s3=s3)
	})
})


###testing fixation of simple driver###############
fixation.grid<-sapply(s3.range,function(s3){
	sapply(d.range,function(d){
		test.fixation.of.simple.driver(d=d,s2=s3,s3=s3)
	})
})

fixation.grid.hwe<-sapply(s3.range,function(s3){
	sapply(d.range,function(D1){

		s.homozyg=s3;d=D1;a=s.homozyg;b=.5-s.homozyg-d;c=d-.5
		quad.soln<- (-b +c(-1,1)*sqrt(b^2-4*a*c))/(2*a)
		return(!any(quad.soln<0.999 & quad.soln >0.001))
	})
})
image(d.range,s3.range,fixation.grid)
image(d.range,s3.range,fixation.grid.hwe)


fixation.lines.s3<-s3.range[apply(fixation.grid,1,function(x){max(which(x))})]
lines(d.range,fixation.lines.s3)
fixation.lines.s3.hwe<-s3.range[apply(fixation.grid.hwe,1,function(x){max(which(x))})]
lines(d.range,fixation.lines.s3)


####bistable invasion of driver against heterozyg. cost
#pdf(file=paste(directory,"bistable_x_vs_d_additive_s.pdf",sep=""))
d.range<-seq(0.5,1,length=200)
x.range<- 10^seq(-7,0,length=200)  #seq(0.001,.99,length=50)

my.s<-c(0.01,0.001,0.0001,0.00001)

if(FALSE){
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
#dev.off()

##############Phase diagram figure for the paper!!!!
if(FALSE){
	D<-matrix(1/2,nrow=2,ncol=3,dimnames=list(c("1 or 2","3"),c("12","13","23"))) ##entries are drive coeff of the 2nd allele listed against the 1st.
	
	d.range<-seq(0.5,1,length=200)
	s3.range<-seq(0,1,length=200)
	
	invasion.grid<-sapply(s3.range,function(s3){
		sapply(d.range,function(d){
			test.invasion.self.prop(d=d,s.het=0,s.hom=s3)
		})
	})
	
	invasion.lines.s3.cutoff<-s3.range[apply(invasion.grid,1,function(x){max(which(x))})]
	
	fixation.grid<-sapply(s3.range,function(s3){
		sapply(d.range,function(d){
			test.fixation.of.simple.driver(d=d,s.het=0,s.hom=s3)
		})
	})
	fixation.lines.s3<-s3.range[apply(fixation.grid,1,function(x){max(which(x))})]
save(file=paste(directory,"invasion_grid_homozy_cost.Robj",sep=""),fixation.lines.s3,fixation.grid,invasion.lines.s3.cutoff,invasion.grid)
}

show(load(paste(directory,"invasion_grid_homozy_cost.Robj",sep="")))

my.cols<-brewer.pal(2,name="Dark2")
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
dev.copy2eps(file=paste(directory,"invasion_space_recessive_driver.eps",sep=""))
##############Phase diagram figure showing bistability 
##done BATCH in homozyg_bistability.R
#
d.range<-seq(0.5,1,length=200)
s3.range<-seq(0,1,length=200)
if(FALSE){
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

load("homozyg_bistability.Robj")
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

#######################################################
#####invasion space for male genotype controlled allele

##fully dominant allele
d.range<-seq(0.5,1,length=200)
s3.range<-seq(0,1,length=200)

invasion.grid.male.control<-sapply(s3.range,function(s3){
	sapply(d.range,function(d){
		test.invasion.self.prop(d=d,s.het=0,s.hom=s3,sperm.or.male="male")
	})
})

fixation.grid.male.control<-sapply(s3.range,function(s3){
	sapply(d.range,function(d){
		test.invasion.self.prop(d=d,s.het=0,s.hom=s3,x=0.999,sperm.or.male="male")
	})
})

invasion.lines.s3.cutoff.male<-s3.range[apply(invasion.grid.male.control,1,function(x){max(which(x))})]
fixation.lines.s3.cutoff.male<-s3.range[apply(fixation.grid.male.control,1,function(x){max(which(x))})]

plot(d.range,invasion.lines.s3.cutoff.male,col="red")
lines(d.range,fixation.lines.s3.cutoff.male,col="green")

if(FALSE){
pdf(file=paste(directory,"effect_of_dominance_on_invasion_space.pdf",sep=""))

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
	
	plot(d.range,invasion.lines.s3.cutoff.male,col="red",ylim=c(0,1),main=paste("dominance of male genotype = ",d.dom),type="l",lwd=2)
	lines(d.range,fixation.lines.s3.cutoff.male,col="green",lwd=2)
}
dev.off()
save(file=paste(directory,"effect_of_dominance_male_control.Robj",sep=""),invasion.w.dom,fixation.w.dom)
}
load(file=paste(directory,"effect_of_dominance_male_control.Robj",sep=""))
my.cols<-brewer.pal(nrow(invasion.w.dom),name="Dark2")
plot(c(0.5,1),c(0,0.5),type="n",,xlab="Drive coefficient",ylab="selection coeff against homozygotes")
sapply(1:nrow(invasion.w.dom),function(i){lines(d.range,invasion.w.dom[i,],lty=1,col=my.cols[i])})
sapply(1:nrow(invasion.w.dom),function(i){lines(d.range,fixation.w.dom[i,],lty=2,col=my.cols[i])})
legend("topleft",lty=c(1,2,rep(NA,nrow(invasion.w.dom))),pch=c(NA,NA,rep(19,nrow(invasion.w.dom))),col=c(rep("black",2),my.cols),legend=c("invasion line","fixation line",paste("dominance = ",seq(0,1,length=5))))
dev.copy2eps(file=paste(directory,"effect_of_dominance_on_invasion_space_one_graph.eps",sep=""))

text(0.6,0.45,expression(d[h]== 0))

layout(matrix(1:4,nrow=2,byrow=TRUE))
par(mar=c(3,3,1.5,0.1))
for(i in c(1,3,4,5)){
	plot(c(0.5,1),c(0,0.5),type="n",ylab="",xlab="")
	lines(d.range,invasion.w.dom[i,],lty=1,col=my.cols[i])
	lines(d.range,fixation.w.dom[i,],lty=2,col=my.cols[i])
	if(i==3){ 
		legend("bottomright",lty=c(1,2),legend=c("invasion line","fixation line")); 
		mtext(side=3,line=0,"h=0.5")

		#text(0.6,0.45,expression(paste(d[h],"=0.5d")))

		}	
	if(i==1){	
		mtext(side=2,line=2,"selection coeff against homozygotes, s",cex=0.8)
		mtext(side=3,line=0,"h=0")
		#text(0.6,0.45,expression(paste(d[h],"=0")))
		}
	if(i==4){
		mtext(side=2,line=2,"selection coeff against homozygotes, s",cex=0.8)
		text(0.67,0.185,"polymorphic",srt=50,col="black")
		mtext(side=3,line=0,"h=0.75")
		#text(0.6,0.45,expression(paste(d[h],"=0.75d")))
		mtext(side=1,line=2,"Drive coefficient,d",cex=0.8)

		}
	if(i==5){
		text(0.78,0.3,"polymorphic",srt=50,col="black",cex=1)
#		text(0.6,0.45,expression(d[h]== d))
		mtext(side=3,line=0,"h=1")		
		mtext(side=1,line=2,"Drive coefficient,d",cex=0.8)

	}
}
lines(d.range, (2*d.range-1)/(2*d.range),col=my.cols[i],lty=3)
legend("bottomright",lty=3,legend=c("HWE invasion line"),col=my.cols[i]); 
dev.copy2eps(file=paste(directory,"effect_of_dominance_on_invasion_space_four_graph.eps",sep=""))


s.array<-rep(0,6)
	names(s.array)<-c("11","12","13","23","22","33")
	s.array["33"]<- 0.35
	
	
	 D["1 or 2","12"]<-0.5
	 D["1 or 2","13"]<-0.5
	 D["3","12"]<-.5
	 D["3","13"]<-0.8
	 D["1 or 2","23"]<- 1-0.5
	 D["3","23"]<-0.5
	female.transmission.probs<-make.male.geno.dep.female.transmission.prob(D) 
old.geno.freqs<-iterate.1.locus.drive(s.array=s.array,num.iterations=4000,female.transmission.probs=female.transmission.probs,initialize.allele.freqs =c(0.99, 0,0.01 ))
plot(old.geno.freqs[["allele.freqs"]][,3],col="green")
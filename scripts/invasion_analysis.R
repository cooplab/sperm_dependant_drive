

test.invasion<-function(d,s2,s3,x=.001,steps=10){
	D["1 or 2","12"]<-0.5
	D["1 or 2","13"]<-0.5
	D["3","12"]<-.5
	D["3","13"]<-d
	D["1 or 2","23"]<- 1-0.5
	D["3","23"]<-0.5
	female.transmission.probs<-make.sperm.dep.female.transmission.prob(D)
	par(mar=c(4,4,1,1))
	old.geno.freqs<-iterate.1.locus.drive(s2=s2,s3=s3,num.iterations=steps,female.transmission.probs=female.transmission.probs,initialize.allele.freqs =c(1-x, 0,x))

	invading<-old.geno.freqs[["allele.freqs"]][steps,3]>old.geno.freqs[["allele.freqs"]][1,3]
	return(invading)
}


###invasion of self-driver against homozy. cost.

d.range<-seq(0.5,1,length=100)
s3.range<-seq(0,1,length=200)

invasion.grid<-sapply(s3.range,function(s3){
	sapply(d.range,function(d){
		test.invasion(d=d,s2=0,s3=s3)
	})
})

invasion.lines.s3.cutoff<-s3.range[apply(invasion.grid,1,function(x){max(which(x))})]

image(d.range, s3.range,invasion.grid,xlab="Drive coeff. D",ylab="Homozygous cost, s")
theory.s3.cutoff<-d.range-0.5
lines(d.range,theory.s3.cutoff)

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

####bistable invasion of driver against heterozyg. cost


s=.01

d.range<-seq(0.5,1,length=10)
x.range<-seq(0.001,.99,length=20)

invasion.grid.bistable<-sapply(x.range,function(x){
	sapply(d.range,function(d){
		test.invasion(d=d,s2=s,s3=0,x=x)
	})
})



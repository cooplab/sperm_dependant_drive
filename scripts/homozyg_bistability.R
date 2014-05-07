 source("iterate_drive_model.R")
 

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


#approx.bistable.freq<-(1-2*d+4*d*s)/(2*s*(2*d-1))
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
	
	
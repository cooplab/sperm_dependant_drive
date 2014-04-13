

#SIMPLE ITERATION IN TERMS OF GENOS

nextGen = function(x.last,d,w){
	makeZygotes = function(geno.freqs,d){
		A.sperm = sum(geno.freqs*c(1,1/2,0))
		B.sperm = sum(geno.freqs*c(0,1/2,1))
		AA.kids = sum(c(geno.freqs["AA"]*A.sperm  , geno.freqs["AB"]*A.sperm/2))
		AB.kids = sum(c(geno.freqs["AA"]*B.sperm  , geno.freqs["AB"]* B.sperm*(1-d)))
		BA.kids = sum(c(geno.freqs["AB"]*A.sperm/2  , geno.freqs["BB"]* A.sperm))
		BB.kids = sum(c(geno.freqs["AB"]*B.sperm*d  , geno.freqs["BB"]* B.sperm))
		c(AA = AA.kids, AB = (AB.kids + BA.kids), BB = BB.kids) 
	}
	indSel = function(geno.freqs,w){
		W = sum(geno.freqs*w)
		apply(cbind(geno.freqs,w),1,prod)/W
	}
	after.drive = makeZygotes(x.last,d=d)
	after.sel   = indSel(after.drive,w=w) 
	return(after.sel)
}




#ITERATE THIS OR PUT IN WHILE LOOP
#nextGen(x.last = c(AA = 1/4, AB = 1/2, BB = 1/4), d=1, w = c(AA=1, AB=1, BB=1))

iterate = function(gen,fB0,d,s){
	temp = numeric(gen)
	temp[1] = fB0
	genos = c(AA = (1-fB0)^2, AB = 2*(fB0)*(1-fB0), BB = (fB0)^2)
	for(g in seq_along(temp)[-1]){
		genos = nextGen(x.last = genos, d=d, w = c(AA=1, AB=1, BB=1-s))
		#print(c(genos,sum(genos)))
		temp[g] = sum(genos*c(0,1/2,1))
	}
	temp[g]
}

k = sapply(seq(0.5,1,.005),function(D){ 
		sapply(seq(0,1,.01),function(S){
			print(paste(D,S))
			iterate(gen = 50,fB0=10^(-5), d = D, s = S)  
	})
})

image(x=seq(0.5,1,.005),y=seq(0,1,.01),t(sign(k-1*10^(-5))),xlab="d",ylab="s",main="Invasion criteria for recessive-lethal, self promoting driver")
curve(  (2*x-1)/(4*x), .5, 1,add=TRUE,lwd=2)
text(.8,.05,"Invasion of self promotor")
text(.6,.8,"Self-promotor cannot invade")
text(.7,.3,"s = 1/2 - 1/(4d)")
arrows(.7,.25,.8, 0.1875,length=.15,lwd=3)

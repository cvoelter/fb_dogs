library(lme4)
lmer.warns<-function(x){
	n.opt.warnings=paste(unlist(x@optinfo$conv$lme4), collapse="@")
	n.fun.warnings=paste(unlist(x@optinfo$warnings), collapse="@")
	return(c(n.opt.warnings=n.opt.warnings, n.fun.warnings=n.fun.warnings))
}
source("/home/roger/r_functions/drop1_para.r")

xdata=data.frame(expand.grid(belief=c("false", "true"), sex=c("F", "M"), baited.fst=c("blue", "grey")))
m.mat=model.matrix(object=~belief+sex+baited.fst, data=xdata)
fe2sim=rep(0, ncol(m.mat))
names(fe2sim)=colnames(m.mat)
fe2sim=c(fe2sim, age=0)
re2sim=fe2sim
fe2sim["(Intercept)"]=qlogis(0.85)
fe2sim["belieftrue"]=qlogis(0.55)-qlogis(0.85)
re2sim["(Intercept)"]=0.5
re2sim["belieftrue"]=0.5

n.per.group=seq(3, 25, by=1)
age.range.t=c(5/12, 13)
age.range.c=c(6, 10)
prob.outs.c=0.4
n.breeds=8
save.image("/home/roger/roger/2020/vienna/done/lucrezia_sim_bw_gr/power_w_one_tr_per_dog_new_feb.RData")
load("/home/roger/roger/2020/vienna/done/lucrezia_sim_bw_gr/power_w_one_tr_per_dog_new_feb.RData")


load("/home/mundryr/transfer/power_w_one_tr_per_dog_new_feb.RData")
library(lme4)
library(parallel)
cl <- makeCluster(getOption("cl.cores", detectCores()))
cl=cl[1:(length(cl)-1)]
parLapply(X=1:length(cl), cl=cl, function(x){library(lme4)})

n.sims=1000
contr=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=10000))
for(n in 1:length(n.per.group)){
	test.data.o=xdata[rep(1:nrow(xdata), each=n.per.group[n]), ]
	clusterExport(cl=cl, 
		varlist=c("test.data.o", "fe2sim", "re2sim", "n.per.group", "age.range.t", "age.range.c", "prob.outs.c", "n.breeds",
			"contr", "drop1p", "lmer.warns"))
	i.res=parLapply(X=1:n.sims, cl=cl, function(s){
		test.data=test.data.o
		test.data$age=runif(n=nrow(test.data), min=min(age.range.c), max=max(age.range.c))
		sel=rbinom(n=nrow(test.data), size=1, prob=prob.outs.c)==1
		age.mod=c(runif(n=sum(sel), min=min(age.range.t), max=min(age.range.c)), runif(n=sum(sel), min=max(age.range.c), max=max(age.range.t)))
		test.data$age[sel]=sample(age.mod, size=sum(sel), replace=F)
		if(as.integer(nrow(test.data)/n.breeds)*n.breeds==nrow(test.data)){
			xx=as.factor(rep(letters[1:n.breeds], each=as.integer(nrow(test.data)/n.breeds)))
		}else{
			xx=as.factor(rep(letters[1:n.breeds], each=as.integer(nrow(test.data)/n.breeds)+1))
		}
		test.data$breed=sample(xx, size=nrow(test.data), replace=F)
		m.mat=model.matrix(object=~belief+sex+baited.fst+age, data=test.data)#xdata
		xx=as.vector(m.mat%*%fe2sim)+rnorm(n=n.breeds, sd=re2sim["(Intercept)"])[as.numeric(test.data$breed)]+
			rnorm(n=n.breeds, sd=re2sim["belieftrue"])[as.numeric(test.data$breed)]*m.mat[, "belieftrue"]
		test.data$choice=rbinom(n=nrow(test.data), size=1, prob=plogis(xx))
		#sum(test.data$choice)
		test.data$z.age=as.vector(scale(test.data$age))
		test.data$belief.code=as.numeric(test.data$belief==levels(test.data$belief)[2])
		test.data$sex.code=as.numeric(test.data$sex==levels(test.data$sex)[2])
		test.data$baited.fst.code=as.numeric(test.data$baited.fst==levels(test.data$baited.fst)[2])
		#test.data$belief.code=test.data$belief.code-mean(test.data$belief.code)
		#test.data$sex.code=test.data$sex.code-mean(test.data$sex.code)
		#test.data$baited.fst.code=test.data$baited.fst.code-mean(test.data$baited.fst.code)
		full=try(glmer(choice~belief+sex+baited.fst+z.age+(1+belief.code+sex.code+baited.fst.code+z.age||breed),
			control=contr, family=binomial, data=test.data), silent=T)
		if(class(full)[[1]]!="try-error"){
			tests=drop1p(model.res=full, para=F, contr=contr)$drop1.res
			return(list(
				fe=summary(full)$coefficients,
				re=as.data.frame(summary(full)$varcor),
				tests=tests,
				rv=table(test.data$choice),
				warns=gsub(x=paste(lmer.warns(full), collapse=""), pattern="boundary (singular) fit: see ?isSingular", replacement="", fixed=T)
			))
		}else{
			return(NULL)
		}
	})
	save(file=paste(c("/home/mundryr/transfer/lucrezia_bwsimwbreed_", n.per.group[n], ".RData"), collapse=""), list="i.res")
	print(paste(c(n.per.group[n], max(n.per.group)), collapse=" out of "))
}


flist=list.files(path="/home/roger/roger/2020/vienna/done/lucrezia_sim_bw_gr/wsps", full.names=T)
all.res=gsub(x=flist, pattern="/home/roger/roger/2020/vienna/done/lucrezia_sim_bw_gr/wsps/lucrezia_bwsimwbreed_", replacement="", fixed=T)
all.res=as.numeric(gsub(x=all.res, pattern=".RData", replacement="", fixed=T))
flist=flist[order(all.res)]
all.res=data.frame(N=all.res[order(all.res)]*8)
all.res$n.emp=NA
all.res$n.sims=NA
all.res$prop.success=NA
all.res$power.belief=NA
all.res$power.sex=NA
all.res$power.baited.fst=NA
all.res$power.z.age=NA
all.warns=matrix(NA, ncol=5, nrow=length(flist))
colnames(all.warns)=c("full", "belief", "sex", "baited.fst", "age")
all.fe.ests=vector("list", nrow(all.res))
all.re.ests=vector("list", nrow(all.res))
for(i in 1:length(flist)){
	load(flist[i])
	N.es=sum(!unlist(lapply(i.res, is.null)))
	all.res$n.emp[i]=sum(unlist(lapply(i.res, function(x){sum(x$rv)})))/N.es
	all.res$n.sims[i]=N.es
	all.res$prop.success[i]=mean(unlist(lapply(i.res, function(x){x$rv[2]/sum(x$rv)})))
	xx=lapply(i.res, function(x){x$tests[c("belief", "sex", "baited.fst", "z.age"), "Pr..Chisq."]})
	xx=matrix(unlist(xx), nrow=N.es, byrow=T)
	xx=apply(xx<=0.05, 2, sum, na.rm=T)
	all.res$power.belief[i]=xx[1]/n.sims
	all.res$power.sex[i]=xx[2]/n.sims
	all.res$power.baited.fst[i]=xx[3]/n.sims
	all.res$power.z.age[i]=xx[4]/n.sims
	all.warns[i, "full"]=mean(nchar(unlist(lapply(i.res, "[[", "warns")))>0)
	xx=unlist(lapply(i.res, is.null))
	xx=lapply(i.res[!xx], function(x){
		xx=apply(x$tests[c("belief", "sex", "baited.fst", "z.age"), c("n.opt.warnings", "n.fun.warnings")], 1, function(x){
			any(nchar(gsub(x=x, pattern="boundary (singular) fit: see ?isSingular", replacement="", fixed=T))>0)
		})
	})
	xx=matrix(unlist(xx), nrow=length(xx), byrow=T)
	all.warns[i, -1]=apply(xx, 2, mean, na.rm=T)
	xx=unlist(lapply(i.res, is.null))
	xx=lapply(i.res[!xx], function(x){x$fe[c("(Intercept)", "belieftrue", "sexM", "baited.fstgrey", "z.age"), "Estimate"]})
	xx=matrix(unlist(xx), nrow=length(xx), byrow=T)
	colnames(xx)=c("(Intercept)", "belieftrue", "sexM", "baited.fstgrey", "z.age")
	all.fe.ests[[i]]=xx
	
	xx=unlist(lapply(i.res, is.null))
	xx=lapply(i.res[!xx], function(x){xx=x$re[, "sdcor"]; names(xx)=x$re[, "var1"]; return(xx[c("(Intercept)", "belief.code", "sex.code", "baited.fst.code", "z.age")])})
	xx=matrix(unlist(xx), nrow=length(xx), byrow=T)
	colnames(xx)=c("(Intercept)", "belief.code", "sex.code", "baited.fst.code", "z.age")
	all.re.ests[[i]]=xx	
}

##plot fitting problems:
dev.off()
X11(height=5)
par(mar=c(3, 3, 0.2, 0.2), mgp=c(1.7, 0.3, 0), tcl=-0.15)
plot(1, 1, type="n", xlab="nr. dogs", ylab="probability", xlim=range(all.res$N), ylim=range(all.warns))
pchs=c(1, 2, 4, 5, 22)
for(i in 1:ncol(all.warns)){
	points(x=all.res$N, y=all.warns[, i], pch=pchs[i])
}
points(x=all.res$N, y=1-all.res$n.sims/n.sims, pch=2)
legend("topright", legend=c("full: error", paste(c("full", "belief", "sex", "baited.fst", "age"), "warrning", sep=": ")),
	pch=c(2, pchs))
savePlot(file="/home/roger/roger/2020/vienna/done/lucrezia_sim_bw_gr/sim_feb_19_plots/problems.png", type="png")

##plot power:
dev.off()
X11(height=5)
par(mar=c(3, 3, 0.2, 0.2), mgp=c(1.7, 0.3, 0), tcl=-0.15)
plot(all.res$N, all.res$power.belief, xlab="nr. dogs", ylab="probability", pch=1, xlim=range(all.res$N), ylim=c(0, 1))
others=c("power.sex", "power.baited.fst", "power.z.age")
pchs=3:5
for(i in 1:length(others)){
	points(all.res$N, all.res[, others[i]], pch=pchs[i])
}
legend("right", legend=c("belief", "sex", "baited.fst", "age"), pch=c(1, pchs))
abline(h=c(0.05, 0.8), lty=3)
savePlot(file="/home/roger/roger/2020/vienna/done/lucrezia_sim_bw_gr/sim_feb_19_plots/power.png", type="png")

all.fe.range=c()
all.re.range=c()
for(i in 1:nrow(all.res)){
	all.fe.range=rbind(all.fe.range, apply(all.fe.ests[[i]], 2, range))
	all.re.range=rbind(all.re.range, apply(all.re.ests[[i]], 2, range))
}
all.fe.range=apply(all.fe.range, 2, range)
all.re.range=apply(all.re.range, 2, range)

transf.est<-function(x){
	y=x
	tr.y=y
	tr.y[abs(tr.y)>1]=(1+log(abs(tr.y[abs(tr.y)>1])))*sign(tr.y[abs(tr.y)>1])
	yy=abs(range(y))
	if(yy[1]>0){
		y.lab=c(-1*rev(2^(0:as.integer(log2(yy[1])))), 0, 2^(0:as.integer(log2(yy[2]))))
	}else{
		y.lab=c(0, 2^(0:as.integer(log2(yy[2]))))
	}
	#browser()
	y.at=y.lab
	y.at[abs(y.at)>1]=(1+log(abs(y.at[abs(y.at)>1])))*sign(y.at[abs(y.at)>1])
	return(list(tr.x=tr.y, x.at=y.at, x.lab=y.lab))
}

fe2sim.r=fe2sim
names(fe2sim.r)[names(fe2sim.r)=="age"]="z.age"
dev.off()
X11(width=12, height=7)
par(mar=c(3, 3, 0.2, 0.2), mgp=c(1.7, 0.3, 0), tcl=-0.15, las=1)#, mfrow=c(2, 3))
for(i in 1:ncol(all.fe.range)){
	xx=transf.est(all.fe.range[, i])
	plot(1, 1, type="n", xlim=range(all.res$n.emp), ylim=range(xx$tr.x), xlab="nr. dogs", ylab="estimate",
		yaxt="n")
	axis(side=2, at=xx$x.at, labels=xx$x.lab, las=1)
	mtext(side=3, line=-1.2, text=colnames(all.fe.range)[i])
	for(j in 1:length(all.fe.ests)){
		points(x=all.res$n.emp[j]+runif(n=nrow(all.fe.ests[[j]]), min=-1.5, max=1.5), y=transf.est(all.fe.ests[[j]][, i])$tr.x,
			cex=0.2, pch=19, col=grey(level=0.25, alpha=0.25))
		to.add=quantile(transf.est(all.fe.ests[[j]][, i])$tr.x, probs=c(0.025, 0.25, 0.5, 0.75, 0.975))
		rect(xleft=all.res$n.emp[j]-1.75, xright=all.res$n.emp[j]+1.5, ybottom=to.add["25%"], ytop=to.add["75%"], border="red", lwd=2)
		segments(x0=all.res$n.emp[j], x1=all.res$n.emp[j], y0=to.add[c("25%", "75%")], y1=to.add[c("2.5%", "97.5%")], col="red", lwd=2)
		segments(x0=all.res$n.emp[j]-1.75, x1=all.res$n.emp[j]+1.5, y0=to.add["50%"], y1=to.add["50%"], lend=1, col="red", lwd=4)
	}
	abline(h=fe2sim.r[colnames(all.fe.range)[i]], col="blue")
	xx=colnames(all.fe.range)[i]
	xx=gsub(xx, pattern="(", replacement="", fixed=T)
	xx=gsub(xx, pattern=")", replacement="", fixed=T)
	xx=gsub(xx, pattern=".", replacement="", fixed=T)
	savePlot(file=paste(c("/home/roger/roger/2020/vienna/done/lucrezia_sim_bw_gr/sim_feb_19_plots/FE_", xx, ".png"), collapse=""), type="png")
}
	
re2sim.r=re2sim
names(re2sim.r)=colnames(all.re.range)
dev.off()
X11(width=12, height=7)
par(mar=c(3, 3, 0.2, 0.2), mgp=c(1.7, 0.3, 0), tcl=-0.15, las=1)#, mfrow=c(2, 3))
for(i in 1:ncol(all.fe.range)){
	xx=transf.est(all.re.range[, i])
	plot(1, 1, type="n", xlim=range(all.res$n.emp), ylim=range(xx$tr.x), xlab="nr. dogs", ylab="estimate",
		yaxt="n")
	axis(side=2, at=xx$x.at, labels=xx$x.lab, las=1)
	mtext(side=3, line=-1.2, text=colnames(all.re.range)[i])
	for(j in 1:length(all.fe.ests)){
		points(x=all.res$n.emp[j]+runif(n=nrow(all.re.ests[[j]]), min=-1.5, max=1.5), y=transf.est(all.re.ests[[j]][, i])$tr.x,
			cex=0.2, pch=19, col=grey(level=0.25, alpha=0.25))
		to.add=quantile(transf.est(all.re.ests[[j]][, i])$tr.x, probs=c(0.025, 0.25, 0.5, 0.75, 0.975))
		rect(xleft=all.res$n.emp[j]-1.75, xright=all.res$n.emp[j]+1.5, ybottom=to.add["25%"], ytop=to.add["75%"], border="red", lwd=2)
		segments(x0=all.res$n.emp[j], x1=all.res$n.emp[j], y0=to.add[c("25%", "75%")], y1=to.add[c("2.5%", "97.5%")], col="red", lwd=2)
		segments(x0=all.res$n.emp[j]-1.75, x1=all.res$n.emp[j]+1.5, y0=to.add["50%"], y1=to.add["50%"], lend=1, col="red", lwd=4)
	}
	abline(h=re2sim.r[colnames(all.re.range)[i]], col="blue")
	xx=colnames(all.re.range)[i]
	xx=gsub(xx, pattern="(", replacement="", fixed=T)
	xx=gsub(xx, pattern=")", replacement="", fixed=T)
	xx=gsub(xx, pattern=".", replacement="", fixed=T)
	savePlot(file=paste(c("/home/roger/roger/2020/vienna/done/lucrezia_sim_bw_gr/sim_feb_19_plots/RE_", xx, ".png"), collapse=""), type="png")
}


save.image("/home/roger/roger/2020/vienna/done/lucrezia_sim_bw_gr/power_w_one_tr_per_dog_new_feb.RData")
load("/home/roger/roger/2020/vienna/done/lucrezia_sim_bw_gr/power_w_one_tr_per_dog_new_feb.RData")

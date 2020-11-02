library(lme4)
library(gdata)
source("/home/mundryr/r_functions/helpers.r")
belief=as.factor(c("false", "true"))
t.order=as.factor(c("f.fst", "t.fst"))
fst.hide.loc=as.factor(c("A", "B"))

base.data=data.frame(expand.grid(belief=belief, t.order=t.order))
##fixed effects:
m.mat=model.matrix(object=~t.order*belief, data=base.data)
coefs=rep(NA, ncol(m.mat))
names(coefs)=colnames(m.mat)
coefs
##we have two scenarios re the order effect...
coefs["(Intercept)"]=qlogis(0.65)
coefs["belieftrue"]=qlogis(0.45)-coefs["(Intercept)"]
coefs["t.ordert.fst"]=qlogis(0.75)-coefs["(Intercept)"]
coefs["t.ordert.fst:belieftrue"]=(qlogis(0.35)-qlogis(0.75))-(qlogis(0.45)-coefs["(Intercept)"])
cbind(base.data, plogis(m.mat[, names(coefs)]%*%coefs))
#coefs=list(s1=coefs, s2=coefs)
#coefs[[2]]["t.ordert.fsf"]=qlogis(0.65)-qlogis(0.25)
#coefs[[2]]["t.ordert.fsf:belieftrue"]=(qlogis(0.45)-qlogis(0.45))-(qlogis(0.65)-qlogis(0.25))
#cbind(base.data, plogis(m.mat[, names(coefs[[2]])]%*%coefs[[2]]))
##random effects:
sd.ind.icpt=0.5
r.sl.order=0.5
r.sl.tf=0.5
##sample sizes:
N=(3:25)
n.subj.per.cell=8
n.tests.per.sub=c(2, 4)

contr=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=10000))
n.sims=1000
library(parallel)
cl <- makeCluster(getOption("cl.cores", detectCores()))
cl=cl[1:(length(cl)-1)]
parLapply(X=1:length(cl), cl=cl, fun=function(i){library(lme4); library(gdata)})
for(ntps in 1:length(n.tests.per.sub)){
  for(n in 1:length(N)){
    xdata=data.frame(expand.grid(t.order=t.order, fst.hide.loc=fst.hide.loc, sex=as.factor(c("F", "M"))))
    xdata=xdata[rep(1:nrow(xdata), times=N[n]), ]
    xdata=xdata[rep(1:nrow(xdata), each=2), ]
    n.subj=N[n]*n.subj.per.cell
    xdata$belief=rep(belief, n.subj)
    xdata$subj=as.factor(rep(paste("s", 1:n.subj, sep="."), each=2))
    xdata$t.order.c=as.numeric(xdata$t.order==levels(xdata$t.order)[2])
    xdata$t.order.c=xdata$t.order.c-mean(xdata$t.order.c)
    xdata$belief.c=as.numeric(xdata$belief==levels(xdata$belief)[2])
    xdata$belief.c=xdata$belief.c-mean(xdata$belief.c)
    if(ntps==2){
      xdata=xdata[rep(1:nrow(xdata), each=2), ]
    }
    m.mat=model.matrix(object=~t.order*belief, data=xdata)
    #all.res=vector("list", n.sims)
    #for(i in 1:n.sims){
    clusterExport(cl=cl, varlist=c("ntps", "xdata", "coefs", "sd.ind.icpt", "r.sl.order", "r.sl.tf", "contr", "lmer.warns", "m.mat", "n.subj"))
    all.res=parLapply(X=1:n.sims, cl=cl, fun=function(i){
      if(ntps==2){
        tr.nr=lapply(unique(as.character(xdata$subj)), function(s){
          res=rep(NA, 4)
          if(unique(as.character(xdata$t.order)[xdata$subj==s])=="t.fst"){
            res[resample(which((as.character(xdata$belief)[xdata$subj==s])=="true"), 1)]=1
          }else{
            res[resample(which((as.character(xdata$belief)[xdata$subj==s])=="false"), 1)]=1
          }
          res[is.na(res)]=sample(2:4, 3, replace=F)
          return(res)
        })
        xdata$tr.nr=unlist(tr.nr)
        #tapply(xdata$tr.nr, list(xdata$subj, xdata$t.order, xdata$belief), min)
        xdata$z.tr.nr=as.vector(scale(xdata$tr.nr))
      }
      xdata$choice=m.mat[, names(coefs)]%*%coefs+##fixed effects
        rnorm(n=n.subj, sd=sd.ind.icpt)[as.numeric(xdata$subj)]+##random intercept
        (rnorm(n=n.subj, sd=r.sl.order)*m.mat[, "t.ordert.fst"])[as.numeric(xdata$subj)]+##random slope of order
        (rnorm(n=n.subj, sd=r.sl.tf)*m.mat[, "belieftrue"])[as.numeric(xdata$subj)]##random slope of true/false belief
      xdata$choice=rbinom(n=nrow(xdata), size=1, prob=plogis(xdata$choice))
      if(ntps==1){
        full=try(glmer(choice~t.order*belief+fst.hide.loc+sex+(1|subj), data=xdata, family=binomial, control=contr), silent=T)
        null=try(glmer(choice~t.order+fst.hide.loc+sex+(1|subj), data=xdata, family=binomial, control=contr), silent=T)
      }else{
        full=try(glmer(choice~t.order*belief+fst.hide.loc+sex+z.tr.nr+(1+t.order.c+belief.c||subj), data=xdata, family=binomial, control=contr), silent=T)
        null=try(glmer(choice~t.order+fst.hide.loc+sex+z.tr.nr+(1+t.order.c+belief.c||subj), data=xdata, family=binomial, control=contr), silent=T)
      }
      xx=try(summary(full)$coefficients, silent=T)
      if(class(full)[[1]]!="try-error" & class(null)[[1]]!="try-error" & class(xx)!="try-error"){
        test.fn=as.data.frame(anova(null, full, test="Chisq"))[2, c("Chisq", "Chi Df", "Pr(>Chisq)")]
        tests=as.data.frame(drop1(full, test="Chisq"))[-1, c("LRT", "Df", "Pr(Chi)")]
        full.warns=paste(lmer.warns(full), collapse="")
        null.warns=paste(lmer.warns(null), collapse="")
        full.warns[full.warns==" "]=""
        null.warns[null.warns==" "]=""
        #all.res[[i]]=list(
        return(list(
          fe.coef=summary(full)$coefficients, re.coef=summary(full)$varcor,
          tests.fn=test.fn, tests=tests, 
          data.fe=aggregate(xdata$choice, xdata[, c("t.order", "belief")], sum),
          data.re=table(apply(table(xdata$subj, xdata$choice), 1, max)==2),
          full.warns=full.warns, null.warns=null.warns, rv.tab=table(xdata$choice),
          mean.per.sex=tapply(xdata$choice, xdata$sex, mean),
          mean.per.fst.hide.loc=tapply(xdata$choice, xdata$fst.hide.loc, mean),
          mean.per.ob=tapply(xdata$choice, list(xdata$t.order, xdata$belief), mean),
          N=nrow(xdata), rv=table(xdata$choice)
        ))
      }else{
        return(NULL)
      }
      #print(paste(c(n.tests.per.sub[ntps], N[n], i), collapse=" "))
    })
    all.res=all.res[!unlist(lapply(all.res, is.null))]
    save(file=paste(c("/home/mundryr/roger/lucretia_power_sim/c_", n.tests.per.sub[ntps], "__N_", N[n], ".RData"), collapse=""), list="all.res")
    print(paste(c(n.tests.per.sub[ntps], N[n]), collapse=" "))
  }
}

parLapply(cl=cl, 1:length(cl), fun=function(x){
  return(rm(list=ls()))
})
stopCluster(cl=cl)

sim.coefs=list(rep(0, ncol(all.fe[[1]])-1), rep(0, ncol(all.fe[[2]])-1))
names(sim.coefs[[1]])=colnames(all.fe[[1]])[-1]
names(sim.coefs[[2]])=colnames(all.fe[[2]])[-1]
sim.coefs[[1]][c(1:3, 6)]=coefs
sim.coefs[[2]][c(1:3, 7)]=coefs
names(sim.coefs[[1]])[1]="(Intercept)"
names(sim.coefs[[2]])[1]="(Intercept)"
names(sim.coefs[[1]])[names(sim.coefs[[1]])=="t.ordert.fst.belieftrue"]="t.ordert.fst:belieftrue"
names(sim.coefs[[2]])[names(sim.coefs[[2]])=="t.ordert.fst.belieftrue"]="t.ordert.fst:belieftrue"

flist=list.files(path="/home/mundryr/roger/lucretia_power_sim/", full.names=T)
eval=data.frame(expand.grid(N=N, trials=c(2, 4)))
eval$n.conv=NA
eval$n.full.probl=NA
eval$n.null.probl=NA
eval$n.fn.probl=NA
eval$power=NA
all.fe=vector("list", 2)
all.re=vector("list", 2)
all.fe.wrong=vector("list", 2)
for(i in 1:nrow(eval)){
  load(paste(c("/home/mundryr/roger/lucretia_power_sim/c_", eval$trials[i], "__N_", eval$N[i], ".RData"), collapse=""))
  eval$n.conv[i]=length(all.res)
  xx=unlist(lapply(all.res, "[[", "full.warns"))
  xx=gsub(xx, pattern="boundary (singular) fit: see ?isSingular", replacement="", fixed=T)
  eval$n.full.probl[i]=sum(xx!="")
  xxn=unlist(lapply(all.res, "[[", "full.warns"))
  xxn=gsub(xxn, pattern="boundary (singular) fit: see ?isSingular", replacement="", fixed=T)
  eval$n.null.probl[i]=sum(xxn!="")
  eval$n.fn.probl[i]=sum(xx!="" | xxn!="")
  xx=lapply(all.res, "[[", "tests.fn")
  xx=matrix(unlist(xx), ncol=3, byrow=T)
  eval$power[i]=sum(xx[, 3]<=0.05)
  
  i.fe=lapply(all.res, "[[", "fe.coef")
  xnames=rownames(i.fe[[1]])
  i.fe=lapply(i.fe, function(x){x[xnames, 1]})
  i.fe=matrix(unlist(i.fe), nrow=length(i.fe), byrow=T)
  colnames(i.fe)=xnames
  all.fe[[eval$trials[i]/2]]=rbind(all.fe[[eval$trials[i]/2]], data.frame(N=eval$N[i], i.fe))

  ##sign. in wrong direction:
  w.i.fe=lapply(all.res, "[[", "fe.coef")
  xnames=rownames(w.i.fe[[1]])
  w.i.fe=lapply(w.i.fe, function(x){x[xnames, "Pr(>|z|)"]})
  w.i.fe=matrix(unlist(w.i.fe), nrow=length(w.i.fe), byrow=T)
  colnames(w.i.fe)=xnames
  all.fe.wrong[[eval$trials[i]/2]]=rbind(all.fe.wrong[[eval$trials[i]/2]], 
    data.frame(N=eval$N[i], t(apply(t(sign(sim.coefs[[eval$trials[i]/2]][xnames])!=t(sign(i.fe))) & w.i.fe<=0.05, 2, sum))))

  i.fe=lapply(all.res, "[[", "re.coef")
  i.fe=lapply(i.fe, function(x){as.data.frame(x)})
  xnames=i.fe[[1]][, "var1"]
  i.fe=lapply(i.fe, function(x){x[, 5]})
  i.fe=matrix(unlist(i.fe), nrow=length(i.fe), byrow=T)
  colnames(i.fe)=xnames
  all.re[[eval$trials[i]/2]]=rbind(all.re[[eval$trials[i]/2]], data.frame(N=eval$N[i], i.fe))
}  

x.at=sort(unique(eval$N))*n.subj.per.cell
x.at=x.at[(1:2)==1]
dev.off()
X11(width=8, height=5)
par(mar=c(3, 3, 1.5, 0.2), mgp=c(1.7, 0.3, 0), tcl=-0.15, las=1)
plot(eval$N*n.subj.per.cell, eval$power/1000, xlab="number individuals", ylab="power", ylim=c(0, 1), pch=c(1, 4)[eval$trials/2], xaxt="n")
axis(side=1, at=x.at)
legend("bottomright", legend=c("2 trials", "4 trials"), pch=c(1, 4), bty="n")
abline(h=0.8, lty=3)
savePlot(file="/home/mundryr/transfer/lucrezia_power_results/power.png", type="png")

dev.off()
X11(width=8, height=5)
par(mar=c(3, 3, 1.5, 0.2), mgp=c(1.7, 0.3, 0), tcl=-0.15, las=1)
plot(eval$N*n.subj.per.cell, eval$n.fn.probl/1000, xlab="number individuals", ylab="proportion convergence problems", 
  ylim=c(0, 1), pch=c(1, 4)[eval$trials/2], xaxt="n")
axis(side=1, at=x.at)
legend("topleft", legend=c("2 trials", "4 trials"), pch=c(1, 4), bty="n")
abline(h=0, lty=3)
savePlot(file="/home/mundryr/transfer/lucrezia_power_results/conv_probl.png", type="png")

all.ip=list(c(), c())
for(i in 1:nrow(eval)){
  load(paste(c("/home/mundryr/roger/lucretia_power_sim/c_", eval$trials[i], "__N_", eval$N[i], ".RData"), collapse=""))
  xx=lapply(all.res, "[[", "tests")
  xnames=rownames(xx[[1]])
  xx=lapply(xx, function(x){x[xnames, "Pr(Chi)"]})
  xx=matrix(unlist(xx), nrow=length(xx), byrow=T)
  colnames(xx)=xnames
  all.ip[[eval$trials[i]/2]]=rbind(all.ip[[eval$trials[i]/2]], data.frame(n=eval$N[i]*n.subj.per.cell, t(apply(xx<=0.05, 2, mean))))
}

dev.off()
X11(width=8, height=5)
par(mar=c(3, 3, 1.5, 0.2), mgp=c(1.7, 0.3, 0), tcl=-0.15, las=1)
xnames=unique(c(names(all.ip[[1]]), names(all.ip[[2]])))[-1]
for(i in 1:length(xnames)){
  plot(all.ip[[2]]$n, all.ip[[2]][, xnames[i]], ylim=c(0, 1), xlab="nr. dogs", ylab="prob. significant", pch=4, xaxt="n")
  axis(side=1, at=x.at)
  if(xnames[i]%in%names(all.ip[[1]])){
    points(all.ip[[1]]$n, all.ip[[1]][, xnames[i]])
  }
  mtext(text=xnames[i], side=3, line=0.2)
  if(xnames[i]!="t.order.belief"){
    abline(h=0.05, lty=3)
  }else{
    abline(h=0.8, lty=3)
  }
  if(xnames[i]!="z.tr.nr"){
    legend("topleft", legend=c("2 trials", "4 trials"), pch=c(1, 4), bty="n")
  }
  xx=gsub(xnames[i], pattern=".", replacement="", fixed=T)
  savePlot(file=paste(c("/home/mundryr/roger/lucretia_power_sim_plots/mean_sig_", xnames[i], ".png"), collapse=""), type="png")
}



sim.coefs=list(rep(0, ncol(all.fe[[1]])-1), rep(0, ncol(all.fe[[2]])-1))
names(sim.coefs[[1]])=colnames(all.fe[[1]])[-1]
names(sim.coefs[[2]])=colnames(all.fe[[2]])[-1]
sim.coefs[[1]][c(1:3, 6)]=coefs
sim.coefs[[2]][c(1:3, 7)]=coefs

hbw=2

for(ss in 1:2){
  sel.fe=all.fe[[ss]]
  sel.fe.w=all.fe.wrong[[ss]]
  for(i in 2:ncol(sel.fe)){
    xx=runif(n=nrow(sel.fe), min=-hbw, max=hbw)
    y=sel.fe[, i]
    tr.y=y
    tr.y[abs(tr.y)>1]=(1+log(abs(tr.y[abs(tr.y)>1])))*sign(tr.y[abs(tr.y)>1])
    yy=abs(range(y))
    y.lab=c(-1*rev(2^(0:yy[1])), 0, 2^(0:yy[2]))
    y.at=y.lab
    y.at[abs(y.at)>1]=(1+log(abs(y.at[abs(y.at)>1])))*sign(y.at[abs(y.at)>1])
    plot(sel.fe$N*n.subj.per.cell+xx, y=tr.y, pch=19, col=grey(alpha=0.3, level=0.1), xlab="nr. dogs", ylab="estimate", yaxt="n", cex=0.3, xaxt="n")
    axis(side=1, at=x.at)
    axis(side=2, at=y.at, labels=y.lab)
    mtext(text=colnames(sel.fe)[i], side=3, line=0.2)
    xx=tapply(tr.y, sel.fe$N*n.subj.per.cell, quantile, probs=c(0.25, 0.5, 0.75))
    xx=matrix(unlist(xx), nrow=length(xx), byrow=T)
    rect(xleft=sort(unique(sel.fe$N))*n.subj.per.cell-hbw, xright=sort(unique(sel.fe$N))*n.subj.per.cell+hbw, ybottom=xx[, 1], ytop=xx[, 3], border="red", lwd=2)
    segments(x0=sort(unique(sel.fe$N))*n.subj.per.cell-hbw, x1=sort(unique(sel.fe$N))*n.subj.per.cell+hbw, y0=xx[, 2], y1=xx[, 2], col="red", lwd=4, lend=1)
    xx=gsub(colnames(sel.fe)[i], pattern=".", replacement="", fixed=T)
    abline(h=sim.coefs[[ss]][i-1], col="blue", lty=2, lwd=2)
    abline(h=0, col="black", lwd=1, lty=3)
    mtext(text=colnames(sel.fe)[i], side=3, line=0.2)
    if(colnames(sel.fe)[i]%in%c("X.Intercept.", "t.ordert.fst", "belieftrue", "t.ordert.fst.belieftrue")){
      mtext(text=sel.fe.w[, colnames(sel.fe)[i]], at=sel.fe.w$N*n.subj.per.cell, line=-1, side=1, cex=0.9)
    }
    xx=gsub(colnames(sel.fe)[i], pattern=".", replacement="", fixed=T)
    savePlot(file=paste(c("/home/mundryr/roger/lucretia_power_sim_plots/fe_", paste(c(xx, ss*2), collapse="_"), ".png"), collapse=""), type="png")
  }  
}

for(ss in 1:2){
  sel.fe=all.re[[ss]]
  for(i in 2:ncol(sel.fe)){
    y=sel.fe[, i]
    tr.y=y
    if(ss==1 & i==2){
      tr.y[abs(tr.y)>1]=(1+log(abs(tr.y[abs(tr.y)>1])))*sign(tr.y[abs(tr.y)>1])
      yy=abs(range(y))
      y.lab=c(0, 2^(0:as.integer(log2(max(y)))))
      y.at=y.lab
      y.at[abs(y.at)>1]=(1+log(abs(y.at[abs(y.at)>1])))*sign(y.at[abs(y.at)>1])
    }else{
      y.lab=pretty(tr.y)
      y.at=y.lab
    }
    xx=runif(n=nrow(sel.fe), min=-2, max=2)
    plot(sel.fe$N*n.subj.per.cell+xx, tr.y, pch=19, cex=0.3, col=grey(alpha=0.3, level=0.1), xlab="nr. dogs", ylab="estimate", yaxt="n", xaxt="n")
    axis(side=1, at=x.at)
    axis(side=2, at=y.at, labels=y.lab)
    mtext(text=colnames(sel.fe)[i], side=3, line=0.2)
    xx=tapply(tr.y, sel.fe$N*n.subj.per.cell, quantile, probs=c(0.25, 0.5, 0.75))
    xx=matrix(unlist(xx), nrow=length(xx), byrow=T)
    rect(xleft=sort(unique(sel.fe$N))*n.subj.per.cell-hbw, xright=sort(unique(sel.fe$N))*n.subj.per.cell+hbw, ybottom=xx[, 1], ytop=xx[, 3], border="red", lwd=2)
    segments(x0=sort(unique(sel.fe$N))*n.subj.per.cell-hbw, x1=sort(unique(sel.fe$N))*n.subj.per.cell+hbw, y0=xx[, 2], y1=xx[, 2], col="red", lwd=4, lend=1)
    abline(h=0.5, col="blue", lty=2, lwd=2)
    xx=gsub(colnames(sel.fe)[i], pattern=".", replacement="", fixed=T)
    savePlot(file=paste(c("/home/mundryr/roger/lucretia_power_sim_plots/re_", paste(c(xx, ss*2), collapse="_"), ".png"), collapse=""), type="png")
  }  
}
  

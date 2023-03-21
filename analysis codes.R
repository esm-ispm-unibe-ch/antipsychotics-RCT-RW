
# load the data -----------------------------------------------------------
rwd=readRDS("rwd.Rdata") ## registry data
rct6=readRDS("rct6.RData") ## registry data


#compare the registries --------------------------------

## compare vs no AP use 
versus.reference=rwd[rwd$treat2=="no antipsychotic",]
dim(versus.reference)
treats.swe=unique(c(rwd$treat1[rwd$stud=="swe"], rwd$treat2[rwd$stud=="swe"]))
treats.fin=unique(c(rwd$treat1[rwd$stud=="fin"], rwd$treat2[rwd$stud=="fin"]))
treats=unique(treats.swe, treats.fin)
not.common.treats=c(treats.swe[!(treats.swe %in% treats.fin)], treats.fin[!(treats.fin %in% treats.swe)])
versus.reference=versus.reference[!(versus.reference$treat1 %in% not.common.treats),]
versus.reference$comp=paste(versus.reference$t1,"-",versus.reference$t2)
comp.dat=versus.reference[,c("stud", "HR", "low", "upp", "comp")]
comp.dat2=reshape(comp.dat, idvar = "comp", timevar="stud", direction = "wide")

library(ggplot2)

xlim=c(0.2,1.2)
ylim=c(0.2,1.2)
drug.comp=unlist(strsplit(comp.dat2$comp, " "))
drug.comp=treats[as.numeric(drug.comp[!(drug.comp %in% c("-","15"))])]
comp.dat2$comp=drug.comp
colnames(comp.dat2)[1]="Drug"
comp.dat2[,2:7]=round(comp.dat2[,2:7],2)

dput(drug.comp)
drug.comp2=c("HAL.O", "HAL.L", "CLO.O", 
             "SER.O", "ZIP.O", "OLA.O", "OLA.L", 
             "QUE.O", "RIS.O", "RIS.L", "ARI.O", 
             "PAL.O", "ARI.L")

## create scatterplot for HR vs no AP use
comp.dat3=comp.dat2
comp.dat3$upp.swe[comp.dat3$upp.swe>xlim[2]]=xlim[2]
comp.dat3$upp.fin[comp.dat3$upp.fin>ylim[2]]=ylim[2]
ggplot(data = comp.dat3,aes(x = HR.swe,y = HR.fin)) + 
  geom_point(size=3, color="Blue" ) + 
  geom_errorbar(aes(ymin = low.fin,ymax = upp.fin)) + 
  geom_errorbarh(aes(xmin = low.swe,xmax = upp.swe))+
  geom_abline(intercept = 0, slope = 1, color="blue", 
              linetype="dashed", size=0.5)+xlim(xlim)+ylim(ylim)+
  geom_text(label=drug.comp2,hjust=1.4, vjust=1.4)+xlab("Sweden")+ylab("Finland")+
  ggtitle("Hazard ratios vs. no antipsychotic use")+
  theme(plot.title = element_text(size=20))


## create forestplot of ratio of hazard ratios vs no AP use
comp.dat2$logRHR=log(comp.dat2$HR.swe/comp.dat2$HR.fin)
comp.dat2$SElogRHR=with(comp.dat2, sqrt((log(upp.swe/low.swe)/3.92)^2+
                                          (log(upp.fin/low.fin)/3.92)^2  )  )

library(meta)
m1=metagen(study=comp.dat2$tr,TE=comp.dat2$logRHR, seTE=comp.dat2$SElogRHR, sm="HR", back = T)
forest(m1, studlab = comp.dat2$Drug, leftcols = c("studlab"), leftlabs = "Drug", 
       backtransf = T, rightlabs = c("RHR", "95% CI", "Weight common", "Weight fixed"), 
       addrows.below.overall = 2, smlab = "Ratio of HR")



## weighted Spearman correlation for HRs vs no AP use
library(wCorr)
comp.dat2$weight=with(comp.dat2, 
                      1/ (((log(upp.swe)-log(low.swe))/3.92)^2+
                            ((log(upp.fin)-log(low.fin))/3.92)^2))
weightedCorr(comp.dat2$HR.swe, comp.dat2$HR.fin,  method = "spearman", 
             weights=comp.dat2$weight)



## compare vs oral haloperidol
versus.reference=rwd[rwd$treat1=="Haloperidol oral",]
versus.reference=versus.reference[!(versus.reference$treat2 %in% not.common.treats),]
dim(versus.reference)
versus.reference$comp=paste(versus.reference$t1,"-",versus.reference$t2)
comp.dat=versus.reference[,c("stud", "HR", "low", "upp", "comp")]
comp.dat2=reshape(comp.dat, idvar = "comp", timevar="stud", direction = "wide")
comp.dat2[,2:7]=1/comp.dat2[,2:7]
colnames(comp.dat2)=c("comp", "HR.swe", "upp.swe", "low.swe", 
                      "HR.fin",        "upp.fin", "low.fin")
xlim=c(0.2,1.6)
ylim=c(0.2,1.6)
drug.comp=unlist(strsplit(comp.dat2$comp, " "))
drug.comp=treats[as.numeric(drug.comp[!(drug.comp %in% c("1","-"))])]
comp.dat2$comp=drug.comp
colnames(comp.dat2)[1]="Drug"
comp.dat2[,2:7]=round(comp.dat2[,2:7],2)
comp.dat2=comp.dat2[comp.dat2$Drug!="no antipsychotic",]
dput(comp.dat2$Drug)
drug.comp2=c("HAL.L", "CLO.O", 
             "SER.O", "ZIP.O", "OLA.O", "OLA.L", 
             "QUE.O", "RIS.O", "RIS.L", "ARI.O", 
             "PAL.L", "ARI.L")

comp.dat3=comp.dat2
comp.dat3$upp.swe[comp.dat3$upp.swe>xlim[2]]=xlim[2]
comp.dat3$upp.fin[comp.dat3$upp.fin>ylim[2]]=ylim[2]
ggplot(data = comp.dat3,aes(x = HR.swe,y = HR.fin)) + 
  geom_point(size=3, color="Blue") + 
  geom_errorbar(aes(ymin = low.fin,ymax = upp.fin), size=0.2) + 
  geom_errorbarh(aes(xmin = low.swe,xmax = upp.swe), size=0.2)+
  geom_abline(intercept = 0, slope = 1, color="black", 
              linetype="dashed", size=0.2)+xlim(xlim)+ylim(ylim)+
  geom_text(label=drug.comp2, size=3,hjust=1.4, vjust=1.4)


## forestplot of RHR vs oral haloperidol
comp.dat2$logRHR=log(comp.dat2$HR.swe/comp.dat2$HR.fin)
comp.dat2$SElogRHR=with(comp.dat2, sqrt((log(upp.swe/low.swe)/3.92)^2+
                                          (log(upp.fin/low.fin)/3.92)^2  )  )

m1=metagen(study=comp.dat2$tr,TE=comp.dat2$logRHR, seTE=comp.dat2$SElogRHR, sm="HR", back = T)
forest(m1, studlab = comp.dat2$Drug, leftcols = c("studlab"), leftlabs = "Drug", 
       backtransf = T, rightlabs = c("RHR", "95% CI", "Weight common", "Weight fixed"), addrows.below.overall = 2, smlab = "Ratio of HR")



#create forestplot of the HRs vs oral haloperidol by registry
drug.comp=unlist(strsplit(comp.dat$comp, " "))
comp.dat$drug=treats[as.numeric(drug.comp[!(drug.comp %in% c("1","-"))])]
comp.dat$HR=log(comp.dat$HR)
comp.dat$se=log(comp.dat$upp/comp.dat$low)/3.92
comp.dat=comp.dat[comp.dat$drug!="no antipsychotic",]
comp.dat$stud=gsub("swe", "Sweden", comp.dat$stud)
comp.dat$stud=gsub("fin", "Finland", comp.dat$stud)
one.stage=metagen(studlab = comp.dat$stud, TE=-comp.dat$HR, seTE=comp.dat$se, 
                  subgroup=comp.dat$drug, 
                  sm="OR")
forest(one.stage, overall = F, subgroup.hetstat = F, subgroup = F, leftcols=	c("studlab"), rightlabs = c("HR", "95% CI") 
       , smlab = "Hazard Ratio",  addrows.below.overall = 2, squaresize = 0.3, xlim=c(0.2, 2.0), overall.hetstat = F, 
       test.subgroup = F,col.square = "BLUE", weight.study = "same", leftlabs = c("Comparison vs. Haloperidol oral"))



### network meta-analysis of RWD only 
library(netmeta)
net.rwd=netmeta(studlab =stud, TE=TE, seTE=seTE,sm="HR", data=rwd, 
                treat1=treat1, treat2=treat2,tol.multiarm = 0.01)

summary(net.rwd)
decomp.design(net.rwd)

netsplit.rct=   netsplit(net.rwd, fixed=F)
print(netsplit.rct, show="both", digits=2)

netrank.rwd=netrank(net.rwd)$ranking.random
forest(net.rwd, sortvar=names(netrank.rwd[order(-netrank.rwd)]),
       reference.group = "no antipsychotic", xlim=c(0.2,2))


# RW comparisons vs. Haloperidol oral
TE.rwd.hal.oral<-net.rwd$TE.random
TE.rwd.hal.oral<-TE.rwd.hal.oral[colnames(TE.rwd.hal.oral)=="Haloperidol oral",]
seTE.rwd.hal.oral<-net.rwd$seTE.random
seTE.rwd.hal.oral<-seTE.rwd.hal.oral[colnames(seTE.rwd.hal.oral)=="Haloperidol oral",]
rwd.hal.oral<-data.frame(cbind(TE.rwd.hal.oral,seTE.rwd.hal.oral ))
rwd.hal.oral$HR=exp(- rwd.hal.oral$TE.rwd.hal.oral)
rwd.hal.oral$HR_low=exp(- rwd.hal.oral$TE.rwd.hal.oral-1.96*rwd.hal.oral$seTE.rwd.hal.oral)
rwd.hal.oral$HR_up=exp(- rwd.hal.oral$TE.rwd.hal.oral+1.96*rwd.hal.oral$seTE.rwd.hal.oral)
rwd.hal.oral[,3:5]=round(rwd.hal.oral[,3:5],2)
rwd.hal.oral$tr=row.names(rwd.hal.oral)
rwd.hal.oral[,3:5]


# NMA of RCTs -------------------------------------------------------------
net.rct=netmeta(studlab =studlab, TE=logHR_6m_selected, 
                seTE=SE.logHR_6m_selected,sm="HR", data=rct6, 
                treat1=treat1, treat2=treat2, reference.group = "placebo",
                tol.multiarm = 0.01)
netsplit.rct=   netsplit(net.rct, fixed=F)
print(netsplit.rct, show="both", digits=2)
netsplit.rct[complete.cases(netsplit.rct),]
decomp.design(net.rct)

netgraph(net.rct, col = "black", plastic = FALSE,
         points = TRUE, pch = 21, cex.points = 3, col.points = "black",
         bg.points = "gray",
         multiarm = FALSE, number = TRUE, pos.number.of.studies = 0.5)

n.summary2<-summary(net.rct)

netrank.rct=netrank(net.rct)$ranking.random
forest(net.rct, sortvar=names(netrank.rct[order(-netrank.rct)]), xlim=c(0.05,5))


# rct comparisons vs. Haloperidol oral
TE.rct.hal.oral<-net.rct$TE.random
TE.rct.hal.oral<-TE.rct.hal.oral[colnames(TE.rct.hal.oral)=="Haloperidol oral",]
seTE.rct.hal.oral<-net.rct$seTE.random
seTE.rct.hal.oral<-seTE.rct.hal.oral[colnames(seTE.rct.hal.oral)=="Haloperidol oral",]
rct.hal.oral<-data.frame(cbind(TE.rct.hal.oral,seTE.rct.hal.oral ))
rct.hal.oral$HR=exp(- rct.hal.oral$TE.rct.hal.oral)
rct.hal.oral$HR_low=exp(- rct.hal.oral$TE.rct.hal.oral-1.96*rct.hal.oral$seTE.rct.hal.oral)
rct.hal.oral$HR_up=exp(- rct.hal.oral$TE.rct.hal.oral+1.96*rct.hal.oral$seTE.rct.hal.oral)
rct.hal.oral[,3:5]=round(rct.hal.oral[,3:5],2)
rct.hal.oral$tr=row.names(rct.hal.oral)
rct.hal.oral[,3:5]=round(rct.hal.oral[,3:5],2)
rct.hal.oral$tr=row.names(rct.hal.oral)
rct.hal.oral[,3:5]


### compare between RCT and RWD ------- 

#compare vs Hal oral
rwd.hal.oral$source="rwd"
rct.hal.oral$source="rct"
not.common.treats=c(rwd.hal.oral$tr[!(rwd.hal.oral$tr %in% rct.hal.oral$tr)],
                    rct.hal.oral$tr[!(rct.hal.oral$tr %in% rwd.hal.oral$tr)])
versus.reference=rbind(rwd.hal.oral[,3:7], rct.hal.oral[,3:7])
versus.reference=versus.reference[!(versus.reference$tr %in% not.common.treats),]
rownames(versus.reference)<-NULL
comp.dat.hal=reshape(versus.reference, idvar = "tr", timevar="source", direction = "wide")


xlim=c(0.01,2.1)
ylim=c(0.01,2.1)
comp.dat3=comp.dat.hal
comp.dat3$HR_up.rct[comp.dat3$HR_up.rct>ylim[2]]=ylim[2]
comp.dat3$HR_up.rwd[comp.dat3$HR_up.rwd>xlim[2]]=xlim[2]
dput(comp.dat3$tr)

comp.dat3=comp.dat3[comp.dat3$tr!="Haloperidol oral",]
lab.RIS=c("ARI.L", "ARI.O", "HAL.L", 
          "OLA.L","OLA.O","PAL.L", "PAL.O","QUE.O" ,
          "RIS.L",  "RIS.O", "SER.O", 
          "ZIP.O")
ggplot(data = comp.dat3,aes(x = HR.rwd,y = HR.rct)) + 
  geom_point(size=3, color="Blue" ) + 
  geom_errorbar(aes(ymin = HR_low.rct,ymax = HR_up.rct)) + 
  geom_errorbarh(aes(xmin = HR_low.rwd,xmax = HR_up.rwd))+
  geom_abline(intercept = 0, slope = 1, color="blue", 
              linetype="dashed", size=0.5)+xlim(xlim)+ylim(ylim)+
  geom_text(label=lab.RIS,hjust=1, vjust=1)



## weighted Speamrman correlation
comp.dat3=comp.dat.hal
comp.dat3=comp.dat3[comp.dat3$tr!="Haloperidol oral",]
comp.dat3$weight=with(comp.dat3, 
                      1/ (
                        ((log(HR_up.rwd)-log(HR_low.rwd))/3.92)^2+
                          ((log(HR_up.rct)-log(HR_low.rct))/3.92)^2         ))
weightedCorr(comp.dat3$HR.rct, comp.dat3$HR.rwd,  method = "spearman", 
             weights=comp.dat3$weight)


#RHR 
comp.dat.hal=comp.dat.hal[comp.dat.hal$tr!="Haloperidol oral",]
comp.dat.hal$logROR=with(comp.dat.hal, log(HR.rwd/HR.rct))
comp.dat.hal$selogROR=with(comp.dat.hal, sqrt((log(HR_up.rwd/HR_low.rwd)/3.92)^2+
                                                (log(HR_up.rct/HR_low.rct)/3.92)^2  )  )
m1=metagen(study=comp.dat.hal$tr,TE=comp.dat.hal$logROR, seTE=comp.dat.hal$selogROR, sm="OR")
forest(m1, studlab = paste(comp.dat.hal$tr), leftcols = c("studlab"), 
       leftlabs = "Comparison vs Haloperidol oral", 
       backtransf = T, rightlabs = c("RHR", "95% CI", "Weight common", "Weight fixed"), 
       addrows.below.overall = 2, smlab = "Ratio of HR")

# compare 
rwd.hal.oral$source="RW"
rct.hal.oral$source="RCT"
colnames(rwd.hal.oral)=colnames(rct.hal.oral)
bind.hal=rbind(rwd.hal.oral, rct.hal.oral)
colnames(bind.hal)[6]="Drug"
bind.hal=bind.hal[!(bind.hal$Drug %in% c("Haloperidol oral", "Placebo", "no antipsychotic")),]
bind.hal=bind.hal[!(bind.hal$Drug %in% not.common.treats),]
one.stage=metagen(studlab = bind.hal$source, TE=-bind.hal$TE.rct.hal.oral, seTE=bind.hal$seTE.rct.hal.oral, subgroup=bind.hal$Drug, 
                  sm="OR")
forest(one.stage, overall = F, subgroup.hetstat = F, subgroup = F, leftcols=	c("studlab"), rightlabs = c("HR", "95% CI") 
       , smlab = "Hazard Ratio",  addrows.below.overall = 2, squaresize = 0.3, xlim=c(0.1, 2), overall.hetstat = F, 
       test.subgroup = F,col.square = "BLUE", weight.study = "same", leftlabs = c("Comparison vs. Haloperidol oral"))






#### joint NMA -------
rwd2=rwd[,c("stud","treat1","treat2","TE","seTE")]
colnames(rct6)=colnames(rwd2)
all.data=rbind(rwd2, rct6[,1:5])
treats=unique(c(all.data$treat1, all.data$treat2))

net.all=netmeta(studlab =stud, TE=TE, seTE=seTE,sm="HR", data=all.data, 
                treat1=treat1, treat2=treat2, reference.group = "Haloperidol oral",
                tol.multiarm = 0.01)

netsplit.all=   netsplit(net.all, fixed=F)
print(netsplit.all, show="both", digits=2)
decomp.design(net.all)
netrank.all=netrank(net.all)$ranking.random

forest(net.all, sortvar=names(netrank.all[order(-netrank.all)]), 
       reference.group = "Haloperidol oral", xlim=c(0.05,5))



#Supplementary Information for Erkenbrack and Thompson, 2019. All directories correspond to my machine (JRT) so you'll need to change to your directories when re-running analyses. 

library("ape")
library("phytools")
library("phyloch")
library("strap")
library("coda")
library("xlsx")

library("devtools")
install_github("rgriff23/btw")
library("btw")


tree<-read.beast("/Users/jeffreythompson/Research/Embryos_and_Ancestors/BEAST_Analysis/Analysis_2_Larger_Final/Output_Tree.tre")
tree$root.time<-tree$height[1]

num_taxa<-length(tree$tip.label)




pdf("/Users/jeffreythompson/Research/Embryos_and_Ancestors/ASR/BayesTraits/geoscaled_tree.pdf", width=10, height=7)

geoscalePhylo(tree=ladderize(tree, right=FALSE), 
              quat.rm=TRUE, units=c("Period", "Epoch"), boxes="Period",
              tick.scale=50, cex.tip=1.2, cex.age=1, cex.ts=1, width=2,
              label.offset=0, x.lim=c(-90, 516), lwd=1.5)

lastPP<-get("last_plot.phylo", envir=.PlotPhyloEnv)


tree$node.label<-tree$posterior
p <- character(length(tree$node.label))
p[tree$node.label >= 0.95] <- "Light Blue"
p[tree$node.label < 0.95 & tree$node.label >= 0.75] <- "gray"
p[tree$node.label < 0.75] <- "white"
nodelabels(pch=21, cex=1.5, bg=p)

HPD_Max<-lastPP$xx[13:23]-(tree$"height_95%_HPD_MAX"-tree$height)

HPD_Min<-lastPP$xx[13:23]+(tree$height-tree$"height_95%_HPD_MIN")

for (n in 1:11){
  
  lines(c(HPD_Max[n], HPD_Min[n]), c(lastPP$yy[num_taxa+n],lastPP$yy[num_taxa+n]), col=rgb(0,0,1,alpha=0.3),lwd=8)
}

legend(40, 12, legend=c("95% Credible Interval"), col=rgb(0,0,1,alpha=0.3), lty=1, cex=1, lwd=8, box.lwd=2, bg="white")

dev.off()



trees_SCM_Full<-read.nexus(file="/Users/jeffreythompson/Research/Embryos_and_Ancestors/BEAST_Analysis/Analysis_2_Larger_Final/18s_Tree.trees") #Import Set of Trees from Posterior of Bayesian Analysis
trees_SCM_Burned<-trees_SCM_Full[20001:length(trees_SCM_Full)]
trees_SCM<-sample(trees_SCM_Burned, 10000)
write.tree(trees_SCM, file="/Users/jeffreythompson/Research/Embryos_and_Ancestors/BEAST_Analysis/Analysis_2_Larger_Final/Pruned_Sample_Trees")

setwd("/Users/jeffreythompson/Research/Embryos_and_Ancestors/ASR/BayesTraits")
#Define variables and states for analyses and plots
Alx1<-c("a","a","a","a","a","a","a","a","a","b","a","a")
names(Alx1)<-trees_SCM$STATE_138272000$tip.label
Alx_names<-matrix(names(Alx1))
Alx_values<-matrix(c(0,0,0,0,0,0,0,0,0,1,0,0), ncol=1, nrow=12)
Alx_col<-character(length(trees_SCM$tip.label))
Alx_col[Alx_values=="0"]<-"blue"
Alx_col[Alx_values=="1"]<-"grey"
Tbr_values<-matrix(c(1,1,1,0,0,0,0,3,"-",4,3,0), ncol=1, nrow=12)
Tbr_col<-character(length(trees_SCM$tip.label))
Tbr_col[Tbr_values=="0"]<-"green"
Tbr_col[Tbr_values=="1"]<-"blue"
Tbr_col[Tbr_values=="2"]<-"cyan"
Tbr_col[Tbr_values=="3"]<-"yellow"
Tbr_col[Tbr_values=="4"]<-"red"
Tbr_col[Tbr_values=="-"]<-"white"
Alx_Data<-cbind(Alx_names, Alx_values)
Ets1_Data<-Alx_Data
Tbr_Data<-cbind(Alx_names, Tbr_values)

Ets1_col<-character(length(trees_SCM$tip.label))
Ets1_col[Alx_values=="0"]<-"green"
Ets1_col[Alx_values=="1"]<-"grey"

#Command Files for all Comparisons of one rate to two rate models
Command_One_Rate<-c("1", "2","Res q10 q01", "ScaleTrees","PriorAll uniform 0 2", "Stones 1000 100000")
Command_Two_Rates<-c("1", "2", "ScaleTrees","PriorAll uniform 0 2", "Stones 1000 100000")

#For All Bayes Factor Tests at Asterozoan MRCA Define Command Files
Command_Fixed_State_Zero<-c("1", "2","Res q10 q01", "ScaleTrees", "Stones 1000 100000",
                            "Iterations 10000000",
                            "Burnin 2000000",
                            "PriorAll uniform 0 1",
                            "AddTag Asterozoans Patiria Amphiura",
                            "Fossil Fossil_Node Asterozoans 0")

Command_Fixed_State_One<-c("1", "2","Res q10 q01", "ScaleTrees", "Stones 1000 100000",
                           "Iterations 10000000",
                           "Burnin 2000000",
                           "PriorAll uniform 0 1",
                           "AddTag Asterozoans Patiria Amphiura",
                           "Fossil Fossil_Node Asterozoans 1")

#Define Command File for Alx1, Ets1 and Tbr
Command<-c("1", "2", "ScaleTrees",
           "Iterations 10000000",
           "Burnin 2000000",
           "Resall q10",
           "PriorAll uniform 0 .2",
           "AddTag Echinozoans Strongylocentrotus Holothuria",
           "AddMRCA MRCAEchinozoans Echinozoans",
           "AddTag Echinoids Strongylocentrotus Prionocidaris",
           "AddMRCA MRCAEchinoids Echinoids",
           "AddTag PurpEchinocardium Strongylocentrotus Echinocardium",
           "AddMRCA MCRCAPurpEchinocardium PurpEchinocardium",
           "AddTag PurpLytechinus Strongylocentrotus Lytechinus",
           "AddMRCA MRCAPurpLytechinus PurpLytechinus",
           "AddTag PurpParacentrotus Paracentrotus Strongylocentrotus",
           "AddMRCA MRCAPurpParacentrotus PurpParacentrotus",
           "AddTag Irregular Echinodiscus Echinocardium",
           "AddMRCA MRCAIrregular Irregular",
           "AddTag Cidaroids Eucidaris Prionocidaris",
           "AddMRCA MRCACidaroids Cidaroids",
           "AddTag Holothurians Holothuria Apostichopus",
           "AddMRCA MRCAHolothurians Holothurians",
           "AddTag Asterozoans Patiria Amphiura",
           "AddMRCA MRCAAsterozoans Asterozoans",
           "AddTag Ophiuroids Amphiura Amphipholis",
           "AddMRCA MRCAOphiuroids Ophiuroids")


#For Alx1
Alx_Dataframe<-as.data.frame(Alx_Data)

#Use Bayes Factor to Compare One Rate to Two Rate Models
BT_One_rate<-bayestraits(Alx_Dataframe, trees_SCM, Command_One_Rate)

BT_Two_rates<-bayestraits( Alx_Dataframe, trees_SCM, Command_Two_Rates)


Log_BF_Rates_Alx1<-2*(BT_Two_rates$Stones$logMarLH-BT_One_rate$Stones$logMarLH)


#Reconstruct Ancestral States Alx1
BT_Run_Alx1<-bayestraits(Alx_Dataframe, trees_SCM, Command)
traceplot(mcmc(BT_Run_Alx1$Log$results$q01))
densplot(mcmc(BT_Run_Alx1$Log$results$q01))
autocorr.plot(mcmc(BT_Run_Alx1$Log$results$q01))

#Plot & Summarize Results for Alx1
Node_Values_Alx<-BT_Run_Alx1$Log[2]
Node_Values_Alx<-Node_Values_Alx$results[,6:ncol(Node_Values_Alx$results)]
Mean_Nodes_Alx<-colMeans(Node_Values_Alx)
Node_List_Alx<-matrix(nrow=length(Mean_Nodes_Alx), ncol=2)
for (a in seq(from=1, to=length(Mean_Nodes_Alx), by=2)){
  Node<-c(Mean_Nodes_Alx[a], Mean_Nodes_Alx[a+1])
  Node_List_Alx[a,]<-c(Mean_Nodes_Alx[a], Mean_Nodes_Alx[a+1])
  pie(Node)
}

Node_List_Alx<-na.omit(Node_List_Alx)
geoscalePhylo(tree=ladderize(tree, right=TRUE), 
              quat.rm=TRUE, units=c("Period"), boxes="Period",
              tick.scale="Period", cex.tip=1, cex.age=.5, cex.ts=.7, width=2,
              label.offset=0, x.lim=c(-90, 516), lwd=1.5)
nodelabels(pie=Node_List_Alx, piecol=c("blue", "grey"))
Numbers<-(1:11)
Numbers<-as.character(Numbers)
nodelabels(Numbers, adj=-2)
tiplabels(pch=21, col="black",  bg=Alx_col, cex=2.3)
write.xlsx(Node_List_Alx, "/Users/jeffreythompson/Research/Embryos_and_Ancestors/ASR/BayesTraits/Real_Final/Sensitivity/0_0.2/Alx1/Alx1_States.xls")

pdf("/Users/jeffreythompson/Research/Embryos_and_Ancestors/ASR/BayesTraits/Real_Final/Sensitivity/0_0.2/Alx1/Alx1_Bayes_Traits.pdf", width=10, height=7)

geoscalePhylo(tree=ladderize(tree, right=TRUE), 
              quat.rm=TRUE, units=c("Period"), boxes="Period",
              tick.scale="Period", cex.tip=1, cex.age=.5, cex.ts=.7, width=2,
              label.offset=0, x.lim=c(-90, 516), lwd=1.5)
nodelabels(pie=Node_List_Alx, piecol=c("blue", "grey"))
Numbers<-(1:11)
Numbers<-as.character(Numbers)
nodelabels(Numbers, adj=-2)
tiplabels(pch=21, col="black",  bg=Alx_col, cex=2.3)

dev.off()


#Use Bayes Factor to Test Hypotheses regarding Alx1 expression at asterozoan MRCA
BT_Asterozoan_Zero<-bayestraits(Alx_Dataframe, trees_SCM, Command_Fixed_State_Zero)

BT_Asterozoan_One<-bayestraits(Alx_Dataframe, trees_SCM, Command_Fixed_State_One)

Log_BF_Asterozoan_Alx1<-2*(BT_Asterozoan_Zero$Stones$logMarLH-BT_Asterozoan_One$Stones$logMarLH)



#For Ets1
Ets1_Dataframe<-as.data.frame(Ets1_Data)
#Use Bayes Factor to Compare One Rate to Two Rate Models
BT_One_rate_Ets<-bayestraits(Ets1_Dataframe, trees_SCM, Command_One_Rate)

BT_Two_rates_Ets<-bayestraits(Ets1_Dataframe, trees_SCM, Command_Two_Rates)

Log_BF_Rates_Ets<-2*(BT_Two_rates_Ets$Stones$logMarLH-BT_One_rate_Ets$Stones$logMarLH)


#Reconstruct Ancestral States Ets1
BT_Run_Ets1<-bayestraits(Ets1_Dataframe, trees_SCM, Command)
traceplot(mcmc(BT_Run_Ets1$Log$results$q01))
densplot(mcmc(BT_Run_Ets1$Log$results$q01))
autocorr.plot(mcmc(BT_Run_Ets1$Log$results$q01))

#Plot & Summarize Ets1
Node_Values_Ets1<-BT_Run_Ets1$Log[2]
Node_Values_Ets1<-Node_Values_Ets1$results[,6:ncol(Node_Values_Ets1$results)]
Mean_Nodes_Ets1<-colMeans(Node_Values_Ets1)
Node_List_Ets1<-matrix(nrow=length(Mean_Nodes_Ets1), ncol=2)
for (a in seq(from=1, to=length(Mean_Nodes_Ets1), by=2)){
  Node_Ets1<-c(Mean_Nodes_Ets1[a], Mean_Nodes_Ets1[a+1])
  Node_List_Ets1[a,]<-c(Mean_Nodes_Ets1[a], Mean_Nodes_Ets1[a+1])
  pie(Node_Ets1)
}

Node_List_Ets1<-na.omit(Node_List_Ets1)
geoscalePhylo(tree=ladderize(tree, right=TRUE), 
              quat.rm=TRUE, units=c("Period"), boxes="Period",
              tick.scale="Period", cex.tip=1, cex.age=.5, cex.ts=.7, width=2,
              label.offset=0, x.lim=c(-90, 516), lwd=1.5)
nodelabels(pie=Node_List_Ets1, piecol=c("green", "grey"))
Numbers<-(1:11)
Numbers<-as.character(Numbers)
nodelabels(Numbers, adj=-2)
tiplabels(pch=21, col="black",  bg=Ets1_col, cex=2.3)
write.xlsx(Node_List_Ets1, "/Users/jeffreythompson/Research/Embryos_and_Ancestors/ASR/BayesTraits/Real_Final/Sensitivity/0_0.2/Ets1/Ets1_States.xls")

pdf("/Users/jeffreythompson/Research/Embryos_and_Ancestors/ASR/BayesTraits/Real_Final/Sensitivity/0_0.2/Ets1/Ets1_Bayes_Traits.pdf", width=10, height=7)
geoscalePhylo(tree=ladderize(tree, right=TRUE), 
              quat.rm=TRUE, units=c("Period"), boxes="Period",
              tick.scale="Period", cex.tip=1, cex.age=.5, cex.ts=.7, width=2,
              label.offset=0, x.lim=c(-90, 516), lwd=1.5)
nodelabels(pie=Node_List_Ets1, piecol=c("green", "grey"))
Numbers<-(1:11)
Numbers<-as.character(Numbers)
nodelabels(Numbers, adj=-2)
tiplabels(pch=21, col="black",  bg=Ets1_col, cex=2.3)
dev.off()


#Use Bayes Factor to Test Hypotheses regarding Ets1 expression at asterozoan MRCA
BT_Asterozoan_Zero_Ets1<-bayestraits(Ets1_Dataframe, trees_SCM, Command_Fixed_State_Zero)

BT_Asterozoan_One_Ets1<-bayestraits(Ets1_Dataframe, trees_SCM, Command_Fixed_State_One)

Log_BF_Asterozoan_Ets1<-2*(BT_Asterozoan_Zero_Ets1$Stones$logMarLH-BT_Asterozoan_One_Ets1$Stones$logMarLH)

#For Tbr
Tbr_Dataframe<-as.data.frame(Tbr_Data)

BT_One_rate_Tbr<-bayestraits(Tbr_Dataframe, trees_SCM, Command_One_Rate)

BT_Two_rates_Tbr<-bayestraits(Tbr_Dataframe, trees_SCM, Command_Two_Rates)

Log_BF_Rates_Tbr<-2*(BT_Two_rates_Tbr$Stones$logMarLH-BT_One_rate_Tbr$Stones$logMarLH)

#Reconstruct Ancestral States for Tbr
BT_Run_Tbr<-bayestraits(Tbr_Dataframe, trees_SCM, Command)
traceplot(mcmc(BT_Run_Tbr$Log$results$q01))
densplot(mcmc(BT_Run_Tbr$Log$results$q01))
autocorr.plot(mcmc(BT_Run_Tbr$Log$results$q01))

#Plot & Summarize Tbr Results
Node_Values_Tbr<-BT_Run_Tbr$Log[2]
Node_Values_Tbr<-Node_Values_Tbr$results[,16:ncol(Node_Values_Tbr$results)]
Mean_Nodes_Tbr<-colMeans(Node_Values_Tbr)
Node_List_Tbr<-matrix(nrow=length(Mean_Nodes_Tbr), ncol=4)
for (a in seq(from=1, to=length(Mean_Nodes_Tbr), by=4)){
  Node_Tbr<-c(Mean_Nodes_Tbr[a], Mean_Nodes_Tbr[a+1], Mean_Nodes_Tbr[a+2], Mean_Nodes_Tbr[a+3])
  Node_List_Tbr[a,]<-c(Mean_Nodes_Tbr[a], Mean_Nodes_Tbr[a+1], Mean_Nodes_Tbr[a+2], Mean_Nodes_Tbr[a+3])
  pie(Node)
}
Node_List_Tbr<-na.omit(Node_List_Tbr)
geoscalePhylo(tree=ladderize(tree, right=TRUE), 
              quat.rm=TRUE, units=c("Period"), boxes="Period",
              tick.scale="Period", cex.tip=1, cex.age=.5, cex.ts=.7, width=2,
              label.offset=0, x.lim=c(-90, 516), lwd=1.5)
nodelabels(pie=Node_List_Tbr, piecol=c("green", "blue", "yellow", "red"))
Numbers<-(1:11)
Numbers<-as.character(Numbers)
nodelabels(Numbers, adj=-2)
tiplabels(pch=21, col="black",  bg=Tbr_col, cex=2.3)
write.xlsx(Node_List_Tbr, "/Users/jeffreythompson/Research/Embryos_and_Ancestors/ASR/BayesTraits/Real_Final/Sensitivity/0_0.2/Tbr/Tbr_States.xls")

pdf("/Users/jeffreythompson/Research/Embryos_and_Ancestors/ASR/BayesTraits/Real_Final/Sensitivity/0_0.2/Tbr/Tbr_Bayes_Traits.pdf", width=10, height=7)
geoscalePhylo(tree=ladderize(tree, right=TRUE), 
              quat.rm=TRUE, units=c("Period"), boxes="Period",
              tick.scale="Period", cex.tip=1, cex.age=.5, cex.ts=.7, width=2,
              label.offset=0, x.lim=c(-90, 516), lwd=1.5)
nodelabels(pie=Node_List_Tbr, piecol=c("green", "blue", "yellow", "red"))
Numbers<-(1:11)
Numbers<-as.character(Numbers)
nodelabels(Numbers, adj=-2)
tiplabels(pch=21, col="black",  bg=Tbr_col, cex=2.3)
dev.off()


#For VegfR
#Prune Trees to Taxa with VegfR
Veg_Tree<-lapply(trees_SCM, drop.tip, tip=c("Echinocardium", "Echinodiscus", "Prionocidaris", "Apostichopus", "Holothuria"))
class(Veg_Tree)<-"multiPhylo"
write.nexus(Veg_Tree, file="/Users/jeffreythompson/Research/Embryos_and_Ancestors/ASR/BayesTraits/Final/Midway")
Veg_Tree<-read.nexus(file="/Users/jeffreythompson/Research/Embryos_and_Ancestors/ASR/BayesTraits/Final/Midway")

VegF_names<-matrix(c("Lytechinus", "Paracentrotus", "Strongylocentrotus", "Eucidaris", "Amphiura", "Amphipholis", "Patiria"))
VegF_values<-matrix(c(0,0,0,0,0,0,1), ncol=1, nrow=7)
VegF_col<-character(length(Veg_Tree$tip.label))
VegF_col[VegF_values=="0"]<-"darkgoldenrod1"
VegF_col[VegF_values=="1"]<-"black"
VegF_Data<-cbind(VegF_names, VegF_values)
VegF_Dataset<-as.data.frame(VegF_Data)

#Use Bayes Factor to Compare One Rate to Two Rate Models
BT_One_rate_VegfR<-bayestraits(VegF_Dataset, Veg_Tree, Command_One_Rate)

BT_Two_rates_VegfR<-bayestraits(VegF_Dataset, Veg_Tree, Command_Two_Rates)

Log_BF_Rates_VegfR<-2*(BT_Two_rates_VegfR$Stones$logMarLH-BT_One_rate_VegfR$Stones$logMarLH)


#Reconstruct Ancestral States VegfR
Command_VegfR<-c("1", "2", "ScaleTrees",
                 "Iterations 100000000",
                 "Burnin 20000000",
                 "Resall q10",
                 "PriorAll uniform 0 .2",
                 "AddTag Echinoids Strongylocentrotus Eucidaris",
                 "AddMRCA MRCAEchinoids Echinoids",
                 "AddTag PurpLytechinus Strongylocentrotus Lytechinus",
                 "AddMRCA MRCAPurpLytechinus PurpLytechinus",
                 "AddTag PurpParacentrotus Paracentrotus Strongylocentrotus",
                 "AddMRCA MRCAPurpParacentrotus PurpParacentrotus",
                 "AddTag Asterozoans Patiria Amphipholis",
                 "AddMRCA MRCAAsterozoans Asterozoans",
                 "AddTag Ophiuroids Amphiura Amphipholis",
                 "AddMRCA MRCAOphiuroids Ophiuroids"
)




BT_Run_VegfR<-bayestraits(VegF_Dataset, Veg_Tree, Command_VegfR)
traceplot(mcmc(BT_Run_VegfR$Log$results$q01))
densplot(mcmc(BT_Run_VegfR$Log$results$q01))
autocorr.plot(mcmc(BT_Run_VegfR$Log$results$q01))

#Plot & Summarize VegfR Results
Node_Values_VegF<-BT_Run_VegfR$Log[2]
Node_Values_VegF<-Node_Values_VegF$results[,6:ncol(Node_Values_VegF$results)]
Mean_Nodes_VegF<-colMeans(Node_Values_VegF)
Node_List_VegF<-matrix(nrow=length(Mean_Nodes_VegF), ncol=2)
for (a in seq(from=1, to=length(Mean_Nodes_VegF), by=2)){
  Node_VegF<-c(Mean_Nodes_VegF[a], Mean_Nodes_VegF[a+1])
  Node_List_VegF[a,]<-c(Mean_Nodes_VegF[a], Mean_Nodes_VegF[a+1])
  pie(Node_VegF)
}
Node_List_VegF<-na.omit(Node_List_VegF)
Tree_VegF<-drop.tip(tree, tip=c("Echinocardium", "Echinodiscus", "Prionocidaris", "Apostichopus", "Holothuria"))
geoscalePhylo(tree=ladderize(Tree_VegF, right=TRUE), 
              quat.rm=TRUE, units=c("Period"), boxes="Period",
              tick.scale="no", cex.tip=1, cex.age=.5, cex.ts=.7, width=2,
              label.offset=0, x.lim=c(-90, 516), lwd=1.5)
nodelabels(pie=Node_List_VegF, piecol=c("darkgoldenrod1", "black"))
Numbers<-(1:6)
Numbers<-as.character(Numbers)
nodelabels(Numbers, adj=-2)
tiplabels(pch=21, col="black",  bg=VegF_col, cex=2.3)
write.xlsx(Node_List_VegF, "/Users/jeffreythompson/Research/Embryos_and_Ancestors/ASR/BayesTraits/Real_Final/Sensitivity/0_0.2/VegfR/VegfR_States.xls")

pdf("/Users/jeffreythompson/Research/Embryos_and_Ancestors/ASR/BayesTraits/Real_Final/Sensitivity/0_0.2/VegfR/VegfR_Bayes_Traits.pdf", width=10, height=7)
geoscalePhylo(tree=ladderize(Tree_VegF, right=TRUE), 
              quat.rm=TRUE, units=c("Period"), boxes="Period",
              tick.scale="no", cex.tip=1, cex.age=.5, cex.ts=.7, width=2,
              label.offset=0, x.lim=c(-90, 516), lwd=1.5)
nodelabels(pie=Node_List_VegF, piecol=c("darkgoldenrod1", "black"), cex=.8)
Numbers<-(1:6)
Numbers<-as.character(Numbers)
nodelabels(Numbers, adj=-2)
tiplabels(pch=21, col="black",  bg=VegF_col, cex=2.3)
dev.off()


#Use Bayes Factor to Test Hypotheses regarding VegfR expression at asterozoan MRCA
BT_Asterozoan_Zero_VegfR<-bayestraits(VegF_Dataset, Veg_Tree, Command_Fixed_State_Zero)

BT_Asterozoan_One_VegfR<<-bayestraits(VegF_Dataset, Veg_Tree, Command_Fixed_State_One)

Log_BF_Asterozoan_VegfR<-2*(BT_Asterozoan_Zero_VegfR$Stones$logMarLH-BT_Asterozoan_One_VegfR$Stones$logMarLH)

#For Erg
#Prune Trees to Taxa with Erg
Erg_Tree<-lapply(trees_SCM, drop.tip, tip=c("Echinocardium", "Paracentrotus", "Echinodiscus", "Prionocidaris", "Apostichopus", "Amphipholis" ))
class(Erg_Tree)<-"multiPhylo"
write.nexus(Erg_Tree, file="/Users/jeffreythompson/Research/Embryos_and_Ancestors/ASR/BayesTraits/Final/Midway")
Erg_Tree<-read.nexus(file="/Users/jeffreythompson/Research/Embryos_and_Ancestors/ASR/BayesTraits/Final/Midway")

Erg_names<-matrix(c("Lytechinus", "Strongylocentrotus", "Eucidaris", "Holothuria", "Amphiura", "Patiria"))
Erg_values<-matrix(c(0,0,0,0,0,1), ncol=1, nrow=6)
Erg_values_for_color<-matrix(c(0,0,0,0,1,0), ncol=1, nrow=6)
Erg_col<-character(length(Erg_Tree$tip.label))
Erg_col[Erg_values_for_color=="0"]<-"green"
Erg_col[Erg_values_for_color=="1"]<-"grey"
Erg_Data<-cbind(Erg_names, Erg_values)
Erg_Dataset<-as.data.frame(Erg_Data)

#Use Bayes Factor to Compare One Rate to Two Rate Models. 
BT_One_rate_Erg<-bayestraits(Erg_Dataset, Erg_Tree, Command_One_Rate)

BT_Two_rates_Erg<-bayestraits(Erg_Dataset, Erg_Tree, Command_Two_Rates)

Log_BF_Rates_Erg<-2*(BT_Two_rates_Erg$Stones$logMarLH-BT_One_rate_Erg$Stones$logMarLH)

#Reconstruct Ancestral States Erg
Command_Erg<-c("1", "2", "ScaleTrees",
                 "Iterations 100000000",
                 "Burnin 20000000",
                 "Resall q10",
                 "PriorAll uniform 0 200",
                 "AddTag Echinozoans Strongylocentrotus Holothuria",
                 "AddMRCA MRCAEchinozoans Echinozoans",
                 "AddTag Echinoids Strongylocentrotus Eucidaris",
                 "AddMRCA MRCAEchinoids Echinoids",
                 "AddTag PurpLytechinus Strongylocentrotus Lytechinus",
                 "AddMRCA MRCAPurpLytechinus PurpLytechinus",
                 "AddTag Asterozoans Patiria Amphiura",
                 "AddMRCA MRCAAsterozoans Asterozoans"
)

BT_Run_Erg<-bayestraits(Erg_Dataset, Erg_Tree, Command_Erg)
traceplot(mcmc(BT_Run_Erg$Log$results$q01))
densplot(mcmc(BT_Run_Erg$Log$results$q01))
autocorr.plot(mcmc(BT_Run_Erg$Log$results$q01))


#Plot & Summarize Erg Results
Node_Values_Erg<-BT_Run_Erg$Log[2]
Node_Values_Erg<-Node_Values_Erg$results[,6:ncol(Node_Values_Erg$results)]
Mean_Nodes_Erg<-colMeans(Node_Values_Erg)
Node_List_Erg<-matrix(nrow=length(Mean_Nodes_Erg), ncol=2)
for (a in seq(from=1, to=length(Mean_Nodes_Erg), by=2)){
  Node_Erg<-c(Mean_Nodes_Erg[a], Mean_Nodes_Erg[a+1])
  Node_List_Erg[a,]<-c(Mean_Nodes_Erg[a], Mean_Nodes_Erg[a+1])
  pie(Node_Erg)
}
Node_List_Erg<-na.omit(Node_List_Erg)
Tree_Erg<-drop.tip(tree, tip=c("Echinocardium", "Paracentrotus", "Echinodiscus", "Prionocidaris", "Apostichopus", "Amphipholis" ))
geoscalePhylo(tree=ladderize(Tree_Erg, right=TRUE), 
              quat.rm=TRUE, units=c("Period"), boxes="Period",
              tick.scale="no", cex.tip=1, cex.age=.5, cex.ts=.7, width=2,
              label.offset=0, x.lim=c(-90, 516), lwd=1.5)
nodelabels(pie=Node_List_Erg, piecol=c("green", "grey"))
Numbers<-(1:5)
Numbers<-as.character(Numbers)
nodelabels(Numbers, adj=-2)
tiplabels(pch=21, col="black",  bg=Erg_col, cex=2.3)
write.xlsx(Node_List_Erg, "/Users/jeffreythompson/Research/Embryos_and_Ancestors/ASR/BayesTraits/Real_Final/Sensitivity/0_200/Erg/Erg_States.xls")

pdf("/Users/jeffreythompson/Research/Embryos_and_Ancestors/ASR/BayesTraits/Real_Final/Sensitivity/0_200/Erg/Erg_Bayes_Traits.pdf", width=10, height=7)
geoscalePhylo(tree=ladderize(Tree_Erg, right=TRUE), 
              quat.rm=TRUE, units=c("Period"), boxes="Period",
              tick.scale="no", cex.tip=1, cex.age=.5, cex.ts=.7, width=2,
              label.offset=0, x.lim=c(-90, 516), lwd=1.5)
nodelabels(pie=Node_List_Erg, piecol=c("green", "grey"), cex=.8)
Numbers<-(1:5)
Numbers<-as.character(Numbers)
nodelabels(Numbers, adj=-2)
tiplabels(pch=21, col="black",  bg=Erg_col, cex=2.3)
dev.off()

#Use Bayes Factor to Test Hypotheses regarding Erg expression at asterozoan MRCA
BT_Asterozoan_Zero_Erg<-bayestraits(Erg_Dataset, Erg_Tree, Command_Fixed_State_Zero)

BT_Asterozoan_One_Erg<<-bayestraits(Erg_Dataset, Erg_Tree, Command_Fixed_State_One)

Log_BF_Asterozoan_Erg<-2*(BT_Asterozoan_Zero_Erg$Stones$logMarLH-BT_Asterozoan_One_Erg$Stones$logMarLH)




#Print Final Figure with all ASRs
pdf("/Users/jeffreythompson/Research/Embryos_and_Ancestors/ASR/BayesTraits/Real_Final/Combo_Tree/Final_Tree_All.pdf", width=10, height=7)
geoscalePhylo(tree=ladderize(tree, right=TRUE), 
              quat.rm=TRUE, units=c("Period"), boxes="Period",
              tick.scale="Period", cex.tip=1, cex.age=.5, cex.ts=.7, width=2,
              label.offset=0, x.lim=c(-90, 516), lwd=1.5)
Numbers<-(1:11)
Numbers<-as.character(Numbers)
nodelabels(Numbers, adj=-4)
nodelabels(pie=Node_List_Tbr, piecol=c("green", "blue", "yellow", "red"), adj=-25, cex=.8)
nodelabels(pie=Node_List_Ets1, piecol=c("green", "grey"), adj=0, cex=.8)
nodelabels(pie=Node_List_Alx, piecol=c("blue", "grey"), adj=25, cex=.8)
tiplabels(pch=21, col="black",  bg=Tbr_col, cex=2.3)
tiplabels(pch=21, col="black",  bg=Ets1_col, cex=2.3, adj=25)
tiplabels(pch=21, col="black",  bg=Alx_col, cex=2.3, adj=50)
dev.off()

#Produce Table showing results from all Hypothesis tests at Asterozoan MRCA

Hypothesis_Tests<-data.frame(row.names=c("Alx1", "Ets1", "VegfR", "Erg"))

Hypothesis_Tests[1,1]<-Log_BF_Asterozoan_Alx1
Hypothesis_Tests[2,1]<-Log_BF_Asterozoan_Ets1
Hypothesis_Tests[3,1]<-Log_BF_Asterozoan_VegfR
Hypothesis_Tests[4,1]<-Log_BF_Asterozoan_Erg

write.xlsx(Hypothesis_Tests, "/Users/jeffreythompson/Research/Embryos_and_Ancestors/ASR/BayesTraits/Real_Final/Bayes_Factors/Hypothesis_Tests.xls")




#Produce Table showing results from model comparisons for 1 or two rate models

Model_Tests<-data.frame(row.names=c("Alx1", "Ets1", "VegfR", "Tbr", "Erg"))

Model_Tests[1,1]<-Log_BF_Rates_Alx1
Model_Tests[2,1]<-Log_BF_Rates_Ets
Model_Tests[3,1]<-Log_BF_Rates_VegfR
Model_Tests[4,1]<-Log_BF_Rates_Tbr
Model_Tests[5,1]<-Log_BF_Rates_Erg
write.xlsx(Model_Tests, "/Users/jeffreythompson/Research/Embryos_and_Ancestors/ASR/BayesTraits/Real_Final/Bayes_Factors/Model_Comparison.xls")




#Run sensitivity analysis pruning taxa to one taxon per family for alx1, ets1, and tbr

One_Per_Tree<-lapply(trees_SCM, drop.tip, tip=c("Echinocardium", "Echinodiscus", "Lytechinus", "Paracentrotus", "Prionocidaris", "Eucidaris", "Holothuria", "Amphipholis"))
class(One_Per_Tree)<-"multiPhylo"

write.nexus(One_Per_Tree, file="/Users/jeffreythompson/Research/Embryos_and_Ancestors/ASR/BayesTraits/Final/One_Per_Tree")
One_Per_Tree<-read.nexus(file="/Users/jeffreythompson/Research/Embryos_and_Ancestors/ASR/BayesTraits/Final/One_Per_Tree")



Command_One_Per<-c("1", "2", "ScaleTrees",
                 "Iterations 100000000",
                 "Burnin 20000000",
                 "Resall q10",
                 "PriorAll uniform 0 2",
                 "AddTag Echinozoans Strongylocentrotus Apostichopus",
                 "AddMRCA MRCAEchinozoans Echinozoans",
                 "AddTag Asterozoans Patiria Amphiura",
                 "AddMRCA MRCAAsterozoans Asterozoans"
)

Tbr_One_Per<-Tbr_Dataframe[-c(1,2,4,5,6,7,9,11),]
Tree_One_Per<-drop.tip(tree, tip=c("Echinocardium", "Echinodiscus", "Lytechinus", "Paracentrotus", "Prionocidaris", "Eucidaris", "Holothuria", "Amphipholis"))

#Alx_One_Per
Alx_One_Per<-Alx_Dataframe[-c(1,2,4,5,6,7,9,11),]
Alx_col_One_Per<-Alx_col[-c(1,2,4,5,6,7,9,11)]
BT_Run_One_Per_Alx1<-bayestraits(Alx_One_Per, One_Per_Tree, Command_One_Per)
traceplot(mcmc(BT_Run_One_Per_Alx1$Log$results$q01))
densplot(mcmc(BT_Run_One_Per_Alx1$Log$results$q01))
autocorr.plot(mcmc(BT_Run_One_Per_Alx1$Log$results$q01))

Node_Values_Alx_One_Per<-BT_Run_One_Per_Alx1$Log[2]
Node_Values_Alx_One_Per<-Node_Values_Alx_One_Per$results[,6:ncol(Node_Values_Alx_One_Per$results)]
Mean_Nodes_Alx_One_Per<-colMeans(Node_Values_Alx_One_Per)
Node_List_Alx_One_Per<-matrix(nrow=length(Mean_Nodes_Alx_One_Per), ncol=2)
for (a in seq(from=1, to=length(Mean_Nodes_Alx_One_Per), by=2)){
  Node<-c(Mean_Nodes_Alx_One_Per[a], Mean_Nodes_Alx_One_Per[a+1])
  Node_List_Alx_One_Per[a,]<-c(Mean_Nodes_Alx_One_Per[a], Mean_Nodes_Alx_One_Per[a+1])
  pie(Node)
}

pdf("/Users/jeffreythompson/Research/Embryos_and_Ancestors/ASR/BayesTraits/Real_Final/Sensitivity/Pruned/Alx1/Alx1_Bayes.pdf", width=10, height=7)
Node_List_Alx_One_Per<-na.omit(Node_List_Alx_One_Per)
geoscalePhylo(tree=ladderize(Tree_One_Per, right=TRUE), 
              quat.rm=TRUE, units=c("Period"), boxes="Period",
              tick.scale="Period", cex.tip=1, cex.age=.5, cex.ts=.7, width=2,
              label.offset=0, x.lim=c(-90, 516), lwd=1.5)
nodelabels(pie=Node_List_Alx_One_Per, piecol=c("blue", "grey"))
Numbers<-(1:3)
Numbers<-as.character(Numbers)
nodelabels(Numbers, adj=-2)
tiplabels(pch=21, col="black",  bg=Alx_col_One_Per, cex=2.3)
dev.off()
write.xlsx(Node_List_Alx_One_Per, "/Users/jeffreythompson/Research/Embryos_and_Ancestors/ASR/BayesTraits/Real_Final/Sensitivity/Pruned/Alx1/Alx1_States.xls")


#Ets1_One_Per
Ets1_One_Per<-Ets1_Dataframe[-c(1,2,4,5,6,7,9,11),]
Ets_col_One_Per<-Ets1_col[-c(1,2,4,5,6,7,9,11)]
BT_Run_One_Per_Ets1<-bayestraits(Ets1_One_Per, One_Per_Tree, Command_One_Per)
traceplot(mcmc(BT_Run_One_Per_Ets1$Log$results$q01))
densplot(mcmc(BT_Run_One_Per_Ets1$Log$results$q01))
autocorr.plot(mcmc(BT_Run_One_Per_Ets1$Log$results$q01))

Node_Values_Ets1_One_Per<-BT_Run_One_Per_Ets1$Log[2]
Node_Values_Ets1_One_Per<-Node_Values_Ets1_One_Per$results[,6:ncol(Node_Values_Ets1_One_Per$results)]
Mean_Nodes_Ets1_One_Per<-colMeans(Node_Values_Ets1_One_Per)
Node_List_Ets1_One_Per<-matrix(nrow=length(Mean_Nodes_Ets1_One_Per), ncol=2)
for (a in seq(from=1, to=length(Mean_Nodes_Ets1_One_Per), by=2)){
  Node<-c(Mean_Nodes_Ets1_One_Per[a], Mean_Nodes_Ets1_One_Per[a+1])
  Node_List_Ets1_One_Per[a,]<-c(Mean_Nodes_Ets1_One_Per[a], Mean_Nodes_Ets1_One_Per[a+1])
  pie(Node)
}


pdf("/Users/jeffreythompson/Research/Embryos_and_Ancestors/ASR/BayesTraits/Real_Final/Sensitivity/Pruned/Ets1/Ets1_Bayes.pdf", width=10, height=7)
Node_List_Ets1_One_Per<-na.omit(Node_List_Ets1_One_Per)
geoscalePhylo(tree=ladderize(Tree_One_Per, right=TRUE), 
              quat.rm=TRUE, units=c("Period"), boxes="Period",
              tick.scale="Period", cex.tip=1, cex.age=.5, cex.ts=.7, width=2,
              label.offset=0, x.lim=c(-90, 516), lwd=1.5)
nodelabels(pie=Node_List_Ets1_One_Per, piecol=c("green", "grey"))
Numbers<-(1:3)
Numbers<-as.character(Numbers)
nodelabels(Numbers, adj=-2)
tiplabels(pch=21, col="black",  bg=Ets_col_One_Per, cex=2.3)
dev.off()
write.xlsx(Node_List_Ets1_One_Per, "/Users/jeffreythompson/Research/Embryos_and_Ancestors/ASR/BayesTraits/Real_Final/Sensitivity/Pruned/Ets1/Ets1_States.xls")


#Tbr_One_Per
BT_Run_One_Per_Tbr<-bayestraits(Alx_One_Per, One_Per_Tree, Command_One_Per)
Tbr_One_Per<-Tbr_Dataframe[-c(1,2,4,5,6,7,9,11),]
Tbr_col_One_Per<-Tbr_col[-c(1,2,4,5,6,7,9,11)]
BT_Run_One_Per_Tbr<-bayestraits(Tbr_One_Per, One_Per_Tree, Command_One_Per)
traceplot(mcmc(BT_Run_One_Per_Tbr$Log$results$q01))
densplot(mcmc(BT_Run_One_Per_Tbr$Log$results$q01))
autocorr.plot(mcmc(BT_Run_One_Per_Tbr$Log$results$q01))

Node_Values_Tbr_One_Per<-BT_Run_One_Per_Tbr$Log[2]
Node_Values_Tbr_One_Per<-Node_Values_Tbr_One_Per$results[,16:ncol(Node_Values_Tbr_One_Per$results)]
Mean_Nodes_Tbr_One_Per<-colMeans(Node_Values_Tbr_One_Per)
Node_List_Tbr_One_Per<-matrix(nrow=length(Mean_Nodes_Tbr_One_Per), ncol=4)
for (a in seq(from=1, to=length(Mean_Nodes_Tbr_One_Per), by=4)){
  Node<-c(Mean_Nodes_Tbr_One_Per[a], Mean_Nodes_Tbr_One_Per[a+1])
  Node_List_Tbr_One_Per[a,]<-c(Mean_Nodes_Tbr_One_Per[a], Mean_Nodes_Tbr_One_Per[a+1], Mean_Nodes_Tbr_One_Per[a+2], Mean_Nodes_Tbr_One_Per[a+3])
  pie(Node)
}

pdf("/Users/jeffreythompson/Research/Embryos_and_Ancestors/ASR/BayesTraits/Real_Final/Sensitivity/Pruned/Tbr/Tbr_Bayes.pdf", width=10, height=7)
Node_List_Tbr_One_Per<-na.omit(Node_List_Tbr_One_Per)
geoscalePhylo(tree=ladderize(Tree_One_Per, right=TRUE), 
              quat.rm=TRUE, units=c("Period"), boxes="Period",
              tick.scale="Period", cex.tip=1, cex.age=.5, cex.ts=.7, width=2,
              label.offset=0, x.lim=c(-90, 516), lwd=1.5)
nodelabels(pie=Node_List_Tbr_One_Per, piecol=c("green", "blue", "yellow", "red"))
Numbers<-(1:3)
Numbers<-as.character(Numbers)
nodelabels(Numbers, adj=-2)
tiplabels(pch=21, col="black",  bg=Tbr_col_One_Per, cex=2.3)
dev.off()
write.xlsx(Node_List_Tbr_One_Per, "/Users/jeffreythompson/Research/Embryos_and_Ancestors/ASR/BayesTraits/Real_Final/Sensitivity/Pruned/Tbr/Tbr_States.xls")



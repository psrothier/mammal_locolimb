
#dat <- read.csv("data_2analyses.csv", sep=",", header=T, na.strings = "NA")
#rownames(dat) <- dat$Species

library(ape); library(geiger)
library(evobiR); library(phytools); library(reshape2); library(mvMORPH)
library(ggplot2); library(dplyr); library(plotly)
library(psych)
library(cowplot)
library(dplyr); library(plotrix)

#mytree <- read.tree("mytree.nex")
#plot(mytree, show.tip.label = FALSE); axisPhylo()
#mytree <- ladderize(mytree)
#plot(mytree, show.tip.label = FALSE); axisPhylo()

#dat <- ReorderData(mytree, dat, taxa.names="rownames")
#mydata <- na.omit(dat) # 797 species
#str(mydata)

#name.check(mytree, mydata)
#to_remove1 <- mytree$tip.label[!mytree$tip.label%in%rownames(mydata)]
#tree <- drop.tip(mytree, to_remove1)
#name.check(tree, mydata)

#write.tree(tree, file="pruned_MCC_tree.nex")
#dataset <- list(mydata=mydata, tree=tree)
#save(dataset, file="dataset.Rdata")



##=========================================================================================================================##
#                                                                                                                           #   
#### Pt I) Dataset preparation                                                                                           ####
#                                                                                                                           #
##=========================================================================================================================##

#load("dataset.Rdata")
#mydata <- dataset$mydata
#tree <- dataset$tree
#name.check(tree, mydata)

#mydata$Transf_Mass <- log10(mydata$Body_Mass^(1/3)) # transform bm in linear scale
#lgm <- log10(apply(mydata[,10:31], 1, geometric.mean))
#full.limb <- as.matrix(log10(mydata[,10:30]))
#limb <- as.matrix(full.limb[,c(1,3,5,6,8,10,11,13,15,16,18,20,21)]) # taking a subset

#bodymass <- mydata$Transf_Mass
#names(bodymass) <- row.names(mydata)



##=========================================================================================================================##
#                                                                                                                           #   
####  Pt II) phyloPCA                                                                                                    ####
#                                                                                                                           #
##=========================================================================================================================##

#data.resid <- phyl.resid(tree, lgm, limb, method='lambda')
#residM<-data.frame(data.resid$resid)
#pca <- phyl.pca(tree, residM, mode='cov') #residual pca
#listdata <- list(tree=tree, bodymass=bodymass, residM=residM, pca=pca, mydata=mydata, limb=limb, lgm=lgm)

#save(listdata, file="listdata.RData")

load("listdata.RData")
tree <- listdata$tree
bodymass <- listdata$bodymass # log 10 body mass transformed into linear scale
limb <- listdata$limb # log 10 limb traits
residM<- listdata$residM
pca <- listdata$pca
summary(pca)
mydata <- listdata$mydata # raw values
locomotion <- mydata$Locomotion
names(locomotion) <- rownames(mydata)

cols<-setNames(c("blue","forestgreen", "azure3","darkred", "black","olivedrab1", "skyblue", "red", "violet", "#FFC055"),
               c("aquatic","arboreal","flight","fossorial","gliding","scansorial","semiaquatic","semifossorial","terrestrialbip","terrestrialquad"))

##==========================================================================================================##
#### phyANOVA                                                                                             ####
##==========================================================================================================##
?phylANOVA
aov <- phylANOVA(tree, locomotion, bodymass, nsim=10000)
bm.loc <- data.frame(bodymass=bodymass, locomotion=mydata$Locomotion)
mean.mass <- as.data.frame(mydata %>% group_by(Locomotion) %>% summarise(mean(Body_Mass)))
mean.mass.kg <- round(mean.mass[,2]/1000,digits=3)
median.mass <- as.data.frame(mydata %>% group_by(Locomotion) %>% summarise(median(Body_Mass)))
median.mass.kg <- round(median.mass[,2]/1000, digits=3)


ggplot(bm.loc, aes(x=locomotion, y=bodymass, fill=locomotion, color=locomotion))+
  theme_classic()+
  geom_point(size=1, shape=21, position=position_jitter(width=0.13, height=0.01), alpha=0.4) +
  scale_fill_manual(values=cols)+
  scale_color_manual(values=cols)+
  geom_boxplot(width=0.4, outlier.colour = NA, fill=cols, color="black", alpha=0.2) +
  xlab(NULL)+
  ylab("log10 transf. body mass(g)")+
  theme(legend.position="none")+
  theme(axis.text.y = element_text(size = 7), axis.title.y = element_text(size=9),
        axis.text.x = element_text(size = 7, angle=45, hjust=1))







##==========================================================================================================##
#### Plot PCA--                                                                                           ####
##==========================================================================================================##


Pcscores <- as.data.frame(pca$S)
Pcscores <- Pcscores[match(tree$tip.label, rownames(Pcscores)),]
Pcscores$Group <- mydata$Group
Pcscores$Locomotion <- mydata$Locomotion


## creating object for BayesTraits --------------- ##
#PC95<-Pcscores[,1:4]
#write.table(PC95*1000,
#            file = "./Pc95_tout.txt",
#            quote = FALSE, col.names = FALSE)


diag(pca$Eval)
pca$L

Pcscores <- as.data.frame(Pcscores)
phylsig.pc <- vector()
for(i in 1:4){
phylsig.pc[i] <- phylosig(tree, Pcscores[,i], method="lambda")
}




####################### !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# include phylogenetic sigma function for body mass and pc scores 
####################################################################


library(ggplot2)

## PC1 x PC2 --------------- ##

plot_tr <- ggplot(Pcscores, aes(PC1, PC2, group=Group, text=rownames(Pcscores))) + 
  geom_point(aes(color=Locomotion), size=2, alpha=0.5) +
  scale_color_manual(values=cols)+
  theme_minimal()+
  xlab("pPC1 (27.4 %)")+
  ylab("pPC2 (17.8 %)")+
  guides(col=guide_legend(ncol=1))+
  theme(axis.text.y = element_text(size = 10), axis.title.y = element_text(size=12),
        axis.text.x = element_text(size = 10), axis.title.x = element_text(size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=10))+
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") 

plotly::ggplotly(plot_tr, tooltip=c("text", "x", "y", "group"))



# with ellipse

ggplot(Pcscores, aes(PC1, PC2, color=Locomotion)) + 
  geom_point(aes(color=Locomotion), size=2, alpha=0.5) +
  scale_color_manual(values=cols)+
  theme_minimal()+
  xlab("pPC1 (27.4 %)")+
  ylab("pPC2 (17.8 %)")+
  guides(col=guide_legend(ncol=1))+
  theme(axis.text.y = element_text(size = 10), axis.title.y = element_text(size=12),
        axis.text.x = element_text(size = 10), axis.title.x = element_text(size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=10))+
  stat_ellipse(geom = "polygon",
             aes(fill = Locomotion),
             alpha = 0.2) + scale_fill_manual(values=cols)+
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") 




## PC1 x bodymass --------------- ##

Pcscores$bodymass <- bodymass

plot_tr2 <- ggplot(Pcscores, aes(bodymass, PC1, group=Group, text=rownames(Pcscores))) + 
  geom_point(aes(color=Locomotion), size=2, alpha=0.5) +
  scale_color_manual(values=cols)+
  theme_classic()+
  xlab("Body mass")+
  ylab("PC1 (24.4 %)")+
  guides(col=guide_legend(ncol=1))+
  theme(axis.text.y = element_text(size = 10), axis.title.y = element_text(size=12),
        axis.text.x = element_text(size = 10), axis.title.x = element_text(size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=10))
plotly::ggplotly(plot_tr2, tooltip=c("text", "x", "y", "group"))





## PC1 x PC3 --------------- ##

plot_tr3 <- ggplot(Pcscores, aes(x=PC1, y=PC3, group=Group, text=rownames(Pcscores))) + # multiplying scores *-1 to help visualization
  geom_point(aes(color=Locomotion), size=2, alpha=0.5) +
  scale_color_manual(values=cols)+
  theme_classic()+
  xlab("PC1")+
  ylab("PC3")+
  guides(col=guide_legend(ncol=1))+
  theme(axis.text.y = element_text(size = 10), axis.title.y = element_text(size=12),
        axis.text.x = element_text(size = 10), axis.title.x = element_text(size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=10))
plotly::ggplotly(plot_tr3, tooltip=c("text", "x", "y", "group"))


ggplot(Pcscores, aes(PC1, color=Locomotion, fill=Locomotion)) +
  theme_classic()+
  geom_density(alpha=0.3)+
  scale_color_manual(values=cols)+
  scale_fill_manual(values=cols)

ggplot(Pcscores, aes(PC2, color=Locomotion, fill=Locomotion)) +
  theme_classic()+
  geom_density(alpha=0.3)+
  scale_color_manual(values=cols)+
  scale_fill_manual(values=cols)

ggplot(Pcscores, aes(bodymass, color=Locomotion, fill=Locomotion)) +
  theme_classic()+
  geom_density(alpha=0.3)+
  scale_color_manual(values=cols)+
  scale_fill_manual(values=cols)


diag(pca$Eval)



## 3D plot --------------- ##


library(plotly)
plot_ly(Pcscores, x=~PC1, y=~PC2, 
            z=~PC3, color=~Locomotion, colors=cols) %>%
  add_markers(size=1) 

library(pals)

unique(Pcscores$Group)

col2 <- cols25(21)
plot_ly(Pcscores, x=~PC1, y=~PC2, 
        z=~PC3, color=~Group, colors=col2) %>%
  add_markers(size=1) 


library(rgl)
library(car)
scatter3d(y=Pcscores$PC1, x=Pcscores$PC2, 
          z=Pcscores$PC3, xlab="pPC1", ylab="pPC2", zlab="pPC3", groups= as.factor(Pcscores$Locomotion),
          surface=FALSE, grid = FALSE, ellipsoid = TRUE, ellipsoid.alpha=0.2,
          surface.col = cols)




##=========================================================================================================================##
#                                                                                                                           #   
#### Pt III) fitMk                                                                                                       ####
#                                                                                                                           #
##=========================================================================================================================##

locomotion <- data.frame(locomotion=mydata$Locomotion)
rownames(locomotion) <- rownames(mydata)
head(locomotion)

locomotion.mode<-setNames(as.factor(locomotion[,1]),rownames(locomotion))
locomotion.mode



#-------------------------------------------------------------------------------------#
## a) comparing fitMk models
#-------------------------------------------------------------------------------------#

## we can load all fitMk simulations below

## Equal-rates (ER) --------------- ##

#fitER<-fitMk(tree,locomotion.mode,model="ER",pi="fitzjohn")
#fitER
#plot(fitER,width=TRUE,offset=0.03,color=TRUE,show.zeros=FALSE)


## All-rates different (ARD) --------------- ##

#fitARD<-fitMk(tree,locomotion.mode,model="ARD",pi="fitzjohn")
#fitARD
#plot(fitARD,width=TRUE,offset=0.03,color=TRUE)

## Symmetric transitions --------------- ##

#fitSYM<-fitMk(tree,locomotion.mode,model="SYM",pi="fitzjohn")
#fitSYM
#plot(fitSYM,width=TRUE,offset=0.03,color=TRUE)


## Directional model 1 --------------- ##

# ordered and directional transitions to specialization, departing from quadrupedal forms and passing by semi-specialized morphologies. return to previous state is not allowed.
#TQ --> semiaquat --> aquat
#TQ --> semifoss --> foss
#TQ --> scansorial --> arboreal --> gliding --> flight
#TQ --> TB


rnames <- c("aquatic","arboreal","flight","fossorial","gliding",
            "scansorial","semiaquatic","semifossorial","terrestrialbip","terrestrialquad")

cnames <- c("aquatic","arboreal","flight","fossorial","gliding",
            "scansorial","semiaquatic","semifossorial","terrestrialbip","terrestrialquad")


directional.model1<-matrix(c(
  0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,6,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,
  0,0,7,0,0,0,0,0,0,0,
  0,9,0,0,0,0,0,0,0,0,
  2,0,0,0,0,0,0,0,0,0,
  0,0,0,4,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,8,1,3,5,0),
  10,10,byrow=TRUE,   ### byrow=TRUE indicates that the matrix is being read in the same direction
  dimnames=list(rnames,cnames))
directional.model1

#fitdir1<-fitMk(tree,locomotion.mode,model=directional.model1,pi="fitzjohn")

#fitdir1
#plot(fitdir1,show.zeros=FALSE)



## Directional model 2 --------------- ##

# as directional model 1, but transitions are bi-directional:
#TQ <--> semiaquat <--> aquat
#TQ <--> semifoss <--> foss
#TQ< --> scansorial <--> arboreal <--> gliding <--> flight
#TQ <--> TB


directional.model2<-matrix(c(
  0,0,0,0,0,0,10,0,0,0,
  0,0,0,0,6,16,0,0,0,0,
  0,0,0,0,14,0,0,0,0,0,
  0,0,0,0,0,0,0,12,0,0,
  0,15,7,0,0,0,0,0,0,0,
  0,9,0,0,0,0,0,0,0,17,
  2,0,0,0,0,0,0,0,0,11,
  0,0,0,4,0,0,0,0,0,13,
  0,0,0,0,0,0,0,0,0,18,
  0,0,0,0,0,8,1,3,5,0),
  10,10,byrow=TRUE,   ### byrow=TRUE indicates that the matrix is being read in the same direction
  dimnames=list(rnames,cnames))
directional.model2

#fitdir2<-fitMk(tree,locomotion.mode,model=directional.model2,pi="fitzjohn")
#plot(fitdir2,show.zeros=FALSE,color=TRUE,width=TRUE)



## Directional model 3 --------------- ##

# as directional model 2, but transitions between semi-specialized states are allowed:
#TQ <--> semiaquat <--> aquat
#TQ <--> semifoss <--> foss
#TQ< --> scansorial <--> arboreal <--> gliding <--> flight
#TQ <--> TB
#Semiaquat <--> Semifoss
#Semiaquat <--> Scansorial
#Semifoss <--> Scansorial

directional.model3<-matrix(c(
  0,0,0,0,0,0,10,0,0,0,
  0,0,0,0,6,16,0,0,0,0,
  0,0,0,0,14,0,0,0,0,0,
  0,0,0,0,0,0,0,12,0,0,
  0,15,7,0,0,0,0,0,0,0,
  0,9,0,0,0,0,22,24,0,17,
  2,0,0,0,0,21,0,19,0,11,
  0,0,0,4,0,23,20,0,0,13,
  0,0,0,0,0,0,0,0,0,18,
  0,0,0,0,0,8,1,3,5,0),
  10,10,byrow=TRUE,   ### byrow=TRUE indicates that the matrix is being read in the same direction
  dimnames=list(rnames,cnames))
directional.model3

#fitdir3<-fitMk(tree,locomotion.mode,model=directional.model3,pi="fitzjohn")
#plot(fitdir3,show.zeros=FALSE,color=TRUE,width=TRUE)


## Directional model 4 --------------- ##

# as directional model 4, but transitions between semi-specialized states are allowed:
#TQ <--> semiaquat <--> aquat
#TQ <--> semifoss <--> foss
#TQ< --> scansorial <--> arboreal <--> gliding <--> flight
#TQ <--> TB
#Semiaquat <--> Semifoss
#Semiaquat <--> Scansorial
#Semifoss <--> Scansorial
#semifoss <--> TB

directional.model4<-matrix(c(
  0,0,0,0,0,0,10,0,0,0,
  0,0,0,0,6,16,0,0,0,0,
  0,0,0,0,14,0,0,0,0,0,
  0,0,0,0,0,0,0,12,0,0,
  0,15,7,0,0,0,0,0,0,0,
  0,9,0,0,0,0,22,24,0,17,
  2,0,0,0,0,21,0,19,0,11,
  0,0,0,4,0,23,20,0,25,13,
  0,0,0,0,0,0,0,26,0,18,
  0,0,0,0,0,8,1,3,5,0),
  10,10,byrow=TRUE,   ### byrow=TRUE indicates that the matrix is being read in the same direction
  dimnames=list(rnames,cnames))
directional.model4

fitdir4<-fitMk(tree,locomotion.mode,model=directional.model4,pi="fitzjohn")
plot(fitdir4,show.zeros=FALSE,color=TRUE,width=TRUE)




#-------------------------------------------------------------------------------------#
## b) comparing all models
#-------------------------------------------------------------------------------------#

#fitMk_list <- list(fitER=fitER,fitARD=fitARD,fitSYM=fitSYM,fitdir1=fitdir1,fitdir2=fitdir2,fitdir3=fitdir3)
#save(fitMk_list, file="fitMk_list.RData")
load("fitMk_list.RData")


fitER <- fitMk_list$fitER
fitARD <- fitMk_list$fitARD
fitSYM <- fitMk_list$fitSYM
fitdir1 <- fitMk_list$fitdir1
fitdir2 <- fitMk_list$fitdir2
fitdir3 <- fitMk_list$fitdir3


lapply(fitMk_list, function(x) logLik(x))

par(mfrow=c(3,2))
plot(fitER,show.zeros=FALSE,color=TRUE)
plot(fitARD,show.zeros=FALSE,color=TRUE,width=TRUE)
plot(fitSYM,show.zeros=FALSE,color=TRUE,width=TRUE)
plot(fitdir1,show.zeros=FALSE,color=TRUE,width=TRUE)
plot(fitdir2,show.zeros=FALSE,color=TRUE,width=TRUE)
plot(fitdir3,show.zeros=FALSE,color=TRUE,width=TRUE)

AOV<-anova(fitER,fitARD,fitSYM,fitdir1,fitdir2,fitdir3, fitdir4) # model with the best support is fitdir2
aicw(setNames(AOV$AIC, rownames(AOV)))




#-------------------------------------------------------------------------------------#
## c) make simmap trees
#-------------------------------------------------------------------------------------#

## mapping character evolution using the best evol fitMk model (fitdir2)

#simmap_Dir2 <- simmap(fitdir2) #n=100
#save(simmap_Dir2, file="simmap_Dir2.RData")


##=========================================================================================================================##
#                                                                                                                           #   
#### Pt IV) GLS model fit                                                                                                ####
#                                                                                                                           #
##=========================================================================================================================##

load("simmap_Dir2.RData")

tr <- as.matrix(residM)
trbm <- as.matrix(cbind(bodymass, limb))

#-------------------------------------------------------------------------------------#
## a) With body mass - understanding the impact of locomotion on body size evolution
#-------------------------------------------------------------------------------------#

## loading the fitted objects calculated below:
load(file="fit_trbm.RData")
fit_BM1=fit_trbm$fit_BM1; fit_OU1=fit_trbm$fit_OU1; fit_EB1=fit_trbm$fit_EB1; fit_BMM1m=fit_trbm$fit_BMM1m; fit_OUM1m=fit_trbm$fit_OUM1m

#fit_BM1 <- mvgls(trbm~1, tree=tree, model="BM", penalty="LASSO", method="LL")
#fit_OU1 <- mvgls(trbm~1, tree=tree, model="OU", penalty="LASSO", method="LL")
#fit_EB1 <- mvgls(trbm~1, tree=tree, model="EB", penalty="LASSO", method="LL")
#fit_BMM1 <- mvgls(trbm~1, tree=simmap_Dir2[[1]], model="BMM", penalty="LASSO", method="LL") # checking with one simmap tree
#fit_OUM1 <- mvgls(trbm~1, tree=simmap_Dir2[[1]], model="OUM", penalty="LASSO", method="LL") # checking with one simmap tree

GIC(fit_BM1); GIC(fit_OU1); GIC(fit_EB1); #GIC(fit_BMM1); GIC(fit_OUM1) #OUM lowest GIC

#fit_BMM1m <- lapply(simmap_Dir2[1:100], function(x) mvgls(trbm~1, tree=x, penalty="LASSO", model="BMM", method = "LL")) # now with 100 simmap trees
gic_bmm1 <- lapply(fit_BMM1m, function(x) GIC(x))
mean_gic_bmm1 <- mean(unlist(lapply(gic_bmm1, function(x) x$GIC)))
mean_ll_bmm1 <- mean(unlist(lapply(fit_BMM1m, function(x) x$logLik)))

#fit_OUM1m <- lapply(simmap_Dir2[1:100], function(x) mvgls(trbm~1, tree=x, penalty="LASSO", model="OUM", method = "LL"))
gic_oum1 <- lapply(fit_OUM1m, function(x) GIC(x))
mean_gic_oum1 <- mean(unlist(lapply(gic_oum1, function(x) x$GIC)))
mean_ll_oum1 <- mean(unlist(lapply(fit_OUM1m, function(x) x$logLik)))

#mean_gic_bmm1; mean_gic_oum1  #OUM with best fit
#fit_OUM1$coefficients


GIC(fit_BM1); GIC(fit_OU1); GIC(fit_EB1); mean_gic_bmm1; mean_gic_oum1
fit_BM1$logLik; fit_OU1$logLik; fit_EB1$logLik; mean_ll_bmm1; mean_ll_oum1

#fit_trbm <- list(fit_BM1=fit_BM1, fit_OU1=fit_OU1, fit_EB1=fit_EB1, fit_BMM1m=fit_BMM1m, fit_OUM1m=fit_OUM1m)
#save(fit_trbm, file="fit_trbm.RData")




#-------------------------------------------------------------------------------------#
## b) body size residuals
#-------------------------------------------------------------------------------------#

load(file="fit_tr.RData")
fit_BM2r=fit_tr$fit_BM2r; fit_OU2r=fit_tr$fit_OU2r; fit_EB2r=fit_tr$fit_EB2r; fit_BMM2mr=fit_tr$fit_BMM2mr; fit_OUM2mr=fit_tr$fit_OUM2mr

#fit_BM2r <- mvgls(tr~1, tree=tree, model="BM", penalty="LASSO", method="LL")
#fit_OU2r <- mvgls(tr~1, tree=tree, model="OU", penalty="LASSO", method="LL")
#fit_EB2r <- mvgls(tr~1, tree=tree, model="EB", penalty="LASSO", method="LL")
#fit_BMM2r <- mvgls(tr~1, tree=simmap_Dir2[[1]], model="BMM", penalty="LASSO", method="LL")
#fit_OUM2r <- mvgls(tr~1, tree=simmap_Dir2[[1]], model="OUM", penalty="LASSO", method="LL")
GIC(fit_BM2r); GIC(fit_OU2r); GIC(fit_EB2r); GIC(fit_BMM2r); GIC(fit_OUM2r) #OUM with best fit

#fit_BMM2mr <- lapply(simmap_Dir2[1:100], function(x) mvgls(tr~1, tree=x, penalty="LASSO", model="BMM", method = "LL"))
gic_bmm2r <- lapply(fit_BMM2mr, function(x) GIC(x))
mean_gic_bmm2r <- mean(unlist(lapply(gic_bmm2r, function(x) x$GIC)))
mean_ll_bmm2r <- mean(unlist(lapply(fit_BMM2mr, function(x) x$logLik)))

#fit_OUM2mr <- lapply(simmap_Dir2[1:100], function(x) mvgls(tr~1, tree=x, penalty="LASSO", model="OUM", method = "LL"))
gic_oum2r <- lapply(fit_OUM2mr, function(x) GIC(x))
mean_gic_oum2r <- mean(unlist(lapply(gic_oum2r, function(x) x$GIC)))
mean_ll_oum2r <- mean(unlist(lapply(fit_OUM2mr, function(x) x$logLik)))


GIC(fit_BM2r); GIC(fit_OU2r); GIC(fit_EB2r); mean_gic_bmm2r; mean_gic_oum2r
fit_BM2r$logLik; fit_OU2r$logLik; fit_EB2r$logLik; mean_ll_bmm2r; mean_ll_oum2r

#fit_tr <- list(fit_BM2r=fit_BM2r, fit_OU2r=fit_OU2r, fit_EB2r=fit_EB2r, fit_BMM2mr=fit_BMM2mr, fit_OUM2mr=fit_OUM2mr)
#save(fit_tr, file="fit_tr.RData")






##==========================================================================================================##
##                                                                                                          ##
#### Pt V) Simulation to verify if simulated models recover empirical data                               ####
##                                                                                                          ##
##==========================================================================================================##


fit_BMM1m
fit_OUM1m
fit_BMM2mr
fit_OUM2mr


##-----------------------------------------#
## I) With body mass
##-----------------------------------------#

# loading results for the following analyses 

load(file="listmultifit1.RData")

fit_OUM1m10=listmultifit1$fit_OUM1m10; 
listsimfit_OUM1=listmultifit1$listsimfit_OUM1; multifitBM1=listmultifit1$multifitBM1; multifitOU1=listmultifit1$multifitOU1; 
multifitEB1=listmultifit1$multifitEB1; multifitBMM1=listmultifit1$multifitBMM1; multifitOUM1=listmultifit1$multifitOUM1;
gic_bm1sim=listmultifit1$gic_bm1sim; gic_ou1sim=listmultifit1$gic_ou1sim; gic_eb1sim=listmultifit1$gic_eb1sim; 
gic_bmm1sim=listmultifit1$gic_bmm1sim; gic_oum1sim=listmultifit1$gic_oum1sim



#fit_OUM1m10 <- sample(fit_OUM1m, size=10)
#listsimfit_OUM1 <- list()
#for(i in 10){
#  for(j in 1:10)
#    listsimfit_OUM1[[j]] <-  lapply(1:10, function(x) fitted(fit_OUM1m10[[i]]) + mvSIM(fit_OUM1m10[[i]]$corrSt$phy, model="BM1", param=list(sigma=fit_OUM1m10[[i]]$sigma$Pinv, theta=rep(0,14)), nsim=1))
#}

#simfit_OUM1 <- do.call(c, listsimfit_OUM1)

#multifitBM1 <- lapply(1:100, function(x) mvgls(simfit_OUM1[[x]]~1, model="BM", tree=tree,  penalty="LASSO", method="LL"))
#multifitOU1 <- lapply(1:100, function(x) mvgls(simfit_OUM1[[x]]~1, model="OU", tree=tree,  penalty="LASSO", method="LL"))
#multifitEB1 <- lapply(1:100, function(x) mvgls(simfit_OUM1[[x]]~1, model="EB", tree=tree,  penalty="LASSO", method="LL"))
#multifitBMM1 <- lapply(1:100, function(x) mvgls(simfit_OUM1[[x]]~1, model="BMM", tree=simmap_Dir2[[1]], penalty="LASSO", method="LL")) 
#multifitOUM1 <- lapply(1:100, function(x) mvgls(simfit_OUM1[[x]]~1, model="OUM", tree=simmap_Dir2[[1]], penalty="LASSO", method="LL"))


#gic_bm1sim <- lapply(multifitBM1, function(x) GIC(x))
#gic_ou1sim <- lapply(multifitOU1, function(x) GIC(x))
#gic_eb1sim <- lapply(multifitEB1, function(x) GIC(x))
#gic_bmm1sim <- lapply(multifitBMM1, function(x) GIC(x))
#gic_oum1sim <- lapply(multifitOUM1, function(x) GIC(x))

mean(unlist(lapply(gic_bm1sim, function(x) x$GIC)))
mean(unlist(lapply(gic_ou1sim, function(x) x$GIC)))
mean(unlist(lapply(gic_eb1sim, function(x) x$GIC)))
mean(unlist(lapply(gic_bmm1sim, function(x) x$GIC)))
mean(unlist(lapply(gic_oum1sim, function(x) x$GIC)))


#listmultifit1 <- list(fit_OUM1m10=fit_OUM1m10, listsimfit_OUM1=listsimfit_OUM1, multifitBM1=multifitBM1, multifitOU1=multifitOU1, multifitEB1=multifitEB1, multifitBMM1=multifitBMM1, multifitOUM1=multifitOUM1,
#                      gic_bm1sim=gic_bm1sim, gic_ou1sim=gic_ou1sim, gic_eb1sim=gic_eb1sim, gic_bmm1sim=gic_bmm1sim, gic_oum1sim=gic_oum1sim)
#save(listmultifit1, file="listmultifit1.RData")

##-----------------------------------------#
## II) Residual
##-----------------------------------------#

load(file="listmultifit2r.RData")

fit_OUM2mr10=listmultifit2r$fit_OUM2mr10; listsimfit_OUM2r=listmultifit2r$listsimfit_OUM2r; multifitBM2r=listmultifit2r$multifitBM2r; 
multifitOU2r=listmultifit2r$multifitOU2r; multifitEB2r=listmultifit2r$multifitEB2r; multifitBMM2r=listmultifit2r$multifitBMM2r; multifitOUM2r=listmultifit2r$multifitOUM2r;
gic_bm2rsim=listmultifit2r$gic_bm2rsim; gic_ou2rsim=listmultifit2r$gic_ou2rsim; gic_eb2rsim=listmultifit2r$gic_eb2rsim; 
gic_bmm2rsim=listmultifit2r$gic_bmm2rsim; gic_oum2rsim=listmultifit2r$gic_oum2rsim


#fit_OUM2mr10 <- sample(fit_OUM2mr, size=10)
#listsimfit_OUM2r <- list()
#for(i in 10){
#  for(j in 1:10)
#    listsimfit_OUM2r[[j]] <-  lapply(1:10, function(x) fitted(fit_OUM2mr10 [[i]]) + mvSIM(fit_OUM2mr10 [[i]]$corrSt$phy, model="BM1", param=list(sigma=fit_OUM2mr10 [[i]]$sigma$Pinv, theta=rep(0,13)), nsim=1))
#}

#simfit_OUM2r <- do.call(c, listsimfit_OUM2r)

#multifitBM2r <- lapply(1:100, function(x) mvgls(simfit_OUM2r[[x]]~1, model="BM", tree=tree,  penalty="LASSO", method="LL"))
#multifitOU2r<- lapply(1:100, function(x) mvgls(simfit_OUM2r[[x]]~1, model="OU", tree=tree,  penalty="LASSO", method="LL"))
#multifitEB2r <- lapply(1:100, function(x) mvgls(simfit_OUM2r[[x]]~1, model="EB", tree=tree,  penalty="LASSO", method="LL"))
#multifitBMM2r <- lapply(1:100, function(x) mvgls(simfit_OUM2r[[x]]~1, model="BMM", tree=simmap_Dir2[[1]], penalty="LASSO", method="LL"))
#multifitOUM2r <- lapply(1:100, function(x) mvgls(simfit_OUM2r[[x]]~1, model="OUM", tree=simmap_Dir2[[1]], penalty="LASSO", method="LL"))


#gic_bm2rsim <- lapply(multifitBM2r, function(x) GIC(x))
#gic_ou2rsim <- lapply(multifitOU2r, function(x) GIC(x))
#gic_eb2rsim <- lapply(multifitEB2r, function(x) GIC(x))
#gic_bmm2rsim <- lapply(multifitBMM2r, function(x) GIC(x))
#gic_oum2rsim <- lapply(multifitOUM2r, function(x) GIC(x))

mean(unlist(lapply(gic_bm2rsim, function(x) x$GIC)))
mean(unlist(lapply(gic_ou2rsim, function(x) x$GIC)))
mean(unlist(lapply(gic_eb2rsim, function(x) x$GIC)))
mean(unlist(lapply(gic_bmm2rsim, function(x) x$GIC)))
mean(unlist(lapply(gic_oum2rsim, function(x) x$GIC)))



#listmultifit2r <- list(fit_OUM2mr10=fit_OUM2mr10, listsimfit_OUM2r=listsimfit_OUM2r, multifitBM2r=multifitBM2r, multifitOU2r=multifitOU2r, multifitEB2r=multifitEB2r, multifitBMM2r=multifitBMM2r, multifitOUM2r=multifitOUM2r,
#                      gic_bm2rsim=gic_bm2rsim, gic_ou2rsim=gic_ou2rsim, gic_eb2rsim=gic_eb2rsim, gic_bmm2rsim=gic_bmm2rsim, gic_oum2rsim=gic_oum2rsim)
#save(listmultifit2r, file="listmultifit2r.RData")





###==========================================================================================================##
##                                                                                                          ##
####  Plot GIC                                                                                            ####
##                                                                                                          ##
##==========================================================================================================##

library(reshape2)

gic1 <- melt(data.frame(
  bm1= as.numeric(unlist(lapply(gic_bm1sim, function(x) x$GIC))),
  ou1= as.numeric(unlist(lapply(gic_ou1sim, function(x) x$GIC))),
  eb1= as.numeric(unlist(lapply(gic_eb1sim, function(x) x$GIC))),
  bmm1= as.numeric(unlist(lapply(gic_bmm1sim, function(x) x$GIC))),
  oum1= as.numeric(unlist(lapply(gic_oum1sim, function(x) x$GIC)))
))


gic2 <- melt(data.frame(
  bm2= as.numeric(unlist(lapply(gic_bm2rsim, function(x) x$GIC))),
  ou2= as.numeric(unlist(lapply(gic_ou2rsim, function(x) x$GIC))),
  eb2= as.numeric(unlist(lapply(gic_eb2rsim, function(x) x$GIC))),
  bmm2= as.numeric(unlist(lapply(gic_bmm2rsim, function(x) x$GIC))),
  oum2= as.numeric(unlist(lapply(gic_oum2rsim, function(x) x$GIC)))
))


p.gic.raw <- ggplot(gic1, aes(x=variable, y=value)) +
  theme_classic()+
  geom_point(size=1.5, shape=21, fill="#0066CC", color="#0066CC", position=position_jitter(width=0.1, height=0.1), alpha=0.3) +
  geom_boxplot(width=0.15, color="black", fill="gray60", outlier.colour = NA,  alpha=0.3) +
  xlab(NULL)+
  ylab("GIC - simulated raw data")+
  theme(legend.position="none")+
  theme(axis.text.y = element_text(size = 10), axis.title.y = element_text(size=12),
        axis.text.x = element_text(size = 10, angle=45, vjust = 0.7))

p.gic.resid <- ggplot(gic2, aes(x=variable, y=value)) +
  theme_classic()+
  geom_point(size=1.5, shape=21, fill="#0066CC", color="#0066CC", position=position_jitter(width=0.1, height=0.1), alpha=0.3) +
  geom_boxplot(width=0.15, color="black", fill="gray60", outlier.colour = NA,  alpha=0.3) +
  xlab(NULL)+
  ylab("GIC -simulated residual data")+
  theme(legend.position="none")+
  theme(axis.text.y = element_text(size = 10), axis.title.y = element_text(size=12),
        axis.text.x = element_text(size = 10, angle=45, vjust = 0.7))


cowplot::plot_grid(p.gic.raw, p.gic.resid, nrow = 1)



##==========================================================================================================##
##                                                                                                          ##
#### Pt VI) Comparing thetas                                                                             ####
##                                                                                                          ##
##==========================================================================================================##



##-----------------------------------------#
## I) Theta raw OUM
##-----------------------------------------#


fit_OUM1m


for(i in 1:100){
  rownames(fit_OUM1m[[i]]$coeff) <- rownames(simmap_Dir2[[1]]$Q)
  colnames(fit_OUM1m[[i]]$coeff) <- colnames(trbm)
}
fit_OUM1m[[1]]$coeff
sumtheta1 <- lapply(fit_OUM1m, function(x) apply(x$coeff, 1, sum))
theta.table <- data.frame(                                                                        
  theta = as.numeric(unlist(lapply(sumtheta1, function(x) x))),
  locomotion = factor(rep(c("aquatic","arboreal","flight","fossorial","gliding","scansorial","semiaquatic","semifossorial","terrestrialbip","terrestrialquad"), 100), 
                      levels = c("aquatic","arboreal","flight","fossorial","gliding","scansorial","semiaquatic","semifossorial","terrestrialbip","terrestrialquad"))
)


boxplot(theta~locomotion, data=theta.table)

aovtheta <- aov(theta~locomotion, data=theta.table)
summary(aovtheta)
tukey <- TukeyHSD(aovtheta)$locomotion
range(tukey[,2])



## calculating mean theta values per trait


fit_OUM1m[[1]]$coeff
fit_OUM1m.thetas <- lapply(fit_OUM1m, function(x) x$coeff)
fit_OUM1m.thetas.df <- as.data.frame(do.call("rbind", fit_OUM1m.thetas))
fit_OUM1m.thetas.df$locomotion <- rep(c("aquatic","arboreal","flight","fossorial","gliding","scansorial","semiaquatic","semifossorial","terrestrialbip","terrestrialquad"), 100)

loc.factor <-  c("aquatic","arboreal","flight","fossorial","gliding","scansorial","semiaquatic","semifossorial","terrestrialbip","terrestrialquad")

list1 <- list()
means.theta.OUM1m <- list()
for(i in 1:length(loc.factor)){
list1[[i]] <- subset(fit_OUM1m.thetas.df, locomotion==loc.factor[i])
means.theta.OUM1m[[i]] <- apply(list1[[i]][,1:14], 2, mean)
}
names(means.theta.OUM1m) <- loc.factor
means.theta.OUM1m <- do.call("rbind", means.theta.OUM1m)


##-----------------------------------------##
## II) Theta resid OUM
##-----------------------------------------##

for(i in 1:100){
  rownames(fit_OUM2mr[[i]]$coefficients) <- rownames(simmap_Dir2[[1]]$Q)
  colnames(fit_OUM2mr[[i]]$coefficients) <- colnames(tr)
}

sumtheta2r <- lapply(fit_OUM2mr, function(x) apply(x$coefficients, 1, sum))
theta.table2r <- data.frame(                                                                        
  theta = as.numeric(unlist(lapply(sumtheta2r, function(x) x))),
  locomotion = factor(rep(c("aquatic","arboreal","flight","fossorial","gliding","scansorial","semiaquatic","semifossorial","terrestrialbip","terrestrialquad"), 100), 
                      levels = c("aquatic","arboreal","flight","fossorial","gliding","scansorial","semiaquatic","semifossorial","terrestrialbip","terrestrialquad"))
)


boxplot(theta~locomotion, data=theta.table2r)


aovtheta2 <- aov(theta~locomotion, data=theta.table2r)
summary(aovtheta2)
tukey2 <- TukeyHSD(aovtheta2)$locomotion
range(tukey2[,2])



## calculating mean theta values per trait

fit_OUM2mr[[1]]$coeff
fit_OUM2mr.thetas <- lapply(fit_OUM2mr, function(x) x$coeff)
fit_OUM2mr.thetas.df <- as.data.frame(do.call("rbind", fit_OUM2mr.thetas))
fit_OUM2mr.thetas.df$locomotion <- rep(c("aquatic","arboreal","flight","fossorial","gliding","scansorial","semiaquatic","semifossorial","terrestrialbip","terrestrialquad"), 100)

loc.factor <-  c("aquatic","arboreal","flight","fossorial","gliding","scansorial","semiaquatic","semifossorial","terrestrialbip","terrestrialquad")

list1 <- list()
means.theta.OUM2mr <- list()
for(i in 1:length(loc.factor)){
  list1[[i]] <- subset(fit_OUM2mr.thetas.df, locomotion==loc.factor[i])
  means.theta.OUM2mr[[i]] <- apply(list1[[i]][,1:13], 2, mean)
}
names(means.theta.OUM2mr) <- loc.factor
means.theta.OUM2mr <- do.call("rbind", means.theta.OUM2mr)

apply(means.theta.OUM2mr, 2, mean)


##==========================================================================================================##
##                                                                                                          ##
####  plot thetas                                                                                         ####
##                                                                                                          ##
##==========================================================================================================##

cols<-setNames(c("blue","forestgreen", "azure3","darkred", "black","olivedrab1", "skyblue", "red", "violet", "#FFC055"),
               c("aquatic","arboreal","flight","fossorial","gliding","scansorial","semiaquatic","semifossorial","terrestrialbip","terrestrialquad"))


##-----------------------------------------#
## I) Density theta
##-----------------------------------------#

theta.table
p.thetaraw <- ggplot(theta.table, aes(x=theta, fill=locomotion, color=locomotion))+
  theme_classic()+
  geom_density(alpha=0.5)+
  scale_fill_manual(values=cols)+
  scale_color_manual(values=cols)+
  theme(legend.position="none")+
  xlab("Sum theta - raw values")+ ylab("")+
  theme(axis.text.y = element_text(size = 10), axis.title.y = element_text(size=12),
        axis.text.x = element_text(size = 10), axis.title.x = element_text(size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=10))


theta.table2r 
p.thetaresid <- ggplot(theta.table2r, aes(x=theta, fill=locomotion, color=locomotion))+
  theme_classic()+
  geom_density(alpha=0.5)+
  scale_fill_manual(values=cols)+ 
  scale_color_manual(values=cols)+
  theme(legend.position="none")+
  xlab("Sum theta - residual values")+ ylab("")+
  theme(axis.text.y = element_text(size = 10), axis.title.y = element_text(size=12),
        axis.text.x = element_text(size = 10), axis.title.x = element_text(size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=10))


leg <- get_legend(ggplot(theta.table2r, aes(x=theta, fill=locomotion, color=locomotion))+
                    theme_classic()+
                    geom_density(alpha=0.5)+
                    scale_fill_manual(values=cols)+ 
                    scale_color_manual(values=cols)+
                    guides(col=guide_legend(ncol=1)))

ptheta1 <- plot_grid(p.thetaraw,p.thetaresid, nrow=1) 
plot_grid(ptheta1, leg, nrow =1, rel_widths = c(3,1))
ptheta2 <- plot_grid(p.thetaraw,p.thetaresid, ncol=1)
plot_grid(ptheta2, leg, nrow =1, rel_widths = c(3,1))


##-----------------------------------------#
## II) 3D Theta
##-----------------------------------------#


## raw

coeff1 <- lapply(fit_OUM1m, function(x) x$coefficients)
for(i in 1:100){
  colnames(coeff1[[i]]) <- colnames(trbm)  
}

dfcoeff1 <- do.call("rbind", coeff1)
dfcoeff1loc <- data.frame(dfcoeff1,
                          locomotion=factor(rep(c("aquatic","arboreal","flight","fossorial","gliding","scansorial","semiaquatic","semifossorial","terrestrialbip","terrestrialquad"), 100), levels = c("aquatic","arboreal","flight","fossorial","gliding","scansorial","semiaquatic","semifossorial","terrestrialbip","terrestrialquad"))
)


pca1 <- prcomp(dfcoeff1loc[1:14])
summary(pca1)
eigenv <- pca1$sd^2
scores <- as.data.frame(pca1$x)
scores$locomotion <-dfcoeff1loc$locomotion

library(rgl)
library(car)
library(plotly)
p <-plot_ly(scores, x=~PC3, y=~PC1, 
            z=~PC2, color=~locomotion, colors=cols) %>%
  add_markers(size=1,showlegend = FALSE  )%>%
layout(scene = list(
       yaxis = list(title = 'PC1 (86.4%)'), 
       zaxis = list(title = "PC2 (10.8%)"), 
       xaxis = list(title="PC3 (1.5%)"), 
       aspectmode='cube'))
print(p)



scatter3d(data=scores, y=scores$PC1, x=scores$PC2, 
          z=scores$PC3, ylab="PC1 - size (86.4%)", xlab="PC2 - hand length (10.8%)", zlab="PC3 - radius length (1.5%)" , groups= scores$locomotion,
          surface=FALSE, grid = FALSE, ellipsoid = TRUE, ellipsoid.alpha=0.2,
          surface.col = cols)

rgl.bbox( shininess=5, alpha=0.2 ) 



## residual

fit_OUM2mr

coeff2 <- lapply(fit_OUM2mr, function(x) x$coefficients)
for(i in 1:100){
  colnames(coeff2[[i]]) <- colnames(tr)  
}

dfcoeff2 <- do.call("rbind", coeff2)
dfcoeff2loc <- data.frame(dfcoeff2,
                        locomotion=factor(rep(c("aquatic","arboreal","flight","fossorial","gliding","scansorial","semiaquatic","semifossorial","terrestrialbip","terrestrialquad"), 100), levels = c("aquatic","arboreal","flight","fossorial","gliding","scansorial","semiaquatic","semifossorial","terrestrialbip","terrestrialquad"))
)

pca2 <- prcomp(dfcoeff2loc[1:13])
summary(pca2)
pca2$sd^2
scores2 <- as.data.frame(pca2$x)
scores2$locomotion <-dfcoeff2loc$locomotion



p2 <-plot_ly(scores2, x=~PC3, y=~PC1, 
            z=~PC2, color=~locomotion, colors=cols) %>%
  add_markers(size=1,showlegend = FALSE )%>%
  layout(scene = list(
    yaxis = list(title = 'PC1 (76.3%)'), 
    zaxis = list(title = "PC2 (12.5%)"), 
    xaxis = list(title="PC3 (7.0%)"), 
    aspectmode='cube'))
print(p2)
require(rgl)

snapshot3d(p2)

##==========================================================================================================##
##                                                                                                          ##
#### Pt VII) Mahalanobis and q                                                                           ####
##                                                                                                          ##
##==========================================================================================================##



##==========================================================================================================##
#### a) peaks with size ####
##==========================================================================================================##

# calculate trait optimal means per category
mat.theta <- cov(dfcoeff1loc[1:14])
means.list <- list()
for(i in 1:length(loc.factor)){
  means.list[[i]] <- colMeans(subset(dfcoeff1loc, locomotion==loc.factor[i] )[1:14])
}

names(means.list) <- loc.factor



# Calculate Mahalanobis distance between each pair of categories
# first creating possible pai combinations
pair <- combn(loc.factor, 2) 
pair.names <- vector()
for(i in 1:45){
pair.names[i] <-  paste(pair[1,i], pair[2,i], sep="_")
}

# now calculating mahalanobis distances between pairs
mah.list1 <- list()
for(i in 1:ncol(pair)){
mah.list1[[i]] <- mahalanobis(means.list[[pair[1,i]]], means.list[[pair[2,i]]], mat.theta)
}

mah.dist1 <- as.data.frame(unlist(mah.list1))
mah.dist1$pair1 <- pair.names
colnames(mah.dist1) <- c("mahalanobis", "pair")


# let's also create new factors with inverted pair, since the sequence is important for q
mah.list2 <- list()
for(i in 1:ncol(pair)){
  mah.list2[[i]] <- mahalanobis(means.list[[pair[2,i]]], means.list[[pair[1,i]]], mat.theta)
}

pair.names2 <- vector()
for(i in 1:45){
  pair.names2[i] <-  paste(pair[2,i], pair[1,i], sep="_")
}

mah.dist2 <- as.data.frame(unlist(mah.list2))
mah.dist2$pair2 <- pair.names2
colnames(mah.dist2) <- c("mahalanobis", "pair")

mah.dist <- rbind(mah.dist1, mah.dist2)


## extracting q
q.mat <- as.Qmatrix(fitdir2)
q.val <- as.vector(q.mat)

aquatic.names <- vector()
arboreal.names <- vector()
flight.names <- vector()
fossorial.names <- vector()
gliding.names <- vector()
scansorial.names <- vector()
semiaquatic.names <- vector()
semifossorial.names <- vector()
tererstrialbip.names <- vector()
tererstrialquad.names <- vector()

for(i in 1:10){
  aquatic.names[i] <-  paste(row.names(q.mat)[i], 'aquatic', sep="_")
  arboreal.names[i] <-  paste(row.names(q.mat)[i],'arboreal', sep="_")
  flight.names[i] <-  paste(row.names(q.mat)[i],'flight', sep="_")
  fossorial.names[i] <-  paste(row.names(q.mat)[i],'fossorial', sep="_")
  gliding.names[i] <-  paste(row.names(q.mat)[i],'gliding', sep="_")
  scansorial.names[i] <-  paste( row.names(q.mat)[i], 'scansorial', sep="_")
  semiaquatic.names[i] <-  paste(row.names(q.mat)[i], 'semiaquatic', sep="_")
  semifossorial.names[i] <-  paste(row.names(q.mat)[i], 'semifossorial', sep="_")
  tererstrialbip.names[i] <-  paste(row.names(q.mat)[i], 'terrestrialbip',sep="_")
  tererstrialquad.names[i] <-  paste(row.names(q.mat)[i],'terrestrialquad', sep="_")
  
  }

names.q <- c(aquatic.names,arboreal.names,flight.names,fossorial.names,gliding.names,scansorial.names,semiaquatic.names, semifossorial.names,tererstrialbip.names,tererstrialquad.names)
names(q.val) <- names.q


match1 <- data.frame(q=q.val[mah.dist$pair])
match.mq1 <- cbind(match1, mah.dist)
lm1 <- lm(mahalanobis~q, data=match.mq1)

plot(mahalanobis~q, data=match.mq1)
abline(lm1)

plot_mq1 <- ggplot(match.mq1, aes(x=q, y=mahalanobis, text=pair)) + 
  geom_point(size=2, alpha=0.5) +
  theme_minimal()+
  xlab("Transition rates (q)")+
  ylab("Mahalanobis distance between peaks")+
  guides(col=guide_legend(ncol=1))+
  theme(axis.text.y = element_text(size = 10), axis.title.y = element_text(size=12),
        axis.text.x = element_text(size = 10), axis.title.x = element_text(size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=10))

plotly::ggplotly(plot_mq1, tooltip=c("text", "x", "y"))



library(ggpubr)


ggplot(match.mq1, aes(x=q, y=mahalanobis)) + 
  geom_smooth(method = "lm", se = TRUE) +
  geom_point(size=2, alpha=0.5) +
  theme_minimal()+
  xlab("Transition rates (q)")+
  ylab("Mahalanobis distance between peaks")+
  theme(axis.text.y = element_text(size = 10), axis.title.y = element_text(size=12),
        axis.text.x = element_text(size = 10), axis.title.x = element_text(size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=10))



##-----------------------------------------#
## separating transitions that q=0 and q>0
##-----------------------------------------#


match.mq1
names.null.q <- rownames(subset(match.mq1, q==0))
null.q <- rep("null", 75)
names(null.q) <- names.null.q

names.allowed.q <- rownames(subset(match.mq1, q>0))
length(names.allowed.q)
allowed.q <- rep("allowed", 15)
names(allowed.q) <- names.allowed.q

possible.q <- as.data.frame(c(null.q, allowed.q))
match.mq1$possible.q <- possible.q[match.mq1$pair,]
match.mq1$possible.q <- factor(match.mq1$possible.q, levels = c("null", "allowed"))


library(ggpubr)
library(rstatix)
stat.test <- match.mq1 %>%
   t_test(mahalanobis ~ possible.q) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()
stat.test 
stat.test <- stat.test %>% add_xy_position(x = "possible.q")

a <- ggplot(match.mq1, aes(x=possible.q, y=mahalanobis)) + 
  geom_violin(fill="#E2E2E2",colour="#E2E2E2")+
  geom_boxplot(width=0.1, alpha=0.3, fill="slategrey")+
  geom_jitter(color="black", size=1.2, alpha=0.2, width=0.05) +
  theme_minimal()+
  xlab("")+
  ylab("")+
  theme(axis.text.y = element_text(size = 10), axis.title.y = element_text(size=11),
        axis.text.x = element_text(size = 10), axis.title.x = element_text(size=11),
        legend.title = element_text(size=11, face="bold"))+
  scale_x_discrete(labels=c("q=0","q>0")) +
  stat_pvalue_manual(stat.test)



allowed.q.df <- subset(match.mq1, q>0)
lm2 <- summary(lm(mahalanobis~q, data=allowed.q.df ))
p <- round(lm2$coefficients[2,4], digits=3)
corr <- round(cor(allowed.q.df$mahalanobis, allowed.q.df$q, method = 'pearson'), digits=2)

b <- ggplot(allowed.q.df , aes(x=q, y=mahalanobis)) + 
  geom_smooth(method = "lm", formula=y~x, se = TRUE, fill='#E2E2E2', alpha=1, color="slategrey") +
  geom_point(size=3.5) +
  theme_minimal()+ xlab("")+
  ylab("")+
  theme(axis.text.y = element_text(size = 10), axis.title.y = element_text(size=11),
        axis.text.x = element_text(size = 10), axis.title.x = element_text(size=11),
       )+
  geom_label(x=0.04, y=18, label="p-value = 0.013
             r = - 0.63   ", 
             label.size = 0.1, fill="white", color="black")




plot.mq <- plot_grid(a,b, ncol=1)

#create common x and y labels

library(grid)
library(gridExtra)

y.lab <- textGrob("Mahalanobis distance (θ)",  gp=gpar(col="black", fontsize=11), rot=90)
x.lab <- textGrob("Transition rates (q)",  gp=gpar(col="black", fontsize=11))
grid.arrange(arrangeGrob(plot.mq, left = y.lab, bottom = x.lab))
  



###==========================================================================================================##
#### b) peaks without size ####
##==========================================================================================================##


# calculate trait optimal means per category
mat.theta.r <- cov(dfcoeff2loc[1:13])
means.list.r <- list()
for(i in 1:length(loc.factor)){
  means.list.r[[i]] <- colMeans(subset(dfcoeff2loc, locomotion==loc.factor[i] )[1:13])
}

names(means.list.r) <- loc.factor



# Calculate Mahalanobis distance between each pair of categories
# first creating possible pai combinations
pair <- combn(loc.factor, 2) 
pair.names <- vector()
for(i in 1:45){
  pair.names[i] <-  paste(pair[1,i], pair[2,i], sep="_")
}

# now calculating mahalanobis distances between pairs
mah.list1.r <- list()
for(i in 1:ncol(pair)){
  mah.list1.r[[i]] <- mahalanobis(means.list.r[[pair[1,i]]], means.list.r[[pair[2,i]]], mat.theta.r)
}

mah.dist1.r <- as.data.frame(unlist(mah.list1.r))
mah.dist1.r$pair1 <- pair.names
colnames(mah.dist1.r) <- c("mahalanobis", "pair")


# let's also create new factors with inverted pair, since the sequence is important for q
mah.list2.r <- list()
for(i in 1:ncol(pair)){
  mah.list2.r[[i]] <- mahalanobis(means.list.r[[pair[2,i]]], means.list.r[[pair[1,i]]], mat.theta.r)
}

pair.names2 <- vector()
for(i in 1:45){
  pair.names2[i] <-  paste(pair[2,i], pair[1,i], sep="_")
}

mah.dist2.r <- as.data.frame(unlist(mah.list2.r))
mah.dist2.r$pair2 <- pair.names2
colnames(mah.dist2.r) <- c("mahalanobis", "pair")

mah.dist.r <- rbind(mah.dist1.r, mah.dist2.r)


## extracting q
q.mat <- as.Qmatrix(fitdir2)
q.val <- as.vector(q.mat)

aquatic.names <- vector()
arboreal.names <- vector()
flight.names <- vector()
fossorial.names <- vector()
gliding.names <- vector()
scansorial.names <- vector()
semiaquatic.names <- vector()
semifossorial.names <- vector()
tererstrialbip.names <- vector()
tererstrialquad.names <- vector()

for(i in 1:10){
  aquatic.names[i] <-  paste(row.names(q.mat)[i], 'aquatic', sep="_")
  arboreal.names[i] <-  paste(row.names(q.mat)[i],'arboreal', sep="_")
  flight.names[i] <-  paste(row.names(q.mat)[i],'flight', sep="_")
  fossorial.names[i] <-  paste(row.names(q.mat)[i],'fossorial', sep="_")
  gliding.names[i] <-  paste(row.names(q.mat)[i],'gliding', sep="_")
  scansorial.names[i] <-  paste( row.names(q.mat)[i], 'scansorial', sep="_")
  semiaquatic.names[i] <-  paste(row.names(q.mat)[i], 'semiaquatic', sep="_")
  semifossorial.names[i] <-  paste(row.names(q.mat)[i], 'semifossorial', sep="_")
  tererstrialbip.names[i] <-  paste(row.names(q.mat)[i], 'terrestrialbip',sep="_")
  tererstrialquad.names[i] <-  paste(row.names(q.mat)[i],'terrestrialquad', sep="_")
  
}

names.q <- c(aquatic.names,arboreal.names,flight.names,fossorial.names,gliding.names,scansorial.names,semiaquatic.names, semifossorial.names,tererstrialbip.names,tererstrialquad.names)
names(q.val) <- names.q


match1.r <- data.frame(q=q.val[mah.dist.r$pair])
match.mq1.r <- cbind(match1.r, mah.dist.r)
lm1.r <- lm(mahalanobis~q, data=match.mq1.r)
summary(lm1.r)

plot(mahalanobis~q, data=match.mq1.r)
abline(lm1.r)

plot_mq1.r <- ggplot(match.mq1.r, aes(x=q, y=mahalanobis, text=pair)) + 
  geom_point(size=2, alpha=0.5) +
  theme_minimal()+
  xlab("Transition rates (q)")+
  ylab("Mahalanobis distance between peaks")+
  guides(col=guide_legend(ncol=1))+
  theme(axis.text.y = element_text(size = 10), axis.title.y = element_text(size=12),
        axis.text.x = element_text(size = 10), axis.title.x = element_text(size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=10))

plotly::ggplotly(plot_mq1.r, tooltip=c("text", "x", "y"))



library(ggpubr)


ggplot(match.mq1.r, aes(x=q, y=mahalanobis)) + 
  geom_smooth(method = "lm", se = TRUE) +
  geom_point(size=2, alpha=0.5) +
  theme_minimal()+
  xlab("Transition rates (q)")+
  ylab("Mahalanobis distance between peaks")+
  theme(axis.text.y = element_text(size = 10), axis.title.y = element_text(size=12),
        axis.text.x = element_text(size = 10), axis.title.x = element_text(size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=10))



##-----------------------------------------#
## separating transitions that q=0 and q>0
##-----------------------------------------#


match.mq1.r
names.null.q <- rownames(subset(match.mq1.r, q==0))
null.q <- rep("null", 75)
names(null.q) <- names.null.q

names.allowed.q <- rownames(subset(match.mq1.r, q>0))
length(names.allowed.q)
allowed.q <- rep("allowed", 15)
names(allowed.q) <- names.allowed.q

possible.q <- as.data.frame(c(null.q, allowed.q))
match.mq1.r$possible.q <- possible.q[match.mq1.r$pair,]
match.mq1.r$possible.q <- factor(match.mq1.r$possible.q, levels = c("null", "allowed"))


library(ggpubr)
library(rstatix)
stat.test.r <- match.mq1.r %>%
  t_test(mahalanobis ~ possible.q) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()
stat.test.r 
stat.test.r <- stat.test.r %>% add_xy_position(x = "possible.q")

a.r <- ggplot(match.mq1.r, aes(x=possible.q, y=mahalanobis)) + 
  geom_violin(fill="#E2E2E2",colour="#E2E2E2")+
  geom_boxplot(width=0.1, alpha=0.3, fill="slategrey")+
  geom_jitter(color="black", size=1.2, alpha=0.2, width=0.05) +
  theme_minimal()+
  xlab("")+
  ylab("")+
  theme(axis.text.y = element_text(size = 10), axis.title.y = element_text(size=11),
        axis.text.x = element_text(size = 10), axis.title.x = element_text(size=11),
        legend.title = element_text(size=11, face="bold"))+
  scale_x_discrete(labels=c("q=0","q>0")) +
  stat_pvalue_manual(stat.test.r)



allowed.q.df.r <- subset(match.mq1.r, q>0)
lm2.r <- summary(lm(mahalanobis~q, data=allowed.q.df.r ))
p.r <- round(lm2.r$coefficients[2,4], digits=3)
corr.r <- round(cor(allowed.q.df.r$mahalanobis, allowed.q.df.r$q, method = 'pearson'), digits=2)

b.r <- ggplot(allowed.q.df.r , aes(x=q, y=mahalanobis)) + 
  geom_smooth(method = "lm", formula=y~x, se = TRUE, fill='#E2E2E2', alpha=1, color="slategrey") +
  geom_point(size=3.5) +
  theme_minimal()+ xlab("")+
  ylab("")+
  theme(axis.text.y = element_text(size = 10), axis.title.y = element_text(size=11),
        axis.text.x = element_text(size = 10), axis.title.x = element_text(size=11),
  )+
  geom_label(x=0.04, y=18, label="p-value = 0.012
             r = - 0.63   ", 
             label.size = 0.1, fill="white", color="black")




plot.mq.r <- plot_grid(a.r,b.r, ncol=1)

#create common x and y labels

library(grid)
library(gridExtra)

y.lab <- textGrob("Mahalanobis distance (θ)",  gp=gpar(col="black", fontsize=11), rot=90)
x.lab <- textGrob("Transition rates (q)",  gp=gpar(col="black", fontsize=11))
grid.arrange(arrangeGrob(plot.mq.r, left = y.lab, bottom = x.lab))




##-----------------------------------------#
## best fit Mk  transitions
##-----------------------------------------#

bestfit <- fitdir2
cols.red<-setNames(c("blue","forestgreen", "azure3","darkred", "black","olivedrab1", "skyblue", "red", "violet", "#FFC055"),
               c("Aq","Arb","Fl","Fo","Gl","Sc","Saq","SFo","TQ","TB"))

cols.lab <- c("white","white", "black","white","white","black", "black", "white", "white", "black")

bestfit$states <- c("Aq","Arb","Fl","Fo","Gl","Sc","SAq","SFo","TB","TQ")

xy <- plot(bestfit, signif=3, show.zeros=FALSE, width=TRUE, spacer=0.15) # cant' figure out how to increase text ; if text=NULL, it disappears
invisible(mapply(draw.circle,xy$x,xy$y,col=cols.red,border = NA,
                 MoreArgs=list(radius=0.13)))
text(xy$x,xy$y,xy$states,cex=1,col=cols.lab,
     font=2)




#####################

# repeating plot with abbreviations
abbrev.transitions <- c("Arb-Gl", "Arb-Sc", "Sc-TQ", "SAq-TQ", "SFo-TQ", "TB-TQ","SAq-Aq", "Gl-Arb", "Sc-Arb", "Gl-Fl", "SFo-Fo", "TQ-Sc","TQ-SAq", "TQ-SFo", "TQ-TB")
allowed.q.df$labels <- abbrev.transitions

lm2 <- summary(lm(mahalanobis~q, data=allowed.q.df ))
p <- round(lm2$coefficients[2,4], digits=3)
corr <- round(cor(allowed.q.df$mahalanobis, allowed.q.df$q, method = 'pearson'), digits=2)

d <- ggplot(allowed.q.df , aes(x=q, y=mahalanobis)) + 
  geom_smooth(method = "lm", formula=y~x, se = TRUE, fill='wheat', alpha=1, color="slategrey") +
  geom_point(size=3, alpha=0.4) +
  theme_minimal()+ xlab("")+
  ylab("")+
  theme(axis.text.y = element_text(size = 12), axis.title.y = element_text(size=14),
        axis.text.x = element_text(size = 12), axis.title.x = element_text(size=14),
        legend.title = element_text(size=12, face="bold"))+
  geom_label(x=0.043, y=18.5, label="p-value = 0.013
             r = -0.63           ", 
             label.size = 0.35, fill="white", color="black")+
  geom_text(label=allowed.q.df$labels, check_overlap = T, nudge_y = -0.35, size=2.7)




plot.mq2 <- plot_grid(a,d, ncol=1)

#create common x and y labels

library(grid)
library(gridExtra)

y.lab <- textGrob("Mahalanobis distance (θ)",  gp=gpar(col="black", fontsize=14), rot=90)
x.lab <- textGrob("Transition rates (q)",  gp=gpar(col="black", fontsize=14))
grid.arrange(arrangeGrob(plot.mq2, left = y.lab, bottom = x.lab))





##=========================================================================================================================##
#                                                                                                                           #   
#### extra                                                                                                               ####
#                                                                                                                           #
##=========================================================================================================================##


##==========================================================================================================##
#### EUCLIDEAN DISTANCES
##==========================================================================================================##

means.list <- list()
for(i in 1:length(loc.factor)){
  means.list[[i]] <- colMeans(subset(dfcoeff1loc, locomotion==loc.factor[i] )[1:14])
}

names(means.list) <- loc.factor






# Calculate Mahalanobis distance between each pair of categories
# first creating possible pai combinations
pair <- combn(loc.factor, 2) 
pair.names <- vector()
for(i in 1:45){
  pair.names[i] <-  paste(pair[1,i], pair[2,i], sep="_")
}

# now calculating euclidean distances
eucl.list1 <- list()
for(i in 1:ncol(pair)){
  eucl.list1[[i]] <-dist(rbind(means.list[[pair[1,i]]], means.list[[pair[2,i]]]), method="euclidean")
}

eucl.dist1 <- as.data.frame(unlist(eucl.list1))
eucl.dist1$pair1 <- pair.names
colnames(eucl.dist1) <- c("euclidean", "pair")


# let's also create new factors with inverted pair, since the sequence is important for q
eucl.list2 <- list()
for(i in 1:ncol(pair)){
  eucl.list2[[i]] <- dist(rbind(means.list[[pair[2,i]]], means.list[[pair[1,i]]]), method="euclidean")
}

pair.names2 <- vector()
for(i in 1:45){
  pair.names2[i] <-  paste(pair[2,i], pair[1,i], sep="_")
}

eucl.dist2 <- as.data.frame(unlist(eucl.list2))
eucl.dist2$pair2 <- pair.names2
colnames(eucl.dist2) <- c("euclidean", "pair")

eucl.dist <- rbind(eucl.dist1, eucl.dist2)


## extracting q
q.mat <- as.Qmatrix(fitdir2)
q.val <- as.vector(q.mat)

aquatic.names <- vector()
arboreal.names <- vector()
flight.names <- vector()
fossorial.names <- vector()
gliding.names <- vector()
scansorial.names <- vector()
semiaquatic.names <- vector()
semifossorial.names <- vector()
tererstrialbip.names <- vector()
tererstrialquad.names <- vector()

for(i in 1:10){
  aquatic.names[i] <-  paste(row.names(q.mat)[i], 'aquatic', sep="_")
  arboreal.names[i] <-  paste(row.names(q.mat)[i],'arboreal', sep="_")
  flight.names[i] <-  paste(row.names(q.mat)[i],'flight', sep="_")
  fossorial.names[i] <-  paste(row.names(q.mat)[i],'fossorial', sep="_")
  gliding.names[i] <-  paste(row.names(q.mat)[i],'gliding', sep="_")
  scansorial.names[i] <-  paste( row.names(q.mat)[i], 'scansorial', sep="_")
  semiaquatic.names[i] <-  paste(row.names(q.mat)[i], 'semiaquatic', sep="_")
  semifossorial.names[i] <-  paste(row.names(q.mat)[i], 'semifossorial', sep="_")
  tererstrialbip.names[i] <-  paste(row.names(q.mat)[i], 'terrestrialbip',sep="_")
  tererstrialquad.names[i] <-  paste(row.names(q.mat)[i],'terrestrialquad', sep="_")
  
}

names.q <- c(aquatic.names,arboreal.names,flight.names,fossorial.names,gliding.names,scansorial.names,semiaquatic.names, semifossorial.names,tererstrialbip.names,tererstrialquad.names)
names(q.val) <- names.q


match1 <- data.frame(q=q.val[eucl.dist$pair])
match.mq1 <- cbind(match1, eucl.dist)
lm1 <- lm(euclidean~q, data=match.mq1)
summary(lm1)
plot(euclidean~q, data=match.mq1)
abline(lm1)

plot_mq1 <- ggplot(match.mq1, aes(x=q, y=euclidean, text=pair)) + 
  geom_point(size=2, alpha=0.5) +
  theme_minimal()+
  xlab("Transition rates (q)")+
  ylab("Euclidean distance between peaks")+
  guides(col=guide_legend(ncol=1))+
  theme(axis.text.y = element_text(size = 10), axis.title.y = element_text(size=12),
        axis.text.x = element_text(size = 10), axis.title.x = element_text(size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=10))

plotly::ggplotly(plot_mq1, tooltip=c("text", "x", "y"))



library(ggpubr)


ggplot(match.mq1, aes(x=q, y=euclidean)) + 
  geom_smooth(method = "lm", se = TRUE) +
  geom_point(size=2, alpha=0.5) +
  theme_minimal()+
  xlab("Transition rates (q)")+
  ylab("euclalanobis distance between peaks")+
  theme(axis.text.y = element_text(size = 10), axis.title.y = element_text(size=12),
        axis.text.x = element_text(size = 10), axis.title.x = element_text(size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=10))



##-----------------------------------------#
## separating transitions that q=0 and q>0
##-----------------------------------------#


match.mq1
names.null.q <- rownames(subset(match.mq1, q==0))
null.q <- rep("null", 75)
names(null.q) <- names.null.q

names.allowed.q <- rownames(subset(match.mq1, q>0))
length(names.allowed.q)
allowed.q <- rep("allowed", 15)
names(allowed.q) <- names.allowed.q

possible.q <- as.data.frame(c(null.q, allowed.q))
match.mq1$possible.q <- possible.q[match.mq1$pair,]
match.mq1$possible.q <- factor(match.mq1$possible.q, levels = c("null", "allowed"))


library(ggpubr)
library(rstatix)
stat.test <- match.mq1 %>%
  t_test(euclidean ~ possible.q) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()
stat.test 
stat.test <- stat.test %>% add_xy_position(x = "possible.q")

a <- ggplot(match.mq1, aes(x=possible.q, y=euclidean)) + 
  geom_violin(fill="#E2E2E2",colour="#E2E2E2")+
  geom_boxplot(width=0.1, alpha=0.3, fill="slategrey")+
  geom_jitter(color="black", size=1.2, alpha=0.2, width=0.05) +
  theme_minimal()+
  xlab("")+
  ylab("")+
  theme(axis.text.y = element_text(size = 10), axis.title.y = element_text(size=11),
        axis.text.x = element_text(size = 10), axis.title.x = element_text(size=11),
        legend.title = element_text(size=11, face="bold"))+
  scale_x_discrete(labels=c("q=0","q>0")) +
  stat_pvalue_manual(stat.test)



allowed.q.df <- subset(match.mq1, q>0)
lm2 <- summary(lm(euclidean~q, data=allowed.q.df ))
p <- round(lm2$coefficients[2,4], digits=3)
corr <- round(cor(allowed.q.df$euclidean, allowed.q.df$q, method = 'pearson'), digits=2)

b <- ggplot(allowed.q.df , aes(x=q, y=euclidean)) + 
  geom_smooth(method = "lm", formula=y~x, se = TRUE, fill='#E2E2E2', alpha=1, color="slategrey") +
  geom_point(size=3.5) +
  theme_minimal()+ xlab("")+
  ylab("")+
  theme(axis.text.y = element_text(size = 10), axis.title.y = element_text(size=11),
        axis.text.x = element_text(size = 10), axis.title.x = element_text(size=11),
  )+
  geom_label(x=0.04, y=18, label="p-value = 0.013
             r = - 0.63   ", 
             label.size = 0.1, fill="white", color="black")




plot.mq <- plot_grid(a,b, ncol=1)

#create common x and y labels

library(grid)
library(gridExtra)

y.lab <- textGrob("euclidean distance (θ)",  gp=gpar(col="black", fontsize=11), rot=90)
x.lab <- textGrob("Transition rates (q)",  gp=gpar(col="black", fontsize=11))
grid.arrange(arrangeGrob(plot.mq, left = y.lab, bottom = x.lab))




##==========================================================================================================##
#### average distance from optima
##==========================================================================================================##


dist.optima <- function(data, mean){
  value <- vector()
  for(i in 1:nrow(data)){
    value[i] <- dist(rbind(data[i,], mean), method="euclidean")
  }
  print(value)
}

## raw

rownames(trbm)==rownames(locomotion)
df.loc <- cbind(trbm, locomotion)

list.sp.loc <- list()
for(i in 1:10){
  list.sp.loc[[i]] <- subset(df.loc, locomotion==loc.factor[i])[,1:14]
  
}

names(list.sp.loc) <- loc.factor


optima.loc <- list()
for (i in 1:10){
optima.loc[[i]] <- (dist.optima(list.sp.loc[[i]], means.list[[i]]))
}

names(optima.loc) <- loc.factor
df.dist.optima <- melt(optima.loc)

p1 <- ggplot(df.dist.optima, aes(y=value, x=L1, fill=L1))+
  geom_boxplot(fill=cols)


### residual

df.loc2 <- cbind(tr, locomotion)

list.sp.loc2 <- list()
for(i in 1:10){
  list.sp.loc2[[i]] <- subset(df.loc2, locomotion==loc.factor[i])[,1:13]
  
}

names(list.sp.loc2) <- loc.factor


optima.loc2 <- list()
for (i in 1:10){
  optima.loc2[[i]] <- (dist.optima(list.sp.loc2[[i]], means.list.r[[i]]))
}

names(optima.loc2) <- loc.factor
df.dist.optima2 <- melt(optima.loc2)

p2 <- ggplot(df.dist.optima2, aes(y=value, x=L1, fill=L1))+
  geom_boxplot(fill=cols)

plot_grid(p1,p2)

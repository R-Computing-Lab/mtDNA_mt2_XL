###code to run first trial of simulation


#source("~/R-Project/BalancedPed/evenInsert.R")
#source("~/R-Project/BalancedPed/famSizeCal.R")
#source("~/R-Project/BalancedPed/SimPed.R")
source("~/R-Project/BalancedPed/HelperFunctions.R")


library(OpenMx)
load("~/R-Project/BalancedPed/FixedPed.RData")
sampleFam <- ped2
Addmat <- as.matrix(ped2add(sampleFam, verbose = TRUE))
Nucmat <- ped2cn(sampleFam)
Extmat <- ped2ce(sampleFam)
Mtdmat <- ped2mt_v3(sampleFam)
Envmat <- diag(1,nrow = nrow(Addmat))
dimnames(Envmat) <- dimnames(Addmat)
Amimat <- Addmat*Mtdmat
Dmgmat <- Addmat*Addmat

# Trial-Combination 11
ad2 <- c( .6, .6, .6, .4, .4, .4, .2, .2, .2, .4, .4, .4)
dd2 <- c( .0, .0, .0, .0, .0, .0, .0, .0, .0, .0, .1, .0)
cn2 <- c( .1, .1, .1, .1, .1, .1, .1, .1, .1, .2, .1, .1)
ce2 <- c( .1, .1, .1, .1, .1, .1, .1, .1, .1, .1, .1, .1)
mt2 <- c( .1,.05,.01, .1,.05,.01, .1,.05,.01,.05,.05, .1)
am2 <- c(.05,.05,.05,.05,.05,.05,.05,.05,.05,.05,.05, .1)
ee2 <- c(.05, .1,.19,.25,.30,.39,.45,.50,.59, .2, .2, .2)


## generate data
library(mvtnorm)

comb <- 1

sumCov <- ad2[comb]*Addmat + dd2[comb]*Addmat*Addmat + cn2[comb]*Nucmat + ce2[comb]*Extmat + mt2[comb]*Mtdmat + am2[comb]*Addmat*Mtdmat + ee2[comb]*Envmat
set.seed(14271)
numfam <- round(10000/nrow(Addmat))
dat <- rmvnorm(numfam, sigma = sumCov)

#setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# fit ML model

totalVar <- 1
totalMea <- 0

ObjectsKeep <- as.character(ls())

## the full model
Model1 <- mxModel(
      "ModelOne",
      mxMatrix(type = "Full", nrow = 1, ncol = 1, free = TRUE, values = ad2[comb]*totalVar, labels = "vad", name = "Vad", lbound = 1e-10),
#      mxMatrix(type = "Full", nrow = 1, ncol = 1, free = TRUE, values = dd2[comb]*totalVar, labels = "vdd", name = "Vdd", lbound = 1e-10),
      mxMatrix(type = "Full", nrow = 1, ncol = 1, free = TRUE, values = cn2[comb]*totalVar, labels = "vcn", name = "Vcn", lbound = 1e-10),
      mxMatrix(type = "Full", nrow = 1, ncol = 1, free = TRUE, values = ce2[comb]*totalVar, labels = "vce", name = "Vce", lbound = 1e-10),
      mxMatrix(type = "Full", nrow = 1, ncol = 1, free = TRUE, values = mt2[comb]*totalVar, labels = "vmt", name = "Vmt", lbound = 1e-10),
      mxMatrix(type = "Full", nrow = 1, ncol = 1, free = TRUE, values = am2[comb]*totalVar, labels = "vam", name = "Vam", lbound = 1e-10),
      mxMatrix(type = "Full", nrow = 1, ncol = 1, free = TRUE, values = ee2[comb]*totalVar, labels = "ver", name = "Ver", lbound = 1e-10)
)

ll <- list()
for (i in 1:numfam){
      ll[[i]] <- dat[i,]
}

modList <- list()
modNames <- paste0("fam", 1:numfam)

for(afam in 1:numfam){
      ytemp <- paste('S', rownames (Addmat))
      fsize <- nrow(Addmat)
      modList[[afam]] <- mxModel(name=modNames[afam],
                                 mxMatrix("Iden", nrow=fsize, ncol=fsize, name="I"), 
                                 mxMatrix("Unit", nrow=fsize, ncol=fsize, name='U'),
                                 mxMatrix("Symm", nrow=fsize, ncol=fsize, values=Addmat, name="A"), 
##                                 mxMatrix("Symm", nrow=fsize, ncol=fsize, values=Dmgmat, name="D"), 
                                 mxMatrix("Symm", nrow=fsize, ncol=fsize, values=Nucmat, name="Cn"), 
                                 mxMatrix("Symm", nrow=fsize, ncol=fsize, values=Extmat, name="Ce"), 
                                 mxMatrix("Symm", nrow=fsize, ncol=fsize, values=Amimat, name="Am"), 
                                 mxMatrix("Symm", nrow=fsize, ncol=fsize, values=Mtdmat, name="Mt"),
                                 mxData(observed = matrix(ll[[afam]], nrow=1, dimnames=list(NULL, ytemp)), type="raw", sort=FALSE),
                                 mxMatrix('Full', nrow=1, ncol=fsize, name='M', free=TRUE, labels='meanLI',
                                          dimnames=list(NULL, ytemp)),
                                 mxAlgebra ((A %x% ModelOne.Vad) 
##                                            + (D %x% ModelOne.Vdd) 
                                            + (Cn %x% ModelOne.Vcn) 
                                            + (U %x% ModelOne.Vce) 
                                            + (Mt %x% ModelOne.Vmt) 
                                            + (Am %x% ModelOne.Vam) 
                                            + (I %x% ModelOne.Ver), 
                                            name="V", dimnames=list(ytemp, ytemp)),
                                 mxExpectationNormal(covariance='V', means='M'), 
                                 mxFitFunctionML()
      )
}
container <- mxModel('Model2', Model1, modList, mxFitFunctionMultigroup(modNames))
container <- mxOption(container, 'Checkpoint Units', 'minutes')
container <- mxOption(container, 'Checkpoint Count', 1)
containerRun <- mxRun(container, intervals=FALSE, checkpoint=TRUE) 

smr1 <- summary(containerRun)

save.image(file = "model1.RData")

rm(list = setdiff(ls(), c(ObjectsKeep, "ObjectsKeep", "smr1")))

##### Model excluding mt and am

Model4 <- mxModel(
      "ModelThree",
      mxMatrix(type = "Full", nrow = 1, ncol = 1, free = TRUE, values = ad2[comb]*totalVar, labels = "vad", name = "Vad", lbound = 1e-10),
#      mxMatrix(type = "Full", nrow = 1, ncol = 1, free = TRUE, values = dd2[comb]*totalVar, labels = "vdd", name = "Vdd", lbound = 1e-10),
      mxMatrix(type = "Full", nrow = 1, ncol = 1, free = TRUE, values = cn2[comb]*totalVar, labels = "vcn", name = "Vcn", lbound = 1e-10),
      mxMatrix(type = "Full", nrow = 1, ncol = 1, free = TRUE, values = ce2[comb]*totalVar, labels = "vce", name = "Vce", lbound = 1e-10),
      #mxMatrix(type = "Full", nrow = 1, ncol = 1, free = TRUE, values = mt2[comb]*totalVar, labels = "vmt", name = "Vmt", lbound = 1e-10),
      #mxMatrix(type = "Full", nrow = 1, ncol = 1, free = TRUE, values = am2[comb]*totalVar, labels = "vam", name = "Vam", lbound = 1e-10),
      mxMatrix(type = "Full", nrow = 1, ncol = 1, free = TRUE, values = ee2[comb]*totalVar, labels = "ver", name = "Ver", lbound = 1e-10)
)

ll2 <- list()
for (i in 1:numfam){
      ll2[[i]] <- dat[i,]
}

modList2 <- list()
modNames2 <- paste0("fam", 1:numfam)

for(afam2 in 1:numfam){
      ytemp2 <- paste('S', rownames (Addmat))
      fsize2 <- nrow(Addmat)
      modList2[[afam2]] <- mxModel(name=modNames2[afam2],
                                   mxMatrix("Iden", nrow=fsize2, ncol=fsize2, name="I"), 
                                   mxMatrix("Unit", nrow=fsize2, ncol=fsize2, name='U'),
                                   mxMatrix("Symm", nrow=fsize2, ncol=fsize2, values=Addmat, name="A"), 
#                                   mxMatrix("Symm", nrow=fsize2, ncol=fsize2, values=Dmgmat, name="D"), 
                                   mxMatrix("Symm", nrow=fsize2, ncol=fsize2, values=Nucmat, name="Cn"), 
                                   mxMatrix("Symm", nrow=fsize2, ncol=fsize2, values=Extmat, name="Ce"), 
                                   #mxMatrix("Symm", nrow=fsize2, ncol=fsize2, values=Amimat, name="Am"), 
                                   #mxMatrix("Symm", nrow=fsize2, ncol=fsize2, values=Mtdmat, name="Mt"),
                                   mxData(observed = matrix(ll2[[afam2]], nrow=1, dimnames=list(NULL, ytemp2)), type="raw", sort=FALSE),
                                   mxMatrix('Full', nrow=1, ncol=fsize2, name='M', free=TRUE, labels='meanLI',
                                            dimnames=list(NULL, ytemp2)),
                                   mxAlgebra ((A %x% ModelThree.Vad) 
#                                              + (D %x% ModelThree.Vdd) 
                                              + (Cn %x% ModelThree.Vcn) 
                                              + (U %x% ModelThree.Vce) 
                                              #+ (Mt %x% ModelThree.Vmt) 
                                              #+ (Am %x% ModelThree.Vam) 
                                              + (I %x% ModelThree.Ver), 
                                              name="V", dimnames=list(ytemp2, ytemp2)),
                                   mxExpectationNormal(covariance='V', means='M'), 
                                   mxFitFunctionML()
      )
}
container2 <- mxModel('Model4', Model4, modList2, mxFitFunctionMultigroup(modNames2))
container2 <- mxOption(container2, 'Checkpoint Units', 'minutes')
container2 <- mxOption(container2, 'Checkpoint Count', 1)
containerRun2 <- mxRun(container2, intervals=FALSE, checkpoint=TRUE) 
smr2 <- summary(containerRun2)

save.image(file = "model2.RData")
rm(list = setdiff(ls(), c(ObjectsKeep,"ObjectsKeep", "smr1", "smr2")))

##### Model excluding mt 

Model5 <- mxModel(
      "ModelFive",
      mxMatrix(type = "Full", nrow = 1, ncol = 1, free = TRUE, values = ad2[comb]*totalVar, labels = "vad", name = "Vad", lbound = 1e-10),
#      mxMatrix(type = "Full", nrow = 1, ncol = 1, free = TRUE, values = dd2[comb]*totalVar, labels = "vdd", name = "Vdd", lbound = 1e-10),
      mxMatrix(type = "Full", nrow = 1, ncol = 1, free = TRUE, values = cn2[comb]*totalVar, labels = "vcn", name = "Vcn", lbound = 1e-10),
      mxMatrix(type = "Full", nrow = 1, ncol = 1, free = TRUE, values = ce2[comb]*totalVar, labels = "vce", name = "Vce", lbound = 1e-10),
      #mxMatrix(type = "Full", nrow = 1, ncol = 1, free = TRUE, values = mt2[comb]*totalVar, labels = "vmt", name = "Vmt", lbound = 1e-10),
      mxMatrix(type = "Full", nrow = 1, ncol = 1, free = TRUE, values = am2[comb]*totalVar, labels = "vam", name = "Vam", lbound = 1e-10),
      mxMatrix(type = "Full", nrow = 1, ncol = 1, free = TRUE, values = ee2[comb]*totalVar, labels = "ver", name = "Ver", lbound = 1e-10)
)

ll3 <- list()
for (i in 1:numfam){
      ll3[[i]] <- dat[i,]
}

modList3 <- list()
modNames3 <- paste0("fam", 1:numfam)

for(afam3 in 1:numfam){
      ytemp3 <- paste('S', rownames (Addmat))
      fsize3 <- nrow(Addmat)
      modList3[[afam3]] <- mxModel(name=modNames3[afam3],
                                   mxMatrix("Iden", nrow=fsize3, ncol=fsize3, name="I"), 
                                   mxMatrix("Unit", nrow=fsize3, ncol=fsize3, name='U'),
                                   mxMatrix("Symm", nrow=fsize3, ncol=fsize3, values=Addmat, name="A"), 
#                                   mxMatrix("Symm", nrow=fsize3, ncol=fsize3, values=Dmgmat, name="D"), 
                                   mxMatrix("Symm", nrow=fsize3, ncol=fsize3, values=Nucmat, name="Cn"), 
                                   mxMatrix("Symm", nrow=fsize3, ncol=fsize3, values=Extmat, name="Ce"), 
                                   mxMatrix("Symm", nrow=fsize3, ncol=fsize3, values=Amimat, name="Am"), 
                                   #mxMatrix("Symm", nrow=fsize3, ncol=fsize3, values=Mtdmat, name="Mt"),
                                   mxData(observed = matrix(ll3[[afam3]], nrow=1, dimnames=list(NULL, ytemp3)), type="raw", sort=FALSE),
                                   mxMatrix('Full', nrow=1, ncol=fsize3, name='M', free=TRUE, labels='meanLI',
                                            dimnames=list(NULL, ytemp3)),
                                   mxAlgebra ((A %x% ModelFive.Vad) 
#                                              + (D %x% ModelFive.Vdd) 
                                              + (Cn %x% ModelFive.Vcn) 
                                              + (U %x% ModelFive.Vce) 
                                              #+ (Mt %x% ModelFive.Vmt) 
                                              + (Am %x% ModelFive.Vam) 
                                              + (I %x% ModelFive.Ver), 
                                              name="V", dimnames=list(ytemp3, ytemp3)),
                                   mxExpectationNormal(covariance='V', means='M'), 
                                   mxFitFunctionML()
      )
}
container3 <- mxModel('Model6', Model5, modList3, mxFitFunctionMultigroup(modNames3))
container3 <- mxOption(container3, 'Checkpoint Units', 'minutes')
container3 <- mxOption(container3, 'Checkpoint Count', 1)
containerRun3 <- mxRun(container3, intervals=FALSE, checkpoint=TRUE) 
smr3 <- summary(containerRun3)

save.image(file = "model3.RData")
rm(list = setdiff(ls(), c(ObjectsKeep,"ObjectsKeep", "smr1", "smr2", "smr3")))


##### Model excluding am 

Model7 <- mxModel(
      "ModelSeven",
      mxMatrix(type = "Full", nrow = 1, ncol = 1, free = TRUE, values = ad2[comb]*totalVar, labels = "vad", name = "Vad", lbound = 1e-10),
#      mxMatrix(type = "Full", nrow = 1, ncol = 1, free = TRUE, values = dd2[comb]*totalVar, labels = "vdd", name = "Vdd", lbound = 1e-10),
      mxMatrix(type = "Full", nrow = 1, ncol = 1, free = TRUE, values = cn2[comb]*totalVar, labels = "vcn", name = "Vcn", lbound = 1e-10),
      mxMatrix(type = "Full", nrow = 1, ncol = 1, free = TRUE, values = ce2[comb]*totalVar, labels = "vce", name = "Vce", lbound = 1e-10),
      mxMatrix(type = "Full", nrow = 1, ncol = 1, free = TRUE, values = mt2[comb]*totalVar, labels = "vmt", name = "Vmt", lbound = 1e-10),
      #mxMatrix(type = "Full", nrow = 1, ncol = 1, free = TRUE, values = am2[comb]*totalVar, labels = "vam", name = "Vam", lbound = 1e-10),
      mxMatrix(type = "Full", nrow = 1, ncol = 1, free = TRUE, values = ee2[comb]*totalVar, labels = "ver", name = "Ver", lbound = 1e-10)
)

ll4 <- list()
for (i in 1:numfam){
      ll4[[i]] <- dat[i,]
}

modList4 <- list()
modNames4 <- paste0("fam", 1:numfam)

for(afam4 in 1:numfam){
      ytemp4 <- paste('S', rownames (Addmat))
      fsize4 <- nrow(Addmat)
      modList4[[afam4]] <- mxModel(name=modNames4[afam4],
                                   mxMatrix("Iden", nrow=fsize4, ncol=fsize4, name="I"), 
                                   mxMatrix("Unit", nrow=fsize4, ncol=fsize4, name='U'),
                                   mxMatrix("Symm", nrow=fsize4, ncol=fsize4, values=Addmat, name="A"), 
#                                   mxMatrix("Symm", nrow=fsize4, ncol=fsize4, values=Dmgmat, name="D"), 
                                   mxMatrix("Symm", nrow=fsize4, ncol=fsize4, values=Nucmat, name="Cn"), 
                                   mxMatrix("Symm", nrow=fsize4, ncol=fsize4, values=Extmat, name="Ce"), 
                                   #mxMatrix("Symm", nrow=fsize4, ncol=fsize4, values=Amimat, name="Am"), 
                                   mxMatrix("Symm", nrow=fsize4, ncol=fsize4, values=Mtdmat, name="Mt"),
                                   mxData(observed = matrix(ll4[[afam4]], nrow=1, dimnames=list(NULL, ytemp4)), type="raw", sort=FALSE),
                                   mxMatrix('Full', nrow=1, ncol=fsize4, name='M', free=TRUE, labels='meanLI',
                                            dimnames=list(NULL, ytemp4)),
                                   mxAlgebra ((A %x% ModelSeven.Vad) 
#                                              + (D %x% ModelSeven.Vdd) 
                                              + (Cn %x% ModelSeven.Vcn) 
                                              + (U %x% ModelSeven.Vce) 
                                              + (Mt %x% ModelSeven.Vmt) 
                                              #+ (Am %x% ModelSeven.Vam) 
                                              + (I %x% ModelSeven.Ver), 
                                              name="V", dimnames=list(ytemp4, ytemp4)),
                                   mxExpectationNormal(covariance='V', means='M'), 
                                   mxFitFunctionML()
      )
}
container4 <- mxModel('Model8', Model7, modList4, mxFitFunctionMultigroup(modNames4))
container4 <- mxOption(container4, 'Checkpoint Units', 'minutes')
container4 <- mxOption(container4, 'Checkpoint Count', 1)
containerRun4 <- mxRun(container4, intervals=FALSE, checkpoint=TRUE) 
smr4 <- summary(containerRun4)

save.image(file = "model4.RData")
rm(list = setdiff(ls(), c("smr1", "smr2", "smr3", "smr4")))
save.image(file = "modelSmr.RData")
rm(list = ls())

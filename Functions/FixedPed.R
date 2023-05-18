### generate 11 fixed ped and save them into a r data file ot used for simulation

source("~/R-Project/BalancedPed/evenInsert.R")
source("~/R-Project/BalancedPed/famSizeCal.R")
source("~/R-Project/BalancedPed/SimPed.R")

set.seed(7)
ped1 <- SimPed(2,4,.5,2/3)
ped2 <- SimPed(3,4,.5,2/3)
ped3 <- SimPed(4,4,.5,2/3)
ped4 <- SimPed(2,5,.5,2/3)
ped5 <- SimPed(3,5,.5,2/3)
ped6 <- SimPed(3,5,.6,2/3)
ped7 <- SimPed(3,6,.5,2/3)
ped8 <- SimPed(3,8,.5,2/3)
ped9 <- SimPed(4,4,.5,.8)
ped10 <- SimPed(8,4,.5,2/3)
ped11 <- SimPed(6,6,.5,2/3)

save.image("~/R-Project/BalancedPed/FixedPed.RData")

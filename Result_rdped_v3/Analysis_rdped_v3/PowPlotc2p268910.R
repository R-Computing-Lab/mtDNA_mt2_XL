library(ggplot2)

powerCal <- function(lamda,df){
    1- pchisq(qchisq(1-.05, df), df, lamda)
    
}
my_palette <- c( "#E41A1C", "#332288", "#E69F00", "#DDCC77", "#377EB8",  "#4DAF4A", "#117A65", "#56B4E9", "#A6CE39", "#A9A9A9","#88CCEE", "#CC6677",  "#AA4499",   "#999933", "#882255", "#984EA3")

demo_n <- 10000
### create a new enviroment for one condition
env_c2p2 <- new.env()
load("~/R-Project/mtDNA_mt2/Result_rdped_v3/c2p2/modelSmr.Rdata", envir = env_c2p2)



meanDiffLL_mtam_c2p2 <- env_c2p2$smr2$Minus2LogLikelihood - env_c2p2$smr1$Minus2LogLikelihood
PN <- ncol(env_c2p2[["smr1"]][["dataSummary"]][["fam1"]])*length(env_c2p2[["smr1"]][["dataSummary"]])
lamdaUnit_c2p2 <- meanDiffLL_mtam_c2p2/PN
SSize_c2p2 <- 1: demo_n
LamdaVec_c2p2 <- lamdaUnit_c2p2*SSize_c2p2
powVec_c2p2 <- as.numeric(lapply(LamdaVec_c2p2,powerCal, df = 1))

df_c2p2 <- data.frame(Nped = SSize_c2p2, 
                      power = powVec_c2p2, 
                      Combination = rep("power2", demo_n))


### create a new enviroment for one condition
env_c2p6 <- new.env()
load("~/R-Project/mtDNA_mt2/Result_rdped_v3/c2p6/modelSmr.Rdata", envir = env_c2p6)


meanDiffLL_mtam_c2p6 <- env_c2p6$smr2$Minus2LogLikelihood - env_c2p6$smr1$Minus2LogLikelihood
PN <- ncol(env_c2p6[["smr1"]][["dataSummary"]][["fam1"]])*length(env_c2p6[["smr1"]][["dataSummary"]])
lamdaUnit_c2p6 <- meanDiffLL_mtam_c2p6/PN
SSize_c2p6 <- 1: demo_n
LamdaVec_c2p6 <- lamdaUnit_c2p6*SSize_c2p6
powVec_c2p6 <- as.numeric(lapply(LamdaVec_c2p6,powerCal, df = 1))

df_c2p6 <- data.frame(Nped = SSize_c2p6, 
                      power = powVec_c2p6, 
                      Combination = rep("power1", demo_n))


### create a new enviroment for one condition
env_c2p8 <- new.env()
load("~/R-Project/mtDNA_mt2/Result_rdped_v3/c2p8/modelSmr.Rdata", envir = env_c2p8)


meanDiffLL_mtam_c2p8 <- env_c2p8$smr2$Minus2LogLikelihood - env_c2p8$smr1$Minus2LogLikelihood
PN <- ncol(env_c2p8[["smr1"]][["dataSummary"]][["fam1"]])*length(env_c2p8[["smr1"]][["dataSummary"]])
lamdaUnit_c2p8 <- meanDiffLL_mtam_c2p8/PN
SSize_c2p8 <- 1: demo_n
LamdaVec_c2p8 <- lamdaUnit_c2p8*SSize_c2p8
powVec_c2p8 <- as.numeric(lapply(LamdaVec_c2p8,powerCal, df = 1))

df_c2p8 <- data.frame(Nped = SSize_c2p8, 
                      power = powVec_c2p8, 
                      Combination = rep("power3", demo_n))

### create a new enviroment for one condition
env_c2p9 <- new.env()
load("~/R-Project/mtDNA_mt2/Result_rdped_v3/c2p9/modelSmr.Rdata", envir = env_c2p9)


meanDiffLL_mtam_c2p9 <- env_c2p9$smr2$Minus2LogLikelihood - env_c2p9$smr1$Minus2LogLikelihood
PN <- ncol(env_c2p9[["smr1"]][["dataSummary"]][["fam1"]])*length(env_c2p9[["smr1"]][["dataSummary"]])
lamdaUnit_c2p9 <- meanDiffLL_mtam_c2p9/PN
SSize_c2p9 <- 1: demo_n
LamdaVec_c2p9 <- lamdaUnit_c2p9*SSize_c2p9
powVec_c2p9 <- as.numeric(lapply(LamdaVec_c2p9,powerCal, df = 1))

df_c2p9 <- data.frame(Nped = SSize_c2p9, 
                      power = powVec_c2p9, 
                      Combination = rep("power4", demo_n))

### create a new enviroment for one condition
env_c2p10 <- new.env()
load("~/R-Project/mtDNA_mt2/Result_rdped_v3/c2p10/modelSmr.Rdata", envir = env_c2p10)


meanDiffLL_mtam_c2p10 <- env_c2p10$smr2$Minus2LogLikelihood - env_c2p10$smr1$Minus2LogLikelihood
PN <- ncol(env_c2p10[["smr1"]][["dataSummary"]][["fam1"]])*length(env_c2p10[["smr1"]][["dataSummary"]])
lamdaUnit_c2p10 <- meanDiffLL_mtam_c2p10/PN
SSize_c2p10 <- 1: demo_n
LamdaVec_c2p10 <- lamdaUnit_c2p10*SSize_c2p10
powVec_c2p10 <- as.numeric(lapply(LamdaVec_c2p10,powerCal, df = 1))

df_c2p10 <- data.frame(Nped = SSize_c2p10, 
                       power = powVec_c2p10, 
                       Combination = rep("power5", demo_n))

### create a data frame for graphs
df_c2 <- rbind(df_c2p2, df_c2p6, df_c2p8,df_c2p9,df_c2p10)
df_c2$Combination <- as.factor(df_c2$Combination)


g1 <-ggplot(data = df_c2)+ geom_line(mapping = aes(x = Nped, y = power, color= Combination), size = 1.5) +
    scale_color_manual(values=my_palette[1:5],
                       name="Pedigree Structures",
                       breaks=c("power1", "power2", "power3", "power4", "power5"),
                       labels=c("k = 3, G = 3, m = 15",
                                "k = 3, G = 4, m = 31",
                                "k = 3, G = 5, m = 67",
                                "k = 3, G = 6, m = 138",
                                "k = 3, G = 8, m = 595")
    )+
    theme(panel.background = element_rect(fill = "transparent"),
          panel.grid = element_line(color = "transparent"),
          axis.line = element_line(size = 1, colour = "black"),
          #axis.line.y = element_blank(),
          axis.text = element_text( color = "black"),
          #axis.text.y = element_blank(),
          #axis.ticks.y = element_blank(),
          text=element_text( family="Calibri",  size = 12),
          legend.spacing = unit(-17,'pt'),
          legend.margin = margin(t=0,b=0,unit='pt'),
          legend.background = element_rect(),
          legend.position=c(.8,.2))+
    xlab("N of Individuals")+
    scale_y_continuous(n.breaks = 6)+
    ylab("Power (mt\u00B2)")+
    geom_hline(yintercept = .8, linetype = 5, size = .8, color = "grey")
# +
# annotate(geom = "text",x = 0.62, y =.92, label = "a\u00B2 = .6", family="Calibri", color = "gray40",size = 3)+
# annotate(geom = "text",x = 0.73, y =.6, label = "a\u00B2 = .4", family="Calibri", color = "gray40",size = 3)+
# annotate(geom = "text",x = 0.9, y =.6, label = "a\u00B2 = .2", family="Calibri", color = "gray40",size = 3)+ 
# annotate(geom = "text",x = 0.9, y =.6, label = "a\u00B2 = .4, d\u00B2 = .1", family="Calibri", color = "gray40",size = 3)
g1

ggsave( "~/R-Project/mtDNA_mt2/Result_rdped_v3/Analysis_rdped_v3/graph_c2p268910.png",g1, width = 6, height = 4.5,  type = "cairo-png", dpi = 900)

#rm(list = ls())

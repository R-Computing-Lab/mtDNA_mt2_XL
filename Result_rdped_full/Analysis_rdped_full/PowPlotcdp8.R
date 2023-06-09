library(ggplot2)

powerCal <- function(lamda,df){
    1- pchisq(qchisq(1-.05, df), df, lamda)
    
}
my_palette <- c( "#E41A1C", "#332288", "#E69F00", "#DDCC77", "#377EB8",  "#4DAF4A", "#117A65", "#56B4E9", "#A6CE39", "#A9A9A9","#88CCEE", "#CC6677",  "#AA4499",   "#999933", "#882255", "#984EA3")

demo_n <- 10000
### create a new enviroment for one condition
env_c2p8 <- new.env()
load("~/R-Project/mtDNA_mt2/Result_rdped_full/c2p8/modelSmr.Rdata", envir = env_c2p8)

meanDiffLL_mtam_c2p8 <- env_c2p8$smr2$Minus2LogLikelihood - env_c2p8$smr1$Minus2LogLikelihood
PN <- ncol(env_c2p8[["smr1"]][["dataSummary"]][["fam1"]])*length(env_c2p8[["smr1"]][["dataSummary"]])
lamdaUnit_c2p8 <- meanDiffLL_mtam_c2p8/PN
SSize_c2p8 <- 1: demo_n
LamdaVec_c2p8 <- lamdaUnit_c2p8*SSize_c2p8
powVec_c2p8 <- as.numeric(lapply(LamdaVec_c2p8,powerCal, df = 1))

df_c2p8 <- data.frame(Nped = SSize_c2p8, 
                      power = powVec_c2p8, 
                      Combination = rep("power1", demo_n))


### create a new enviroment for one condition
env_c13p8 <- new.env()
load("~/R-Project/mtDNA_mt2/Result_rdped_full/c13p8/modelSmr.Rdata", envir = env_c13p8)


meanDiffLL_mtam_c13p8 <- env_c13p8$smr2$Minus2LogLikelihood - env_c13p8$smr1$Minus2LogLikelihood
PN <- ncol(env_c13p8[["smr1"]][["dataSummary"]][["fam1"]])*length(env_c13p8[["smr1"]][["dataSummary"]])
lamdaUnit_c13p8 <- meanDiffLL_mtam_c13p8/PN
SSize_c13p8 <- 1: demo_n
LamdaVec_c13p8 <- lamdaUnit_c13p8*SSize_c13p8
powVec_c13p8 <- as.numeric(lapply(LamdaVec_c13p8,powerCal, df = 1))

df_c13p8 <- data.frame(Nped = SSize_c13p8, 
                       power = powVec_c13p8, 
                       Combination = rep("power2", demo_n))



### create a data frame for graphs
df_c2 <- rbind(df_c2p8, df_c13p8)
df_c2$Combination <- as.factor(df_c2$Combination)


g1 <-ggplot(data = df_c2)+ geom_line(mapping = aes(x = Nped, y = power, color= Combination), size = 1.5) +
    scale_color_manual(values=my_palette[1:2],
                       name="Pedigree Structures",
                       breaks=c("power1", "power2"),
                       labels=c("d\u00B2 = 0",
                                "d\u00B2 = .1")
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
    xlab("n of Individuals")+
    scale_y_continuous(n.breaks = 6)+
    ylab("Power (mt\u00B2)")+
    geom_hline(yintercept = .8, linetype = 5, size = .8, color = "grey")
# +
# annotate(geom = "text",x = 0.62, y =.92, label = "a\u00B2 = .6", family="Calibri", color = "gray40",size = 3)+
# annotate(geom = "text",x = 0.73, y =.6, label = "a\u00B2 = .4", family="Calibri", color = "gray40",size = 3)+
# annotate(geom = "text",x = 0.9, y =.6, label = "a\u00B2 = .2", family="Calibri", color = "gray40",size = 3)+ 
# annotate(geom = "text",x = 0.9, y =.6, label = "a\u00B2 = .4, d\u00B2 = .1", family="Calibri", color = "gray40",size = 3)
g1

ggsave( "~/R-Project/mtDNA_mt2/Result_rdped_full/Analysis_rdped_full/graph_cdp8.png",g1, width = 6, height = 4.5,  type = "cairo-png", dpi = 900)

#rm(list = ls())

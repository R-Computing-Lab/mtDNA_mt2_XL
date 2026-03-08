library(ggplot2)

powerCal <- function(lamda,df){
    1- pchisq(qchisq(1-.05, df), df, lamda)
    
}

# Define warm colors
warm_colors <-colorRampPalette(c("darkred", "lightsalmon", "#A6CE39"))(5)

#my_palette <- c( "#E41A1C", "#332288", "#E69F00", "#DDCC77", "#377EB8",  "#4DAF4A", "#117A65", "#56B4E9", "#A6CE39", "#A9A9A9","#88CCEE", "#CC6677",  "#AA4499",   "#999933", "#882255", "#984EA3")

my_palette <- warm_colors

demo_n <- 10000

### create a new enviroment for one condition
env_c1p10 <- new.env()
load("~/R-Project/mtDNA_mt2/Result_rdped20k/c1p10/modelSmr.Rdata", envir = env_c1p10)



meanDiffLL_mtam_c1p10 <- env_c1p10$smr2$Minus2LogLikelihood - env_c1p10$smr1$Minus2LogLikelihood
PN <- ncol(env_c1p10[["smr1"]][["dataSummary"]][["fam1"]])*length(env_c1p10[["smr1"]][["dataSummary"]])
lamdaUnit_c1p10 <- meanDiffLL_mtam_c1p10/PN
SSize_c1p10 <- 1: demo_n
LamdaVec_c1p10 <- lamdaUnit_c1p10*SSize_c1p10
powVec_c1p10 <- as.numeric(lapply(LamdaVec_c1p10,powerCal, df = 1))

df_c1p10 <- data.frame(Nped = SSize_c1p10, 
                      power = powVec_c1p10, 
                      Combination = rep("power1", demo_n))


### create a new enviroment for one condition
env_c2p10 <- new.env()
load("~/R-Project/mtDNA_mt2/Result_rdped20k/c2p10/modelSmr.Rdata", envir = env_c2p10)



meanDiffLL_mtam_c2p10 <- env_c2p10$smr2$Minus2LogLikelihood - env_c2p10$smr1$Minus2LogLikelihood
PN <- ncol(env_c2p10[["smr1"]][["dataSummary"]][["fam1"]])*length(env_c2p10[["smr1"]][["dataSummary"]])
lamdaUnit_c2p10 <- meanDiffLL_mtam_c2p10/PN
SSize_c2p10 <- 1: demo_n
LamdaVec_c2p10 <- lamdaUnit_c2p10*SSize_c2p10
powVec_c2p10 <- as.numeric(lapply(LamdaVec_c2p10,powerCal, df = 1))

df_c2p10 <- data.frame(Nped = SSize_c2p10, 
                      power = powVec_c2p10, 
                      Combination = rep("power2", demo_n))

### create a new enviroment for one condition
env_c3p10 <- new.env()
load("~/R-Project/mtDNA_mt2/Result_rdped20k/c3p10/modelSmr.Rdata", envir = env_c3p10)



meanDiffLL_mtam_c3p10 <- env_c3p10$smr2$Minus2LogLikelihood - env_c3p10$smr1$Minus2LogLikelihood
PN <- ncol(env_c3p10[["smr1"]][["dataSummary"]][["fam1"]])*length(env_c3p10[["smr1"]][["dataSummary"]])
lamdaUnit_c3p10 <- meanDiffLL_mtam_c3p10/PN
SSize_c3p10 <- 1: demo_n
LamdaVec_c3p10 <- lamdaUnit_c3p10*SSize_c3p10
powVec_c3p10 <- as.numeric(lapply(LamdaVec_c3p10,powerCal, df = 1))

df_c3p10 <- data.frame(Nped = SSize_c3p10, 
                      power = powVec_c3p10, 
                      Combination = rep("power3", demo_n))

### create a new enviroment for one condition
env_c4p10 <- new.env()
load("~/R-Project/mtDNA_mt2/Result_rdped20k/c4p10/modelSmr.Rdata", envir = env_c4p10)



meanDiffLL_mtam_c4p10 <- env_c4p10$smr2$Minus2LogLikelihood - env_c4p10$smr1$Minus2LogLikelihood
PN <- ncol(env_c4p10[["smr1"]][["dataSummary"]][["fam1"]])*length(env_c4p10[["smr1"]][["dataSummary"]])
lamdaUnit_c4p10 <- meanDiffLL_mtam_c4p10/PN
SSize_c4p10 <- 1: demo_n
LamdaVec_c4p10 <- lamdaUnit_c4p10*SSize_c4p10
powVec_c4p10 <- as.numeric(lapply(LamdaVec_c4p10,powerCal, df = 1))

df_c4p10 <- data.frame(Nped = SSize_c4p10, 
                      power = powVec_c4p10, 
                      Combination = rep("power4", demo_n))

### create a new enviroment for one condition
env_c5p10 <- new.env()
load("~/R-Project/mtDNA_mt2/Result_rdped20k/c5p10/modelSmr.Rdata", envir = env_c5p10)



meanDiffLL_mtam_c5p10 <- env_c5p10$smr2$Minus2LogLikelihood - env_c5p10$smr1$Minus2LogLikelihood
PN <- ncol(env_c5p10[["smr1"]][["dataSummary"]][["fam1"]])*length(env_c5p10[["smr1"]][["dataSummary"]])
lamdaUnit_c5p10 <- meanDiffLL_mtam_c5p10/PN
SSize_c5p10 <- 1: demo_n
LamdaVec_c5p10 <- lamdaUnit_c5p10*SSize_c5p10
powVec_c5p10 <- as.numeric(lapply(LamdaVec_c5p10,powerCal, df = 1))

df_c5p10 <- data.frame(Nped = SSize_c5p10, 
                      power = powVec_c5p10, 
                      Combination = rep("power5", demo_n))


# Overall
df_p10 <- rbind(df_c1p10, df_c2p10, df_c3p10, df_c4p10, df_c5p10)
df_p10$Combination <- as.factor(df_p10$Combination)

#graph
g1 <-ggplot(data = df_p10)+ geom_line(mapping = aes(x = Nped, y = power, color= Combination), size = 1.5) +
    scale_color_manual(values=my_palette[1:5],
                       name="Variance Combinations",
                       breaks=c("power1", "power2", "power3", "power4", "power5"),
                       labels=c("mt\u00B2 = .100",
                                "mt\u00B2 = .050", 
                                "mt\u00B2 = .010",
                                "mt\u00B2 = .005",
                                "mt\u00B2 = .001")
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
    xlab("N of Pedigrees")+
    scale_y_continuous(n.breaks = 6)+
    ylab("Power (mt\u00B2+j\u00B2)")+
    geom_hline(yintercept = .8, linetype = 5, size = .8, color = "grey")
# +
# annotate(geom = "text",x = 0.62, y =.92, label = "a\u00B2 = .6", family="Calibri", color = "gray40",size = 3)+
# annotate(geom = "text",x = 0.73, y =.6, label = "a\u00B2 = .4", family="Calibri", color = "gray40",size = 3)+
# annotate(geom = "text",x = 0.9, y =.6, label = "a\u00B2 = .2", family="Calibri", color = "gray40",size = 3)+ 
# annotate(geom = "text",x = 0.9, y =.6, label = "a\u00B2 = .4, d\u00B2 = .1", family="Calibri", color = "gray40",size = 3)
g1


ggsave( "~/R-Project/mtDNA_mt2/Result_rdped20k/Analysis_rdped20k/graph_p10c12345.png",g1, width = 6, height = 4.5,  type = "cairo-png", dpi = 900)

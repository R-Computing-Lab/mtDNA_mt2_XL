library(ggplot2)

powerCal <- function(lamda,df){
    1- pchisq(qchisq(1-.05, df), df, lamda)
    
}

# Define warm colors
warm_colors <-colorRampPalette(c("darkred", "lightsalmon", "#A6CE39"))(5)

#my_palette <- c( "#E41A1C", "#332288", "#E69F00", "#DDCC77", "#377EB8",  "#4DAF4A", "#117A65", "#56B4E9", "#A6CE39", "#A9A9A9","#88CCEE", "#CC6677",  "#AA4499",   "#999933", "#882255", "#984EA3")

my_palette <- warm_colors

demo_n <- 20000

### create a new enviroment for one condition
env_c7p6 <- new.env()
load("~/R-Project/mtDNA_mt2/Result_rdped20k/c7p6/modelSmr.Rdata", envir = env_c7p6)



meanDiffLL_mtam_c7p6 <- env_c7p6$smr2$Minus2LogLikelihood - env_c7p6$smr1$Minus2LogLikelihood
PN <- ncol(env_c7p6[["smr1"]][["dataSummary"]][["fam1"]])*length(env_c7p6[["smr1"]][["dataSummary"]])
lamdaUnit_c7p6 <- meanDiffLL_mtam_c7p6/PN
SSize_c7p6 <- 1: demo_n
LamdaVec_c7p6 <- lamdaUnit_c7p6*SSize_c7p6
powVec_c7p6 <- as.numeric(lapply(LamdaVec_c7p6,powerCal, df = 1))

df_c7p6 <- data.frame(Nped = SSize_c7p6, 
                      power = powVec_c7p6, 
                      Combination = rep("power1", demo_n))


### create a new enviroment for one condition
env_c3p6 <- new.env()
load("~/R-Project/mtDNA_mt2/Result_rdped20k/c3p6/modelSmr.Rdata", envir = env_c3p6)



meanDiffLL_mtam_c3p6 <- env_c3p6$smr2$Minus2LogLikelihood - env_c3p6$smr1$Minus2LogLikelihood
PN <- ncol(env_c3p6[["smr1"]][["dataSummary"]][["fam1"]])*length(env_c3p6[["smr1"]][["dataSummary"]])
lamdaUnit_c3p6 <- meanDiffLL_mtam_c3p6/PN
SSize_c3p6 <- 1: demo_n
LamdaVec_c3p6 <- lamdaUnit_c3p6*SSize_c3p6
powVec_c3p6 <- as.numeric(lapply(LamdaVec_c3p6,powerCal, df = 1))

df_c3p6 <- data.frame(Nped = SSize_c3p6, 
                      power = powVec_c3p6, 
                      Combination = rep("power2", demo_n))

### create a new enviroment for one condition
env_c10p6 <- new.env()
load("~/R-Project/mtDNA_mt2/Result_rdped20k/c10p6/modelSmr.Rdata", envir = env_c10p6)



meanDiffLL_mtam_c10p6 <- env_c10p6$smr2$Minus2LogLikelihood - env_c10p6$smr1$Minus2LogLikelihood
PN <- ncol(env_c10p6[["smr1"]][["dataSummary"]][["fam1"]])*length(env_c10p6[["smr1"]][["dataSummary"]])
lamdaUnit_c10p6 <- meanDiffLL_mtam_c10p6/PN
SSize_c10p6 <- 1: demo_n
LamdaVec_c10p6 <- lamdaUnit_c10p6*SSize_c10p6
powVec_c10p6 <- as.numeric(lapply(LamdaVec_c10p6,powerCal, df = 1))

df_c10p6 <- data.frame(Nped = SSize_c10p6, 
                       power = powVec_c10p6, 
                       Combination = rep("power3", demo_n))

### create a new enviroment for one condition
env_c4p6 <- new.env()
load("~/R-Project/mtDNA_mt2/Result_rdped20k/c4p6/modelSmr.Rdata", envir = env_c4p6)


# Overall
df_p6 <- rbind(df_c7p6, df_c3p6, df_c10p6)
df_p6$Combination <- as.factor(df_p6$Combination)

#graph
g1 <-ggplot(data = df_p6)+ geom_line(mapping = aes(x = Nped, y = power, color= Combination), size = 1.5) +
    scale_color_manual(values=my_palette[1:3],
                       name="Variance Combinations",
                       breaks=c("power3", "power2", "power1"),
                       labels=c("a\u00B2 = .6",
                                "a\u00B2 = .4", 
                                "a\u00B2 = .2")
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
    ylab("Power (mt\u00B2)")+
    geom_hline(yintercept = .8, linetype = 5, size = .8, color = "grey")
# +
# annotate(geom = "text",x = 0.62, y =.92, label = "a\u00B2 = .6", family="Calibri", color = "gray40",size = 3)+
# annotate(geom = "text",x = 0.73, y =.6, label = "a\u00B2 = .4", family="Calibri", color = "gray40",size = 3)+
# annotate(geom = "text",x = 0.9, y =.6, label = "a\u00B2 = .2", family="Calibri", color = "gray40",size = 3)+ 
# annotate(geom = "text",x = 0.9, y =.6, label = "a\u00B2 = .4, d\u00B2 = .1", family="Calibri", color = "gray40",size = 3)
g1


ggsave( "~/R-Project/mtDNA_mt2/Result_rdped20k/Analysis_rdped20k/graph_p6c7310.png",g1, width = 6, height = 4.5,  type = "cairo-png", dpi = 900)

## Run the simulations

source("~/R-Project/mtDNA_mt2/InitialData/InitialDataPrep_R2.R")
source("~/R-Project/mtDNA_mt2/Functions/RunSim_rdped.R")

for(i in 6:11){
    
    for(j in 1: 10){
        print(Sys.time())
        cat(paste("start c",i,"p",j, "\n"))
        target_folder <- paste0("~/R-Project/mtDNA_mt2/Result_rdped","/","c",i,"p",j) 
        if (!dir.exists(target_folder)){
            dir.create(target_folder)
        }
        # do.call(RunSim, as.list(df_var[i,],
        #                         paste0("ped",j),
        #                         target_folder))
        RunSim_rd(Var = df_var[i,],
                  kpc = df_ped$k[j],
                  Ngen = df_ped$G[j],
                  sexR = df_ped$p[j],
                  marR = df_ped$r[j],
                  n = 10000,
                  path_results = target_folder)
        print(Sys.time())
        cat(paste("end c",i,"p",j, "\n"))
    }
}

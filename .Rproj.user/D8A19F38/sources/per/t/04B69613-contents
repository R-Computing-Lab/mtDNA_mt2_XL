## Run the simulations

source("~/R-Project/mtDNA_mt2/InitialData/InitialDataPrep.R")
source("~/R-Project/mtDNA_mt2/Functions/RunSim.R")

for(i in 1: nrow(df_var)){
    
    for(j in 1: nrow(df_ped)){
        
        target_folder <- paste0("~/R-Project/mtDNA_mt2/Results","/","c",i,"p",j)
        if (!dir.exists(target_folder)){
            dir.create(target_folder)
        } else {
            
        }
        RunSim(Var = df_var[i,],
               Ped = paste0("ped",j),
               path_results = target_folder )
    }
}

# dir.create(paste0("~/R-Project/mtDNA_mt2/Results","/","c1p1"))
# path_save <- ""
# RunSim(Var = df_var[1,],
#        Ped = ped1,
#        path_results = )

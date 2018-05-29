#' Code to build emulator based on Atlantis runs for
#' GOM sensitivity analysis
#' @author Hem Nalini Morzaria Luna
#' @date June 2017

#library("devtools")
#devtools::install_github('topepo/caret/pkg/caret')

# List of packages for session
.packages = c(
  "readr",
  "stringi",
  "data.table",
  "plyr",
  "tidyverse",
  "magrittr",
  "R.utils",
  "future",
  "doSNOW",
  "dismo",
  "gbm",
  "ggthemes",
  "caret",
  "stringr"
)

# Install CRAN packages (if not already installed)
.inst <- .packages %in% installed.packages()
if (length(.packages[!.inst]) > 0)
  install.packages(.packages[!.inst], dependencies = TRUE)

# Load packages into session
lapply(.packages, require, character.only = TRUE)

# clean up the space
rm(list = ls())
#specify directories
work.folder <- "~/GOMcode"
data.folder <- "~/fullemulator"
model.folder <- "~/modelres"
save.folder <- "~/modelpreds"
data.folder.single <- "~/fullemulatorsingle"

dir.create(save.folder)
dir.create(data.folder.single)

setwd(work.folder)
target.time <- 365
#possible timesteps, for each year. Original output is 30 day intervals
#c(365,730,1095,1460,1825,2190,2555,2920,3285,3650,4015,4380,4745,5110,5474.5)

# species names and guilds
species.names <-
  fread("Group_Names.csv", header = TRUE) %>% tbl_df %>%
  mutate(`predator name` = tolower(`Long Name`)) %>% 
  filter(isFish==1 | `Group Description`=="Squid" | str_detect(`Group Description`, "Shrimp"), Guild!="Other fish") %>% 
  dplyr::select(`predator name`, Code, Guild) %>%
  setNames(c("group_name", "group", "guild")) 

species.groups <- species.names %>% 
  .$group

get.variables <- function(x){
  setwd(model.folder)
  
  print(x)
  model.fits <- c("glm", "elm", "bcart", "rf", "brt", "bmars")
  
  catch.file <-
    paste("*",
          x,
          model.fits[1],
          "catch",
          "model.rda",
          sep = "_")
  #test for catch models for first model
  
  catch.test <- list.files(model.folder, pattern = catch.file, full.names = TRUE)
  
  if (length(catch.test) == 0) {
    variables <- c(paste("/home/atlantis/modelres/",paste(target.time,x,model.fits,paste("biomass",sep="_"),"model.rda",sep="_"),sep=""))
  } else{
    variables <- c((paste("/home/atlantis/modelres/",paste(target.time,x,model.fits,"catch","model.rda",sep="_"),sep="")),paste("/home/atlantis/modelres/",paste(target.time,x,model.fits,"biomass","model.rda",sep="_"),sep=""))
  }
  
  variables
  
}

var.list <- lapply(species.groups,get.variables)

setwd(save.folder)
save(var.list, file=paste(target.time,"varlist.RData",sep="_"))

load.models <- function(thislist, thisvar){
  
  setwd(model.folder)
  
  thesefiles <- thislist %>% unlist %>% 
    grep(thisvar,.,value=TRUE)
  
  if(length(thesefiles)!=0){
    print(thesefiles)
    this.species  <- thesefiles[1] %>% str_split(., "_") %>% unlist %>% .[2]
    
    
    base::load(thesefiles[1])
    base::load(thesefiles[2])
    base::load(thesefiles[3])
    base::load(thesefiles[4])
    base::load(thesefiles[5])
    base::load(thesefiles[6])
    
    model.ensemble <-  list(
      glm = fit.glm.sel,
      elm = fit.elm.sel,
      bcart = fit.bcart.sel,
      rf = fit.rf.sel,
      brt = fit.brt.sel,
      bmars = fit.bmars.sel
    )
    
    setwd(save.folder)
    save(
      model.ensemble,
      file = paste(
        thisvar,
        target.time,
        this.species,
        "model_ensemble.rda",
        sep = "_"
      ))
    print(this.species)
    rm(model.ensemble)
    this.species
    
  } else {
    
    a <- NA
    a
  }
  
}

catch.models <- mclapply(var.list,load.models,thisvar="catch")
biomass.models <- mclapply(var.list,load.models,thisvar="biomass")

data.set <- fread(paste(work.folder,"/emulator_set.csv",sep=""), verbose=TRUE) %>% 
  as.data.frame() %>% 
  as.list

model.files <- list.files(save.folder, pattern = "*.model_ensemble.rda", full.names = TRUE)

predict.models <- function(this.file){
  
  # print(head(data.set))
  
  this.species  <- this.file %>% str_split(., "_") %>% unlist %>% .[3]
  this.var  <- this.file %>% str_split(., "_") %>% unlist %>% .[1] %>% 
    str_split(., "[/]") %>% unlist %>% .[5]
 
  setwd(save.folder)
  load(this.file)
  
  preds.log <-
    predict(model.ensemble, data.set) %>% 
    bind_rows() %>% 
    tbl_df %>% 
    mutate(species=this.species,time=target.time, var=this.var)
  
  rm(model.ensemble)
  preds.log			
}

model.preds <- mclapply(model.files, predict.models)

save(model.preds, file=paste(target.time,"catchpreds.RData",sep="_"))

model.list <- list.files(pattern= "catchpreds.RData")

pred.fun <- function(thismodel){
  
  load(thismodel)
  
  target.time <- thismodel %>% str_split(.,"_") %>% 
    unlist %>% .[1]
  
  model.vc <- model.preds %>% 
    rbindlist() %>%
    tbl_df %>%
    gather(model, pred, -species, -time, -var)
  
  setwd(save.folder)
  fwrite(model.vc,paste(target.time,"model_preds.csv",sep="_"))
}

mclapply(model.list,pred.fun)




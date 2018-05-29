#' Code to build emulator based on Atlantis runs for
#' GOM sensitivity analysis
#' @author Hem Nalini Morzaria Luna
#' @date June 2017

# List of packages for session
.packages = c("readr","stringi","data.table",
              "tidyverse","stringr","R.utils","magrittr","future",
              "parallel","doSNOW","dismo","gbm","plyr","dplyr","sparklyr",
              "ggthemes")

# Install CRAN packages (if not already installed)
.inst <- .packages %in% installed.packages()
if(length(.packages[!.inst]) > 0) install.packages(.packages[!.inst], dependencies = TRUE)

# Load packages into session 
lapply(.packages, require, character.only=TRUE)


# clean up the space

rm(list=ls())

result.folders <- "~/biomassres"
work.folder <- "~/GOMcode"
prm.folder <- "~/GOM_Diets/prms"
save.folder <- "~/atlantisres"

#result.folders <- "G:/Climate_data_2017/GOMAtlantis_Sensitivity_results/Results"
#work.folder <- "E:/Archivos/1Archivos/Articulos/En preparacion/GOM/Code"
#prm.folder <- "E:/Archivos/1Archivos/Articulos/En preparacion/GOM/Code/prms"

setwd(work.folder)

prey.names <- fread("prey_names.csv", header=TRUE) %>% dplyr::select(prey) %>% .$prey %>% as.character() %>% tolower()

functional.groups <- fread("Group_Names.csv", header=TRUE) %>% tbl_df %>% 
  mutate(`predator name`=tolower(`Long Name`)) %>% 
  dplyr::select(`predator name`, Code) %>% 
  setNames(c("group_name", "group"))

#get names of predators that were randomized and use those as response variables in emulator

predator.names <- fread("predator_names.csv", header= TRUE) %>% 
  dplyr::select(predator) %>% .$predator %>% tolower()

#will use last 5 years of biomass
#target.time <- c(5474.5)
target.time <- c(365,730,1095,1460,1825,2190,2555,2920,3285,3650,4015,4380,4745,5110,5474.5)

train.files <- 1:1000

pathfile <- "emulator_data.csv" # 

for(eachindex in train.files) {

  
  #get biomass file
  this.biom <- paste("GOM_OUTBiomIndx_",eachindex,".txt",sep="")
  
  print(paste(this.biom, eachindex)) 
  
  #read biomass
  setwd(result.folders)
  run.biom <- try(read_delim(this.biom, delim=" ",na=c("NA","N/A",""," "),col_names=TRUE)) %>% 
    as.data.frame() %>% .[,1:92] %>% tbl_df %>% 
    filter(Time %in% target.time) %>% 
    mutate(run=eachindex) %>% 
    dplyr::select(Time, run,everything()) %>% 
    gather(group,biomass,-Time, -run)
  
  #get catch file
  this.catch <- paste("GOM_OUTCatch_",eachindex,".txt",sep="")
  
  print(this.catch) 
  
  #read biomass
  setwd(result.folders)
  run.catch <- fread(this.catch, header=TRUE) %>% 
    dplyr::select(Time:SQU) %>% 
    filter(Time %in% target.time) %>% 
    mutate(run=eachindex) %>% 
    dplyr::select(Time, run,everything()) %>% 
    gather(group,catch, -Time, -run)
  
  all.data <- run.catch %>% 
    full_join(run.biom) %>% 
    mutate(catch=(ifelse(is.na(catch),0,catch))) %>% 
    tbl_df
  
  #read prm for availability matrix
  this.prm <- paste("GOM_PRM_2015_",eachindex,".prm",sep="")
  
  setwd(prm.folder)
  
  #read prm file
  prm.data <-  scan(this.prm, character(0), sep = "\n", skip = 7550, nlines=525)
  
  group.no <- "94"
  
  prey.test <- regexpr(pattern = "pPREY", text = prm.data) %>% 
    gsub("-1",FALSE,.) %>% 
    gsub("1",TRUE,.) %>% 
    as.logical() 
  
  prey.headers <- prm.data[prey.test] %>% 
    as.data.frame() %>% tbl_df %>% 
    na.omit %>% 
    setNames("prey") %>% 
    filter(!grepl("is the availability",prey)) %>% 
    separate(prey,c("name","groups")," ") %>% 
    mutate(name=gsub(",",'',name)) %>% 
    mutate(name=gsub("94",'',name)) %>% 
    dplyr::select(name)
  
  numeric.test <- regexpr(pattern = "pPREY", text = prm.data) %>% 
    gsub("-1",TRUE,.) %>% 
    gsub("1",FALSE,.) %>%
    as.logical() 
  
  all.prey.data <- prm.data[numeric.test] %>% 
    as.data.frame() %>% tbl_df %>% 
    na.omit %>% 
    setNames("prey") %>% 
    separate(prey,prey.names," ") %>% 
    mutate_all(funs(as.numeric)) %>% 
    na.omit %>% 
    bind_cols(prey.headers) %>% 
    dplyr::select(name, everything()) %>% 
    gather(group_name, availability,-name) %>% 
    filter(group_name %in% predator.names) %>% 
    left_join(functional.groups,by="group_name") %>% 
    mutate(prey_group = paste(name,group, sep="_")) %>% 
    dplyr::select(prey_group, availability) %>% 
    filter(availability!=0) %>% 
    spread(prey_group, availability) 
  
  setwd(work.folder)
  
  names.brt.frame <- names(all.prey.data)
  
  if(file.exists("availability_names.csv")==FALSE){
       write.csv(names.brt.frame, "availability_names.csv")
  }
  
  prey.data.brt <- all.prey.data %>% 
    as.matrix
  names(prey.data.brt) <- NULL
  
  system.time(prey.data.mat <- matrix(rep(t(prey.data.brt), 1365) , ncol = ncol(prey.data.brt) , byrow = TRUE) %>% 
    as.data.frame() %>% setNames(names.brt.frame) %>% tbl_df)
  
  #combine biomass and availability matrix
  brt.frame <- all.data %>% 
    bind_cols(prey.data.mat) %>% 
    tbl_df
  
  if(eachindex == 1){
    
    setwd(work.folder)
    write_delim(brt.frame, path=pathfile, delim=",", append = FALSE, col_names = TRUE)
  } else{
    
    setwd(work.folder)
    write_csv(brt.frame, path=pathfile,append = TRUE, col_names = FALSE)
    
  }
  
}


# extract biomass from selected cells, using BoxBiomass

result.folders <- "~/biomasscell/biomasscell"
oil.folder <- "~/GOMAtlantisOil"

# will use polygons most affected by oil, same as Ainsworth et al. 2017
target.cell <- c(20, 21, 23, 31, 38, 39, 43, 52, 53, 54, 56, 63)
target.time <- c(5370.0,5400.0,5430.0,5460.0,5474.5)

setwd(work.folder)

species.names <- fread("Group_Names.csv", header=TRUE) %>% tbl_df %>% 
  mutate(`predator name`=tolower(`Long Name`)) %>% 
  dplyr::select(`predator name`, Code, Guild) %>% 
  setNames(c("group_name", "group", "guild"))


setwd(result.folders) 
test.files <- list.files(pattern="BoxBiomass_")

index.files <- 1:length(test.files)

pathfile <- "cell_biomass.csv"

setwd(nooil.folder)

nooil.cell.biomass <- try(read_delim("GOM_OUTBoxBiomass.txt", delim=" ",na=c("NA","N/A",""," "),col_names=TRUE)) %>% 
  as.data.frame() %>% 
  tbl_df %>% 
  filter(Box %in% target.cell) %>% 
  # filter(Time %in% target.time) %>% 
  dplyr::select(-Box) %>% 
  group_by(Time) %>% 
  summarise_all(funs(sum)) %>% 
  mutate(months = 1:nrow(.)) %>% 
  dplyr::select(months, everything(), -Time) %>% 
  gather(group,nooil_biomass,-months) %>% 
  left_join(species.names,by="group") %>% 
  dplyr::select(months,guild,nooil_biomass) %>% 
  group_by(guild, months) %>% 
  summarise_all(funs(sum))


for(eachindex in index.files){
  
  this.biom <- test.files[eachindex]
  print(paste(this.biom, eachindex, ":",length(test.files)))
        
  run.num <- this.biom %>% 
    str_split(.,"_") %>% unlist() %>% 
    str_split(.,"[.]") %>% unlist() %>% .[3]
    
  #read biomass
  setwd(result.folders)
  run.biom <- try(read_delim(this.biom, delim=" ",na=c("NA","N/A",""," "),col_names=TRUE)) %>% 
    as.data.frame() %>% 
    tbl_df %>% 
    filter(Box %in% target.cell) %>% 
   # filter(Time %in% target.time) %>% 
    dplyr::select(-Box) %>% 
    group_by(Time) %>% 
    summarise_all(funs(sum)) %>% 
    mutate(run=run.num, months = 1:nrow(.)) %>% 
    dplyr::select(run,months, Time, everything()) %>% 
    gather(group,cell_biomass,-run, -months, -Time) %>% 
    rename(time = Time) %>% 
    mutate(time = ifelse(time == 365, 2011, ifelse(
      time == 730, 2012, ifelse(time == 1095, 2013,
                                ifelse(
                                  time == 1460, 2014, ifelse(time == 1825, 2015, ifelse(
                                    time == 2190, 2016, ifelse(time == 2555, 2017,
                                                               ifelse(time == 2920, 2018, ifelse(time == 3285, 2019, ifelse(
                                                                 time == 3650, 2020, ifelse(time == 4015, 2021, ifelse(
                                                                   time == 4380, 2022, ifelse(time == 4745, 2023, ifelse(
                                                                     time == 5110, 2024, ifelse(time == 5474.5, 2025, NA)
                                                                   ))
                                                                 ))
                                                               ))
                                                               ))
                                  ))
                                ))
    ))
    ) %>% 
    left_join(species.names,by="group") %>% 
    dplyr::select(run,months,guild,cell_biomass) %>% 
    group_by(run, guild, months) %>% 
    summarise_all(funs(sum)) %>% 
    left_join(nooil.cell.biomass, by=c("guild", "months")) %>% 
    mutate(rel_biomass = cell_biomass/nooil_biomass)

  if(eachindex == 1){
    
    setwd(work.folder)
    write_delim(run.biom, path=pathfile, delim=",", append = FALSE, col_names = TRUE)
  } else {
    
    setwd(work.folder)
    write_csv(run.biom, path=pathfile, append = TRUE, col_names = FALSE)
    
  }
  
}



setwd(work.folder)

cell.biomass <- fread("cell_biomass.csv", header= TRUE, verbose = TRUE) %>% 
  tbl_df  

max.cell.biomass <- cell.biomass %>% 
  dplyr::group_by(guild, months) %>% 
  dplyr::summarise(max_biomass=max(rel_biomass)) %>% 
  na.omit()
  
min.cell.biomass <- cell.biomass %>% 
  dplyr::group_by(guild, months) %>% 
  dplyr::summarise(min_biomass=min(rel_biomass)) %>% 
  na.omit()

frame.cell.biomass <- cell.biomass %>% tbl_df %>% 
  dplyr::group_by(guild, months) %>% 
  dplyr::summarise(mean_biomass=mean(rel_biomass)) %>% 
  na.omit() %>% 
  left_join(max.cell.biomass, by=c("guild","months")) %>% 
  left_join(min.cell.biomass, by=c("guild","months"))


ribbon.plot <- ggplot(frame.cell.biomass, aes(months))+
  geom_ribbon(aes(ymin = min_biomass, ymax = max_biomass), fill = "darkslateblue") +
  geom_line(aes(y = mean_biomass)) +
  facet_wrap(~ guild) +
  theme_tufte(base_size = 12, base_family='GillSans')+
  ylab("Biomass relative to base model")+
  xlab("Months")
  
setwd(work.folder)
ggsave("relative_biomass.png", ribbon.plot, device = "png", width = 7, height = 5.5, dpi=350)



#' Code to build emulator based on Atlantis runs for
#' GOM sensitivity analysis
#' @author Hem Nalini Morzaria Luna
#' @date June 2017

# List of packages for session
.packages = c(
  "readr",
  "stringi",
  "data.table",
  "tidyverse",
  "R.utils",
  "future",
  "doSNOW",
  "parallel",
  "ggthemes",
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

result.folders <- "~/emulatorpreds"
work.folder <- "~/GOMcode"
data.folder <- "~/emulatorprms"
save.folder <- "~/modelres"
dir.create(save.folder)
dir.create(data.folder)

fit.results <- list.files(save.folder, pattern="*fitresults.csv$")

file.index <- 1:length(fit.results)

nThreads <- detectCores() - 1

target.time <- 365

#c(365,730,1095,1460,1825,2190,2555,2920,3285,3650,4015,4380,4745,5110,5474.5)


cl <- makeSOCKcluster(nThreads)
registerDoSNOW(cl)

rmse.results <- foreach(eachresult=file.index, .verbose=TRUE) %dopar% {
  
  .packages = c(
    "readr",
    "stringi",
    "data.table",
    "tidyverse",
    "R.utils",
    "ggthemes",
    "stringr"
  )
  
  
  # Install CRAN packages (if not already installed)
  .inst <- .packages %in% installed.packages()
  if (length(.packages[!.inst]) > 0)
    install.packages(.packages[!.inst], dependencies = TRUE)
  
  # Load packages into session
  lapply(.packages, require, character.only = TRUE)
  
  this.variable <- fit.results[eachresult] %>% str_split(.,"_") %>% 
    unlist %>% .[3]
  
  this.timestep <- fit.results[eachresult]  %>% str_split(.,"_") %>% 
    unlist %>% .[1]
  
 
  setwd(save.folder)
  this.result <- read_csv(fit.results[eachresult]) %>% tbl_df %>% 
    filter(metric=="RMSE") %>% 
    dplyr::group_by(group,model_fit) %>% 
    dplyr::summarise(RMSE=min(values)) %>% 
    mutate(variable=this.variable, time = this.timestep, index = eachresult)
  
}

rmse.table <- rmse.results %>% rbindlist()

this.file <- paste(work.folder,"/",target.time,"_varfit.csv",sep="")
fwrite(rmse.table,this.file)


#' Summarise predicted vs observed

# clean up the space
rm(list = ls())

result.folders <- "~/emulatorpreds"
work.folder <- "~/GOMcode"
save.folder <- "~/modelres"

fit.results <- list.files(save.folder, pattern="*pred.csv$")

file.index <- 1:length(fit.results)

target.time <- 365

#c(365,730,1095,1460,1825,2190,2555,2920,3285,3650,4015,4380,4745,5110,5474.5)


nThreads <- detectCores() - 1
cl <- makeSOCKcluster(nThreads)
registerDoSNOW(cl)

rmse.results <- foreach(eachresult=file.index, .verbose=TRUE) %dopar% {
  
  .packages = c(
    "readr",
    "stringi",
    "data.table",
    "tidyverse",
    "R.utils",
    "ggthemes",
    "stringr"
  )
  
  
  # Install CRAN packages (if not already installed)
  .inst <- .packages %in% installed.packages()
  if (length(.packages[!.inst]) > 0)
    install.packages(.packages[!.inst], dependencies = TRUE)
  
  # Load packages into session
  lapply(.packages, require, character.only = TRUE)
  
  this.variable <- fit.results[eachresult] %>% str_split(.,"_") %>% 
    unlist %>% .[3]
  
  this.timestep <- fit.results[eachresult]  %>% str_split(.,"_") %>% 
    unlist %>% .[1]
  
  setwd(save.folder)
  
  this.result <- read_csv(fit.results[eachresult]) %>% tbl_df
  
  
}

rmse.table <- rmse.results %>% rbindlist()

this.file <- paste(work.folder,"/",target.time,"_predobs.csv",sep="")
fwrite(rmse.table,this.file)



#' Code to build point figure of RMSE values for
#' GOM sensitivity analysis
#' @author Hem Nalini Morzaria Luna
#' @date June 2017

# List of packages for session
.packages = c(
  "readr",
  "stringi",
  "data.table",
  "tidyverse",
  "R.utils",
  "ggthemes",
  "stringr",
  "ggpubr",
  "plotrix",
  "openair",
  "png",
  "gridExtra")


# Install CRAN packages (if not already installed)
.inst <- .packages %in% installed.packages()
if (length(.packages[!.inst]) > 0)
  install.packages(.packages[!.inst], dependencies = TRUE)

# Load packages into session
lapply(.packages, require, character.only = TRUE)

# clean up the space
rm(list = ls())

result.folders <- "~/emulatorpreds"
work.folder <- "~/GOMcode"
data.folder <- "~/emulatorprms"
save.folder <- "~/modelres"
emulator.folder <- "~/modelpreds"
dir.create(save.folder)
dir.create(data.folder)

species.groups <-
  fread(paste(work.folder, "Group_Names.csv", sep = "/")) %>%
  tbl_df %>%
  rename(group = Code)

#get names of predators that were randomized and use those as response variables in emulator
predator.names <-
  fread(paste(work.folder, "predator_names.csv", sep = "/"), header = TRUE) %>%
  dplyr::select(predator) %>% .$predator %>% tolower()

# concatenate files from all timesteps

fit.results <-
  list.files(work.folder, pattern = "*varfit.csv$", full.names = TRUE)

file.index <- 1:length(fit.results)

#rmse.results <- foreach(eachresult=file.index, .verbose=TRUE) %dopar% {

rmse.results <- function(eachresult) {
  this.timestep <- eachresult  %>% str_split(., "_") %>%
    unlist %>% .[1]
  
  setwd(save.folder)
  
  this.result <- fread(eachresult) %>% tbl_df
  
  
}

rmse.data <- lapply(fit.results, rmse.results) %>%
  rbindlist() %>% tbl_df %>%
  rename(model = model_fit) %>%
  mutate(model = ifelse(
    model == "BaggedCART",
    "BCART",
    ifelse(
      model == "BaggedMultivariateAdaptiveRegressionSplines",
      "BMARS",
      ifelse(
        model == "BoostedGeneralizedAdditiveModel",
        "ELM",
        ifelse(
          model == "BoostedRegressionTrees",
          "BRT",
          ifelse(
            model == "RandomForest",
            "RF",
            ifelse(model == "GeneralizedLinearModel", "GLM", NA)
          )
        )
      )
    )
  )) %>%
  mutate(group = ifelse(group == "Inf", "INF", group)) %>%
  left_join(species.groups, by = "group") %>%   # ifelse(time==3285,8,
  mutate(`description` = tolower(`Group Description`)) %>%
  filter(description %in% predator.names) %>%
  mutate(time = as.factor(time)) %>%
  mutate(cellindex = 1:nrow(.)) %>% 
  rename(guild=Guild)


#get models with best fit and lowest rmse

min.rmse <- rmse.data %>%
  group_by(group, time, variable) %>%
  slice(which.min(RMSE))

# make plot of RMSE values
rmse.plot <- function(thisvariable) {
  rmse.data %>%
    filter(variable == thisvariable) %>%
    filter(guild!="Filter feeders",guild!="Invertebrates",
           guild!="Other fish",guild!="Epibenthos") %>%   #  filter(model != "ELM") %>%
    ggplot(aes(x = model, y = RMSE)) +
    geom_point(size=0.7) +
    #  facet_grid(time ~ ., scales = "free_y", margins = TRUE) +
    labs(x = "Model",
         y = "RMSE")+
    theme(legend.position = "bottom") +
    #  scale_color_ptol("Model") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))+
    theme(axis.text.y = element_text(size = rel(0.85)))
  
  
}

catch.plot <- rmse.plot("catch")
biomass.plot <- rmse.plot("biomass")

rmse.plot <- ggarrange(biomass.plot,catch.plot + rremove("ylab"), 
                       labels = c("Biomass", "Catch"),
                       ncol = 2, nrow = 1)

setwd(work.folder)
ggsave(
  "rmse_plot.png",
  rmse.plot,
  device = "png",
  height = 4,
  width = 6,
  dpi = 150
)
rm(rmse.plot)

# save combined data
this.file <- paste(work.folder,"/","all_pred_obs.csv",sep="")
fwrite(rmse.data,this.file)

#' Taylor diagrams
#' combine observed vs predicted and plot in taylor diagrams to assess model performance

#read files that have predicted and observed data, these are based on test data

fit.results <- list.files(work.folder, pattern="*predobs.csv$", full.names = TRUE)

#these are Taylor diagrams for all models, grouped by functional groups

rmse.results <-function(eachresult){  
  
  this.timestep <- eachresult %>% str_split(.,"/") %>% 
    unlist %>% .[5] %>% 
    str_split(.,"_") %>% 
    unlist %>% .[1]
  
  print(eachresult)
  
  setwd(save.folder)
  
  this.result <- fread(eachresult) %>% tbl_df %>% 
    mutate(model=toupper(model)) %>% 
    mutate(time=as.factor(time)) %>% 
    filter(group_name %in% predator.names) %>% 
    na.omit()
  
  this.min.data <- min.rmse %>% 
    filter(time==this.timestep)
  
  #use set with min rmse value to select only obs and pred values for best model for each species
  min.rmse.data <- this.result %>% 
    left_join(this.min.data, by=c("group","model","variable","time","guild")) %>% 
    na.omit()
  
  #min.rmse.data 
  this.result
  
}

rmse.data <- lapply(fit.results,rmse.results) %>% 
  rbindlist() %>% tbl_df

taylor.plot <- function(thisvariable){
  options(digits=1)
  
  this.data1 <- rmse.data %>% 
    filter(variable==thisvariable) %>% 
    filter(guild!="Filter feeders",guild!="Invertebrates",
           guild!="Other fish",guild!="Epibenthos") %>% 
    mutate(model=as.factor(model)) %>% 
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
    mutate(time=as.factor(time)) %>% 
    as.data.frame()
  
  setwd(work.folder)  
  png(paste(thisvariable,"all_taylordiagram.png",sep="_"),width = 16 * 150, height = 8 * 150, res = 150)
  
  t.dia1 <- TaylorDiagram(this.data1, obs="obs", mod="pred", group="model", 
                          type = "group", normalise = TRUE, cols="Dark2", 
                          rms.col = "darkorange4", key.title = "Model", key.columns = 1,
                          ylab = "Normalized standard deviation",xlab = "Normalized standard deviation",
                          key.position="bottom")
  
  dev.off()
}

taylor.plot("biomass")
taylor.plot("catch")

# save combined data
this.file <- paste(work.folder,"/","mall_pred_obs.csv",sep="")
fwrite(rmse.data,this.file)

#these are Taylor diagrams for all models, grouped by time and guild

rmse.results <-function(eachresult){  
  
  this.timestep <- eachresult %>% str_split(.,"/") %>% 
    unlist %>% .[5] %>% 
    str_split(.,"_") %>% 
    unlist %>% .[1]
  
  setwd(save.folder)
  
  this.result <- fread(eachresult) %>% tbl_df %>% 
    mutate(model=toupper(model)) %>% 
    mutate(time=as.factor(time)) %>% 
    filter(group_name %in% predator.names) %>% 
    na.omit()
  
  this.min.data <- min.rmse %>% 
    filter(time==this.timestep)
  
  #use set with min rmse value to select only obs and pred values for that model for each species
  min.rmse.data <- this.result %>% 
    left_join(this.min.data, by=c("group","model","variable","time","guild")) %>% 
    select(RMSE,everything()) %>% 
    na.omit()
  
  #min.rmse.data 
  this.result
  
}

rmse.data <- lapply(fit.results,rmse.results) %>% 
  rbindlist() %>% tbl_df 

taylor.plot <- function(thisvariable){
  options(digits=1)
  
  this.data1 <- rmse.data %>% 
    filter(variable==thisvariable) %>% 
    filter(guild!="Filter feeders",guild!="Invertebrates",
           guild!="Other fish",guild!="Epibenthos") %>% 
    mutate(model=as.factor(model)) %>% 
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
    mutate(time=as.factor(time)) %>% 
    as.data.frame()
  
  setwd(work.folder)  
  png(paste(thisvariable,"taylordiagram_time.png",sep="_"),width = 16 * 150, height = 8 * 150, res = 150)
  
  t.dia1 <- TaylorDiagram(this.data1, obs="obs", mod="pred", group="time", 
                          type = "guild", normalise = TRUE, cols="Dark2", 
                          rms.col = "darkorange4", key.title = "Year", 
                          key.columns = 1, ylab = "Normalized standard deviation",
                          xlab = "Normalized standard deviation",
                          key.position="bottom")
  
  dev.off()
}

taylor.plot("biomass")
taylor.plot("catch")

#these are Taylor diagrams for all models, grouped by model and guild

taylor.plot <- function(thisvariable){
  options(digits=1)
  
  this.data1 <- rmse.data %>% 
    filter(variable==thisvariable) %>% 
    filter(guild!="Filter feeders",guild!="Invertebrates",
           guild!="Other fish",guild!="Epibenthos") %>% 
    mutate(model=as.factor(model)) %>% 
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
    mutate(time=as.factor(time)) %>% 
    as.data.frame()
  
  var.name <- str_to_title(thisvariable)
  
  setwd(work.folder)  
  png(paste(thisvariable,"taylordiagram_model.png",sep="_"),width = 16 * 150, height = 8 * 150, res = 150)
  
  t.dia1 <- TaylorDiagram(this.data1, obs="obs", mod="pred", group="model", 
                          type = "guild", normalise = TRUE, cols="Dark2", 
                          rms.col = "darkorange4", key.title = "Model", 
                          key.columns = 1, ylab = "Normalized standard deviation",
                          xlab = "Normalized standard deviation",
                          key.pos="right", main=var.name)
  
  dev.off()
}

taylor.plot("biomass")
taylor.plot("catch")

# save combined data
this.file <- paste(work.folder,"/","all_pred_obs.csv",sep="")
fwrite(rmse.data,this.file)

setwd(work.folder)

png("taylordiagram_model.png",width = 8*400, height = 8*400, res = 400)

rl <- lapply(list("biomass_taylordiagram_model.png", "catch_taylordiagram_model.png"), png::readPNG)
gl <- lapply(rl, grid::rasterGrob)
do.call(gridExtra::grid.arrange, gl)
dev.off()

rm(rl, gl)

#' Sensitivity analysis
#' Plot trajectory of biomass, catch and catch per fishery for 1000 Atlantis runs


biomass.catch.res <- "~/biomasscatchRes/biomasscatchOUT"
fishery.catch <- "~/catchfishery"

biomass.res <- list.files(path=biomass.catch.res,pattern="GOM_OUTBiomIndx_*.*$", full.names = TRUE)
catch.res <- list.files(path=biomass.catch.res,pattern="GOM_OUTCatch_*.*$", full.names = TRUE)
fishery.res <- list.files(path=fishery.catch,pattern="GOM_OUTCatchPerFishery_*.*$", full.names = TRUE)

setwd(work.folder)
fleet.names <-  fread("fleetnames.csv", header = TRUE)


target.time <- c(365,730,1095,1460,1825,2190,2555,2920,3285,3650,4015,4380,4745,5110,5474.5)

nooil.biomass <- fread("~/GOMAtlantisNoOil/GOM_OUTBiomIndx.txt", header = TRUE, select=1:91) %>% tbl_df %>% 
  dplyr::filter(Time %in% target.time) %>% 
  dplyr::rename(time=Time) %>% 
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
  gather(group,variable,-time) %>% 
  mutate(var = "biomass")

nooil.catch <- fread("~/GOMAtlantisNoOil/GOM_OUTCatch.txt", header = TRUE, select=1:70) %>% tbl_df %>% 
  dplyr::filter(Time %in% target.time) %>% 
  dplyr::rename(time=Time) %>% 
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
  gather(group,variable,-time) %>% 
  mutate(var="catch")  


base.biomass <- fread("~/GOMAtlantisOil_00/GOM_OUTBiomIndx_00.txt", header = TRUE, select=1:91) %>% tbl_df %>% 
  dplyr::filter(Time %in% target.time) %>% 
  dplyr::rename(time=Time) %>% 
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
  gather(group,variable,-time) %>% 
  mutate(run=0, var = "biomass")


base.catch <- fread("~/GOMAtlantisOil_00/GOM_OUTCatch_00.txt", header = TRUE, select=1:70) %>% tbl_df %>% 
  dplyr::filter(Time %in% target.time) %>% 
  dplyr::rename(time=Time) %>% 
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
  gather(group,variable,-time) %>% 
  mutate(run=0, var="catch")  

base.fishery <- fread("~/GOMAtlantisNoOil/GOM_OUTCatchPerFishery.txt", header = TRUE, select=1:70) %>% tbl_df %>% 
  dplyr::filter(Time %in% target.time) %>% 
  dplyr::rename(time=Time) %>% 
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
  gather(group,variable,-time, -Fishery) %>% 
  mutate(run=0, var="fisherycatch") 

base.data <- bind_rows(base.biomass,base.catch) %>% 
  left_join(nooil.biomass, by=(c("group","time","var"))) %>% 
  rename(variable = variable.x, nooil=variable.y)

get.data <- function(eachfile, endlimit){
  
  run <- eachfile %>% str_split(., pattern = "/") %>%
    unlist() %>% .[6] %>%
    str_split(., pattern = "_") %>%
    unlist() %>% .[3] %>% 
    str_split(., pattern = "[.]") %>% 
    unlist() %>% .[1] %>% as.numeric()
  
  print(run)
  
  this.data <- fread(eachfile, header = TRUE, select=1:endlimit) %>% tbl_df %>% 
    dplyr::filter(Time %in% target.time) %>% 
    dplyr::rename(time=Time) %>% 
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
    gather(group,variable,-time) %>% 
    mutate(run=run)
}

biomass.data <- lapply(biomass.res,get.data, endlimit=91) %>% 
  rbindlist() %>% tbl_df %>% 
  left_join(species.groups, by="group") %>% 
  mutate(group_name=tolower(`Long Name`)) %>% 
  filter(group_name %in% predator.names) %>% 
  rename(guild=Guild)


catch.data <- lapply(catch.res,get.data, endlimit=70) %>% 
  rbindlist() %>% 
  left_join(species.groups, by="group")%>% 
  mutate(group_name=tolower(`Long Name`)) %>% 
  filter(group_name %in% predator.names) %>% 
  rename(guild=Guild)

select.fisheries <- fleet.names %>% 
  distinct(Fishery) %>% .$Fishery

get.data.fish <- function(eachfile){
  
  run <- eachfile %>% str_split(., pattern = "/") %>%
    unlist() %>% .[5] %>%
    str_split(., pattern = "_") %>%
    unlist() %>% .[3] %>% 
    str_split(., pattern = "[.]") %>% 
    unlist() %>% .[1] %>% as.numeric()
  
  print(run)
  
  this.data <- fread(eachfile, header = TRUE) %>% tbl_df %>% 
    dplyr::filter(Time %in% target.time) %>% 
    dplyr::rename(time=Time) %>% 
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
    gather(group,variable,-time, -Fishery) %>% 
    mutate(run=run, var="fisherycatch") %>% 
    filter(Fishery %in% select.fisheries) %>% 
    filter(time==2011)
    
}

options(digits = 4)

fishery.data <- lapply(fishery.res,get.data.fish) %>% 
  rbindlist()
  
sum.fishery.base <- base.fishery %>% 
  filter(time==2011) %>% 
  filter(Fishery %in% select.fisheries) %>% 
  rename(base=variable) %>% 
  select(-run) %>% 
  group_by(time,Fishery) %>% 
  summarise(sum_base = sum(base))
            
sum.fishery <- fishery.data %>% 
  filter(run!=0) %>% 
  group_by(time,Fishery,run) %>% 
  summarise(sum_bio = sum(variable)) %>% 
  left_join(sum.fishery.base, by= c("Fishery","time")) %>% 
   mutate(catch_change = ((sum_bio/sum_base)*100)-100) %>% 
  group_by(time,Fishery) %>% 
  summarise(sum_group = mean(catch_change), min_val = min(catch_change), max_val = max(catch_change)) %>% 
  left_join(fleet.names, by = "Fishery") %>% 
  mutate(ymin_val=min_val - sum_group , ymax_val =sum_group - max_val)
  
sim.hist.plot <- sum.fishery %>%
  ggplot(., aes(`full name`, sum_group))+
  geom_bar(stat="identity", fill = "deepskyblue3")+
  geom_errorbar(aes(ymin=min_val, ymax=max_val), colour="black", width=.5)+
  theme_classic(base_size = 12, base_family='GillSans')+
  ylab(paste("Change in catch (%)"))+
  xlab("Fleet")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

setwd(work.folder)
ggsave("fishery_catch_simulations.png", sim.hist.plot, device = "png", width = 6, height = 4, dpi=300)


plot.trajectory <- function(data, thisvariable, thisbase, thisnooil){
  
  options(digits=4)
  
  thisdata <- data %>% left_join(thisnooil, by=c("group","time")) %>% 
    rename(variable=variable.x, nooil = variable.y) %>% 
    mutate(ratio = variable/nooil) %>% 
    dplyr::select(time,group,ratio,everything())
  
  thisbasedata <- thisbase %>% left_join(thisnooil, by=c("group","time")) %>% 
    rename(variable=variable.x, nooil = variable.y, var = var.x) %>% 
    mutate(ratio = variable/nooil) %>% 
    dplyr::select(time,group,ratio,everything(), -var.y)
  
  
  max.val.guild <- thisdata %>% 
    dplyr::group_by(guild, time) %>% 
    mutate(ratio=ifelse(is.nan(ratio),NA,ratio)) %>%
    na.omit %>% 
    dplyr::summarise(max_val=max(ratio)) %>% 
    na.omit()
  
  min.val.guild <- thisdata %>% 
    dplyr::group_by(guild, time) %>% 
    mutate(ratio=ifelse(is.nan(ratio),NA,ratio)) %>% 
    na.omit %>% 
    dplyr::summarise(min_val=min(ratio)) %>% 
    na.omit()
  
  std.val.guild <- thisdata %>% 
    dplyr::group_by(guild, time) %>% 
    mutate(ratio=ifelse(is.nan(ratio),NA,ratio)) %>% 
    na.omit %>% 
    dplyr::summarise(std_val=std.error(ratio)) %>% 
    na.omit()
  
  bad.guilds <- c("Epibenthos","Filter feeders","Invertebrates","Other fish")
  
  frame.data.guild <- thisdata %>% tbl_df %>% 
    dplyr::group_by(guild, time) %>% 
    filter(!guild %in% bad.guilds) %>% 
    mutate(ratio=ifelse(is.nan(ratio),NA,ratio)) %>% 
    na.omit %>% 
    dplyr::summarise(mean_val=mean(ratio)) %>% 
    na.omit() %>% 
    left_join(max.val.guild, by=c("guild","time")) %>% 
    left_join(min.val.guild, by=c("guild","time")) %>% 
    left_join(std.val.guild, by=c("guild","time"))
  
  # summarize by group
  
  max.val.group <- thisdata %>% 
    mutate(ratio=ifelse(is.nan(ratio),NA,ratio)) %>% 
    na.omit %>% 
    dplyr::group_by(group, time) %>% 
    dplyr::summarise(max_val=max(ratio)) %>% 
    na.omit()
  
  min.val.group <- thisdata %>% 
    mutate(ratio=ifelse(is.nan(ratio),NA,ratio)) %>% 
    na.omit %>% 
    dplyr::group_by(group, time) %>% 
    dplyr::summarise(min_val=min(ratio)) %>% 
    na.omit()
  
  std.val.group <- thisdata %>% 
    mutate(ratio=ifelse(is.nan(ratio),NA,ratio)) %>% 
    na.omit %>% 
    dplyr::group_by(group, time) %>% 
    dplyr::summarise(std_val=std.error(ratio)) %>% 
    na.omit()
  
  frame.data.group <- thisdata %>% tbl_df %>% 
    filter(!guild %in% bad.guilds) %>% 
    mutate(ratio=ifelse(is.nan(ratio),NA,ratio)) %>% 
    na.omit %>% 
    dplyr::group_by(group, time) %>% 
    dplyr::summarise(mean_val=mean(ratio)) %>% 
    na.omit() %>% 
    left_join(max.val.group, by=c("group","time")) %>% 
    left_join(min.val.group, by=c("group","time")) %>% 
    left_join(std.val.group, by=c("group","time"))
  
  base.data.guild <- thisbasedata %>%  
    left_join(species.groups, by="group") %>% 
    mutate(group_name=tolower(`Long Name`)) %>% 
    filter(group_name %in% predator.names) %>% 
    rename(guild=Guild) %>% 
    dplyr::group_by(guild, time) %>% 
    mutate(ratio=ifelse(is.nan(ratio),NA,ratio)) %>% 
    na.omit %>% 
    dplyr::summarise(base_val=mean(ratio)) %>% 
    filter(!guild %in% bad.guilds)
  
  base.data.group <- thisbasedata %>%  
    left_join(species.groups, by="group") %>% 
    mutate(group_name=tolower(`Long Name`)) %>% 
    filter(group_name %in% predator.names) %>% 
    rename(guild=Guild) %>% 
    filter(!guild %in% bad.guilds) %>% 
    dplyr::group_by(group, time) %>% 
    mutate(ratio=ifelse(is.nan(ratio),NA,ratio)) %>% 
    na.omit %>% 
    dplyr::summarise(base_val=mean(ratio))
  
  var.name <- str_to_title(thisvariable)
  
  options(digits=1)
  
  ribbon.plot.guild <- frame.data.guild %>%
    ggplot(., aes(time))+
    geom_ribbon(aes(ymin = min_val, ymax = max_val), fill = "deepskyblue3") +
    geom_line(data = base.data.guild, aes(y = base_val)) +
    geom_hline(yintercept=1, color="grey35", linetype = 2, size = 0.8)+
    facet_wrap(~ guild, scales = "free_y") +
    theme_tufte(base_size = 12, base_family='GillSans')+
    ylab(paste(var.name,"relative to no spill"))+
    xlab("Year")
  
  ribbon.plot.group <- frame.data.group %>%
    ggplot(., aes(time))+
    geom_ribbon(aes(ymin = min_val, ymax = max_val), fill = "deepskyblue3") +
    geom_line(data = base.data.group, aes(y = base_val)) +
    geom_hline(yintercept=1, color="grey35", linetype = 2, size = 0.8)+
    facet_wrap(~ group, scales = "free_y") +
    theme_tufte(base_size = 12, base_family='GillSans')+
    ylab(paste(var.name,"relative to no spill"))+
    xlab("Year")
  
  setwd(work.folder)
  ggsave(paste("simulated",thisvariable,"guild.png",sep="_"), ribbon.plot.guild, device = "png", width = 10, height = 8, dpi=300)
  
  ggsave(paste("simulated",thisvariable,"group.png",sep="_"), ribbon.plot.group, device = "png", width = 18, height = 8, dpi=300)
  
}

plot.trajectory (biomass.data, thisvariable="Biomass",thisbase=base.biomass, thisnooil = nooil.biomass)
plot.trajectory (catch.data, thisvariable="Catch",thisbase=base.catch, thisnooil = nooil.catch)

options(digits=4)

# combine emulated data and plot

fit.results <- list.files(emulator.folder, pattern="*model_preds.csv$", full.names = TRUE)

emu.preds <-function(eachresult){  
  
  this.timestep <- eachresult %>% str_split(.,"/") %>% 
    unlist %>% .[5] %>% 
    str_split(.,"_") %>% 
    unlist %>% .[1]
  
  print(this.timestep)
  
  setwd(emulator.folder)
  
  this.result <- fread(eachresult, verbose = TRUE) %>% tbl_df %>% 
    mutate(model=toupper(model)) %>% 
    mutate(time=as.factor(time)) %>% 
    rename(group=species) %>% 
    left_join(species.groups,by="group") %>% 
    mutate(group_name=tolower(`Group Description`)) %>% 
    filter(group_name %in% predator.names) %>% 
    rename(variable=var) %>% 
    rename(guild=Guild)
  
  print(this.timestep)
  
  this.min.data <- min.rmse %>% 
    filter(time==this.timestep)
  
  print(this.timestep)
  #use set with min rmse value to select only obs and pred values for that model for each species
  min.rmse.data <- this.result %>% 
    left_join(this.min.data, by=c("group","model","variable","time","guild", "GroupID","isFish",
                                  "VirgBToCurrentB", "Diag Name","Long Name","Group Description","Name")) %>% 
    na.omit() %>% 
    dplyr::select(group,time,variable,model,pred,group_name,guild)
  
  print(this.timestep)
  min.rmse.data 
  
}

emupred.data <- lapply(fit.results,emu.preds) %>% 
  rbindlist() %>% tbl_df 

target.time <- c(365,730,1095,1460,1825,2190,2555,2920,3285,3650,4015,4380,4745,5110,5474.5)
bad.guilds <- c("Epibenthos","Filter feeders","Invertebrates","Other fish")

pred.trajectory <- function(thisdata, thisvariable, thisbase, thisnooil){
  
  
  
  this.preddata <- thisdata %>% 
    filter(variable==thisvariable) %>% 
    mutate(model=as.factor(model)) %>% 
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
    mutate(simulation=rep(1:10000,times=(nrow(.)/10000))) %>% 
    mutate(pred=ifelse(pred<0,0,pred)) %>% 
    left_join(thisnooil, by=c("group","time")) %>% 
    rename(variable=variable.x, nooil = variable.y) %>% 
    mutate(ratio = pred/nooil) %>% 
    dplyr::select(time,group,ratio,everything(), -var) %>% 
    mutate(ratio=ifelse(is.nan(ratio),NA,ratio)) %>% 
    na.omit
  
  thisbasedata <- thisbase %>% left_join(thisnooil, by=c("group","time")) %>% 
    rename(variable=variable.x, nooil = variable.y, var = var.x) %>% 
    mutate(ratio = variable/nooil) %>% 
    dplyr::select(time,group,ratio,everything(), -var.y)
  
  max.val.guild <- this.preddata %>% 
    dplyr::group_by(guild, time) %>% 
    dplyr::summarise(max_val=max(ratio)) %>% 
    na.omit()
  
  min.val.guild <- this.preddata %>% 
    dplyr::group_by(guild, time) %>% 
    dplyr::summarise(min_val=min(ratio)) %>% 
    na.omit()
  
  std.val.guild <- this.preddata %>% 
    dplyr::group_by(guild, time) %>% 
    dplyr::summarise(std_val=std.error(ratio)) %>% 
    na.omit()
  
  frame.data.guild <- this.preddata %>% tbl_df %>% 
    dplyr::group_by(guild, time) %>% 
    dplyr::summarise(mean_val=mean(ratio)) %>% 
    left_join(max.val.guild, by=c("guild","time")) %>% 
    left_join(min.val.guild, by=c("guild","time")) %>% 
    left_join(std.val.guild, by=c("guild","time")) %>% 
    filter(!guild %in% bad.guilds)
  
  
  # summarize by group
  
  max.val.group <- this.preddata %>% 
    filter(!guild %in% bad.guilds) %>% 
    dplyr::group_by(group, time) %>% 
    dplyr::summarise(max_val=max(ratio)) %>% 
    na.omit()
  
  min.val.group <- this.preddata %>% 
    filter(!guild %in% bad.guilds) %>% 
    dplyr::group_by(group, time) %>% 
    dplyr::summarise(min_val=min(ratio)) %>% 
    na.omit()
  
  std.val.group <- this.preddata %>% 
    filter(!guild %in% bad.guilds) %>% 
    dplyr::group_by(group, time) %>% 
    dplyr::summarise(std_val=std.error(ratio)) %>% 
    na.omit()
  
  frame.data.group <- this.preddata %>% tbl_df %>% 
    filter(!guild %in% bad.guilds) %>% 
    dplyr::group_by(group, time) %>% 
    dplyr::summarise(mean_val=mean(ratio)) %>% 
    na.omit() %>% 
    left_join(max.val.group, by=c("group","time")) %>% 
    left_join(min.val.group, by=c("group","time")) %>% 
    left_join(std.val.group, by=c("group","time"))
  
  base.data.guild <- thisbasedata %>%  
    left_join(species.groups, by="group") %>% 
    mutate(group_name=tolower(`Long Name`)) %>% 
    filter(group_name %in% predator.names) %>% 
    rename(guild=Guild) %>% 
    filter(!guild %in% bad.guilds) %>% 
    dplyr::group_by(guild, time) %>% 
    mutate(ratio=ifelse(is.nan(ratio),NA,ratio)) %>% 
    na.omit %>% 
    dplyr::summarise(base_val=mean(ratio))
  
  base.data.group <- thisbasedata %>%  
    left_join(species.groups, by="group") %>% 
    mutate(group_name=tolower(`Long Name`)) %>% 
    filter(group_name %in% predator.names) %>% 
    rename(guild=Guild) %>% 
    filter(!guild %in% bad.guilds) %>% 
    dplyr::group_by(group, time) %>% 
    mutate(ratio=ifelse(is.nan(ratio),NA,ratio)) %>% 
    na.omit %>% 
    dplyr::summarise(base_val=mean(ratio))
  
  var.name <- str_to_title(thisvariable)
  
  ribbon.plot.guild <- frame.data.guild %>%
    ggplot(., aes(time))+
    geom_ribbon(aes(ymin = min_val, ymax = max_val), fill = "darkgoldenrod2") +
    geom_line(data = base.data.guild, aes(y = base_val)) +
    geom_hline(yintercept=1, color="grey35", linetype = 2, size = 0.8)+
    facet_wrap(~ guild, scales = "free_y") +
    theme_tufte(base_size = 12, base_family='GillSans')+
    ylab(paste(var.name,"relative to no spill"))+
    xlab("Year")
  
  ribbon.plot.group <- frame.data.group %>%
    ggplot(., aes(time))+
    geom_ribbon(aes(ymin = min_val, ymax = max_val), fill = "darkgoldenrod2") +
    geom_line(data = base.data.group, aes(y = base_val)) +
    geom_hline(yintercept=1, color="grey35", linetype = 2, size = 0.8)+
    facet_wrap(~ group, scales = "free_y") +
    theme_tufte(base_size = 12, base_family='GillSans')+
    ylab(paste(var.name,"relative to no spill"))+
    xlab("Year")
  
  setwd(work.folder)
  ggsave(paste("emulated",thisvariable,"guild.png",sep="_"), ribbon.plot.guild, device = "png", width = 10, height = 8, dpi=300)
  
  ggsave(paste("emulated",thisvariable,"group.png",sep="_"), ribbon.plot.group, device = "png", width = 18, height = 8, dpi=300)
  
  
group.max.bio <- frame.data.group %>% 
    group_by(group) %>% 
    summarise(max_group=max(max_val)) 
 
   # plot for number of samples per predator 
  samples.frame <- fread("Final_CSV_For_Dirchelet.csv", header = TRUE) %>% 
  mutate(group_name=tolower(Predator)) %>%
  left_join(species.groups, by="group_name") %>% 
  mutate(index = 1) %>% 
  filter(group_name %in% predator.names) %>% 
  group_by(group_name) %>% 
  summarise(sum_prey = sum(index)) %>% 
  left_join(species.groups, by="group_name") %>% 
  filter(!Guild %in% bad.guilds) %>% 
  left_join(group.max.bio, by = "group")
  
  var.name <- tolower(thisvariable)
  
  sample.max <- samples.frame %>%
    ggplot(., aes(sum_prey, max_group, color = Guild))+
    geom_point() +
    facet_wrap(~ Guild, scales = "free_y")+
    theme_tufte(base_size = 12, base_family='GillSans')+
    ylab(paste("Maximum", var.name,"relative to no spill"))+
    xlab("Number of samples")
    #stat_smooth(method = "lm", col = "red")
  
  ggsave(paste("group_sample",thisvariable,".png",sep="_"), sample.max, device = "png", width = 12, height = 8, dpi=300)
  
  }


pred.trajectory(thisdata=emupred.data,thisvariable="biomass",thisbase = base.biomass,thisnooil=nooil.biomass)

pred.trajectory(thisdata=emupred.data,thisvariable="catch",thisbase = base.catch,thisnooil=nooil.catch)


#set decimal output back to 4
options(digits=4)

#' Recovery time, estimate time by 90% of original biomass
#' 
setwd("~/GOMAtlantisOil_00")

bad.guilds <- c("Epibenthos","Filter feeders","Invertebrates","Other fish", "Zooplankton", "Infauna", "Primary producers",
                "Bacteria & detritus","Marine mammals", "Seabirds","Sea turtles")

virgin.biomass <- scan("log.txt",character(0), sep = "\n", skip=2362, nlines=92) %>% 
  as.data.frame() %>% setNames("text") %>% 
  separate("text", c("text1","text2"),sep="species") %>% 
  separate("text2", c("group","vir_biom"),sep="virgin biomass is") %>%
  mutate(vir_biom=gsub(" t","",vir_biom), vir_biom=as.numeric(vir_biom), group=gsub(" ","",group)) %>% 
  dplyr::select(group,vir_biom) %>% tbl_df %>% 
  left_join(species.groups, by=c("group"))
  
min.base.biomass <- base.biomass %>% 
  left_join(species.groups, by=c("group")) %>% 
  filter(!Guild %in% bad.guilds) %>% 
  left_join(virgin.biomass, by="group") %>% 
  mutate(ratio_v_b = variable/vir_biom) %>% 
  group_by(group) %>% 
  dplyr::summarise(min_biom_r=min(ratio_v_b)) 

min.time.base.biomass <- base.biomass %>% 
  left_join(species.groups, by=c("group")) %>% 
  filter(!Guild %in% bad.guilds) %>%  
  left_join(virgin.biomass, by="group") %>% 
  mutate(ratio_v_b = variable/vir_biom) %>% 
  mutate(recovered=ratio_v_b>=0.99) %>% 
  dplyr::select(recovered,everything()) %>% 
  filter(recovered==TRUE) %>% 
  group_by(recovered,group) %>% 
  dplyr::summarise(min_time=min(time))

min.base.frame <- min.time.base.biomass %>% 
  left_join(min.base.biomass, by=c("group")) %>% 
  mutate(min_months=(min_time-2010)*12) %>% 
  left_join(species.groups, by="group")%>% 
  mutate(class="simulator") %>% 
  filter(!Guild %in% bad.guilds) 
  

min.sim.biomass <- biomass.data %>% tbl_df %>% 
   filter(!guild %in% bad.guilds) %>% 
  left_join(virgin.biomass, by="group") %>% 
  mutate(ratio_v_b = variable/vir_biom) %>% 
  group_by(group,run) %>% 
  dplyr::summarise(min_biom_r=min(ratio_v_b)) 

min.sim.bio.sum <- min.sim.biomass %>% 
  group_by(group) %>% 
  summarise(min_biom_g=min(min_biom_r),max_biom_g=max(min_biom_r),mean_biom_g=max(min_biom_r))
  
min.time.sim.biomass <- biomass.data %>% tbl_df %>% 
  filter(!guild %in% bad.guilds) %>% 
  left_join(virgin.biomass, by="group") %>% 
  mutate(ratio_v_b = variable/vir_biom) %>% 
  mutate(recovered=ratio_v_b>=0.99) %>% 
  dplyr::select(recovered,everything()) %>% 
  filter(recovered==TRUE) %>% 
  group_by(recovered,group,run) %>% 
  dplyr::summarise(min_time=min(time))

min.sim.time.sum <- min.time.sim.biomass %>% 
  group_by(group) %>% 
  summarise(min_time_g=min(min_time),max_time_g=max(min_time),mean_time_g=mean(min_time))

min.biom.frame <- min.time.sim.biomass %>% 
  left_join(min.sim.biomass, by=c("group","run")) %>% 
  mutate(min_months=(min_time-2010)*12) %>% 
  left_join(species.groups, by="group")%>% 
  mutate(class="simulator")

sum.biom.frame <- min.sim.time.sum %>% 
  left_join(min.sim.bio.sum, by=c("group")) %>% 
  mutate(min_months=(min_time_g-2010)*12) %>%
  mutate(max_months=(max_time_g-2010)*12) %>%
  mutate(mean_months=(mean_time_g-2010)*12) %>% 
  left_join(species.groups, by="group")

not.recovered <- anti_join(min.sim.bio.sum,min.sim.time.sum, by="group") %>% 
  left_join(species.groups, by="group") %>% 
  dplyr::select(Guild, everything()) %>% 
  distinct(Guild, group)
print("groups not recovered")
print(not.recovered)
          
min.emu.biomass <- emupred.data %>% 
  filter(variable=="biomass") %>% 
  filter(!guild %in% bad.guilds) %>% 
  mutate(model=as.factor(model)) %>% 
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
  mutate(run=rep(1001:11000,times=(nrow(.)/10000))) %>% 
  mutate(pred=ifelse(pred<0,0,pred))  %>% 
  tbl_df %>%
  left_join(virgin.biomass, by="group") %>% 
  mutate(ratio_v_b = pred/vir_biom) %>% 
  group_by(group,run) %>% 
  dplyr::summarise(min_biom_r=min(ratio_v_b)) 

min.time.emu.biomass <- emupred.data %>% 
  filter(variable=="biomass") %>% 
  filter(!guild %in% bad.guilds) %>% 
  mutate(model=as.factor(model)) %>% 
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
  mutate(run=rep(1001:11000,times=(nrow(.)/10000))) %>% 
  mutate(pred=ifelse(pred<0,0,pred))  %>% 
  tbl_df %>% 
  left_join(virgin.biomass, by="group") %>% 
  mutate(ratio_v_b = pred/vir_biom) %>% 
  mutate(recovered=ratio_v_b>=0.99) %>% 
  dplyr::select(recovered,everything()) %>% 
  filter(recovered==TRUE) %>% 
  group_by(recovered,group,run) %>% 
  dplyr::summarise(min_time=min(time))

min.emu.frame <- min.time.emu.biomass %>% 
  left_join(min.emu.biomass, by=c("group","run")) %>% 
  mutate(min_months=(min_time-2010)*12) %>% 
  left_join(species.groups, by="group") %>% 
  mutate(class="emulator")

result.frame <- bind_rows(min.biom.frame,min.emu.frame)

pal.sim <- c("cadetblue3","chocolate2","darkolivegreen3","turquoise3","darkgoldenrod1","springgreen3","royalblue2","indianred3","darkseagreen3")

sim.recovery.plot <- result.frame %>%
  ggplot(., aes(min_months, min_biom_r,color=Guild,fill = Guild))+ #, shape = Guild
  geom_point(size=1.5)+
  geom_point(aes(min_months, min_biom_r), data = min.base.frame, color = "black", shape = 17, size = 1.75)+
  #geom_jitter(width = 0.25, height = 0.25)+
  facet_wrap(~Guild)+ 
  xlim(1, 190)+
  scale_color_manual(values=pal.sim)+
  #scale_shape_manual(values=c(9,1,5,7,6,3,4,8,2))+
  geom_hline(yintercept=1, color="royalblue4", size = 0.65)+
  geom_vline(xintercept=50, color="royalblue4", linetype = 3, size = 0.65)+
  theme_tufte(base_size = 14, base_family='GillSans')+
  #  theme_classic(base_size = 12, base_family='GillSans')+
  ylab("Biomass minima")+
  xlab("Months to recovery")+
  theme(legend.position="none")

pal.hist <- c("darkslategray4","darkgoldenrod1")
sim.hist.plot <- result.frame %>%
  ggplot(., aes(min_months, color = class, fill = class))+
  geom_bar(stat="count")+
  scale_color_manual(values=pal.hist, guide = guide_legend(title = "Model"))+
  scale_fill_manual(values=pal.hist, guide = guide_legend(title = "Model"))+
  xlim(1, 190)+
  theme_classic(base_size = 14, base_family='GillSans')+
  ylab(paste("Count"))+
  xlab("")
  

recovery.plot <- ggarrange(sim.hist.plot, sim.recovery.plot, 
                       ncol = 1, nrow = 2, heights = c(0.85,2.15), align = "v") #widths = c(1,2), 

setwd(work.folder)
ggsave(
  "recovery_plot.png",
  recovery.plot,
  device = "png",
  height = 10,
  width = 14,
  dpi = 400
)


  
  

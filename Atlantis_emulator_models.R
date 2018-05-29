#' Code to build emulator based on Atlantis runs for
#' GOM sensitivity analysis
#' @author Hem Nalini Morzaria Luna
#' @date June 2017

# List of packages for session
.packages = c("readr",
              "stringi",
              "data.table",
              "tidyverse",
              "stringr",
              "R.utils")

# Install CRAN packages (if not already installed)
.inst <- .packages %in% installed.packages()
if (length(.packages[!.inst]) > 0)
  install.packages(.packages[!.inst], dependencies = TRUE)

# Load packages into session
lapply(.packages, require, character.only = TRUE)

# clean up the space
rm(list = ls())

result.folders <- "~/biomassres"
work.folder <- "~/GOMcode"
data.folder <- "~/emulatorprms"
save.folder <- "~/modelres"
dir.create(save.folder)
dir.create(data.folder)

#result.folders <- "G:/Climate_data_2017/GOMAtlantis_Sensitivity_results/Results"
#work.folder <- "E:/Archivos/1Archivos/Articulos/En preparacion/GOM/Code"
#prm.folder <- "E:/Archivos/1Archivos/Articulos/En preparacion/GOM/Code/prms"

setwd(work.folder)

#' Machine learning emulator models
setwd(work.folder)
target.time <- 365

#c(365,730,1095,1460,1825,2190,2555,2920,3285,3650,4015,4380,4745,5110,5474.5)

this.data <- paste(target.time, "_emulator_data.csv", sep = "")

all.biomass <- fread(this.data, header = TRUE, verbose = TRUE) %>%
  tbl_df

# These are all functional groups in the model
func.names <- all.biomass %>%
  distinct(group) %>%
  .$group %>%
  as.character()

# Use this to create individual files for each species:time combination
for (eachspecies in func.names) {
  this.species.data <- all.biomass %>% filter(group == eachspecies)
  print(paste("Analyzing species", eachspecies))
  
  write_csv(this.species.data, file.path(
    data.folder,
    paste(eachspecies, target.time, "emulator_data.csv", sep = "_")
  ))
}


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
data.folder <- "~/emulatorprms"
save.folder <- "~/modelres"
dir.create(save.folder)
dir.create(data.folder)

setwd(work.folder)

# species names and guilds
species.names <-
  fread("Group_Names.csv", header = TRUE) %>% tbl_df %>%
  mutate(`predator name` = tolower(`Long Name`)) %>%
  dplyr::select(`predator name`, Code, Guild) %>%
  setNames(c("group_name", "group", "guild"))

#get names of predators that were randomized and use those as response variables in emulator
predator.names <- fread("predator_names.csv", header = TRUE) %>%
  dplyr::select(predator) %>% .$predator %>% tolower()

setwd(work.folder)
target.time <- 365
#possible timesteps, for each year. Original output is 30 day intervals
#c(365,730,1095,1460,1825,2190,2555,2920,3285,3650,4015,4380,4745,5110,5474.5)

#list individual catch and biomass files
list.data.files <-
  list.files(data.folder, pattern = "emulator", full.names = TRUE)

list.index <- 1:length(list.data.files)

nThreads <- detectCores() - 1

# Initiate cluster

cl <- makeSOCKcluster(nThreads)
registerDoSNOW(cl)
clusterCall(cl, function(x)
  .libPaths(x), .libPaths())

#' Machine learning emulator models

glm.model <- function(eachspecies) {
  #read datafile

  #get species name
  group.name <-
    str_split(list.data.files[eachspecies], pattern = "/") %>%
    unlist() %>% .[5] %>%
    str_split(., pattern = "_") %>%
    unlist() %>% .[1]
  
  print(paste("Analyzing this species", group.name, eachspecies))
  
  #create training and testing partitions
  # set seed to make repeatable
  
  set.seed(235)
  
  this.dat <- future({
    fread(list.data.files[eachspecies]) %>%
      tbl_df %>%
      dplyr::select(-group)
    
  }) %plan% multiprocess
  
  this.data <- value(this.dat)
  
  
  inTraining <- 
    createDataPartition(this.data$run, p = .75, list = FALSE)
  #train data
  this.train.data <- this.data[inTraining, ] %>%
    dplyr::select(-Time, -run)
  #test data
  this.test.data  <- this.data[-inTraining,]
  
  #' Analyze catch data
  print(paste("Analyzing catch data", group.name))
  
  #remove biomass
  this.train.data.v <- this.train.data %>%
    dplyr::select(-biomass)
  
  #test that species has catch
  catch.test <- this.train.data.v %>%
    dplyr::select(catch) %>%
    summarise(sum_c = sum(catch)) %>%
    .$sum_c
  
  biomass.test <- this.train.data %>%
    dplyr::select(biomass) %>%
    distinct(biomass) %>%
    NROW()
  
  if (biomass.test > 1) {
    #if dataset has catch then build emulators
    if (catch.test != 0) {
      #specify variable
      this.variable <- "catch"
      
      # fit Generalized Linear Model using package caret
      #specify model
      this.model <- "glm"
      
      print(paste("Analyzing ", this.model, this.variable, group.name))
      #set seed for repeatability
      #fit model
      f.glm <- future({
       
        set.seed(235)
        caret::train(catch ~ .,
                              data = this.train.data.v,
                              method = 'glm')
    }) %plan% multiprocess

      fit.glm <- value(f.glm)
      
      #summary of model fit
      summary.fit <- try(summary(fit.glm))
      
      #this separates variables that are singularities, and are therefore not independent variables
      
      coeff.frame <- summary.fit$coefficients %>%
        as.data.frame() %>% rownames() %>% .[-1]
      
      this.train.b.sel <- this.train.data.v %>%
        dplyr::select(catch, one_of(coeff.frame)) %>%
        as.data.frame()
      
      # use K-fold cross-validation to prevent overfitting
      # fit generalized linear model
      # removes availability variables that are not randomized
      
      f.glm.sel <- future({
      fitControl <- trainControl(
        ## 10-fold CV
        method = "repeatedcv",
        number = 10,
        ## repeated ten times
        repeats = 10,
        savePredictions = "final"
      )
      #set seed for repeatability
      set.seed(235)
      
      print(paste("Analyzing ", this.model, this.variable, group.name))
      # fit model
      caret::train(
        catch ~ .,
        data = this.train.b.sel,
        method = 'glm',
        trControl = fitControl,
        tuneGrid = expand.grid(.parameter = sqrt(ncol(
          this.train.b.sel
        )))
      )
      
      }) %plan% multiprocess
      
      fit.glm.sel <- value(f.glm.sel)
      
      #summary model fit
      summary.fit.glm <- try(summary(fit.glm.sel))
      
      #extract variable importance
      lmVarImp  <- varImp(fit.glm.sel)
      imp.var.glm <- lmVarImp$importance %>%
        tbl_df %>%
        mutate(var = rownames(lmVarImp$importance))
      
      #save model output
      setwd(save.folder)
      capture.output(
        summary.fit.glm ,
        file = paste(
          target.time,
          group.name,
          this.model,
          this.variable,
          "fit.txt",
          sep = "_"
        )
      )
      
      save(
        fit.glm.sel,
        file = paste(
          target.time,
          group.name,
          this.model,
          this.variable,
          "model.rda",
          sep = "_"
        )
      )
      write.csv(
        imp.var.glm,
        file = paste(
          target.time,
          group.name,
          this.model,
          this.variable,
          "varimp.csv",
          sep = "_"
        )
      )
      
    }
  }
  
  
  
    #fit models using futures
  
  #' fit Extreme Learning Machine
  
  #specify model name
  this.model <- "elm"
  print(paste("Analyzing ", this.model, this.variable, group.name))
  
  #fit model
  fit.elm.sel <- caret::train(
    catch ~ .,
    data = this.train.b.sel,
    method = "elm",
    tuneGrid = expand.grid(
      nhid = 5:50,
      actfun = c('sig', 'sin', 'purelin', 'radbas', 'poslin')
    ),
    trControl = fitControl,
    metric = "RMSE",
    type = 'Regression'
  )
  
  #summary model fit
  summary.fit.elm <- try(summary(fit.elm.sel))
  #extract variable importance
  elmVarImp  <- varImp(fit.elm.sel)
  imp.var.elm <- elmVarImp$importance %>%
    tbl_df %>%
    mutate(var = rownames(elmVarImp$importance))
  #save model output
  setwd(save.folder)
  capture.output(
    summary.fit.elm ,
    file = paste(
      target.time,
      group.name,
      this.model,
      this.variable,
      "fit.txt",
      sep = "_"
    )
  )
  save(
    fit.elm.sel,
    file = paste(
      target.time,
      group.name,
      this.model,
      this.variable,
      "model.rda",
      sep = "_"
    )
  )
  write.csv(
    imp.var.elm,
    file = paste(
      target.time,
      group.name,
      this.model,
      this.variable,
      "varimp.csv",
      sep = "_"
    )
  )
  
  #' Fit bagged CART
  
  #specify model
  this.model <- "bcart"
  print(paste("Analyzing ", this.model, this.variable, group.name))
  #set seed for repeatability
  set.seed(235)
  #fit model
  fit.bcart.sel <- caret::train(
    catch ~ .,
    data = this.train.b.sel,
    method = 'treebag',
    trControl = fitControl,
    metric = 'RMSE'
  )
  
  #summary model fit
  summary.fit.bcart <- try(summary(fit.bcart.sel))
  
  #extract variable importance
  bamVarImp  <- varImp(fit.bcart.sel)
  imp.var.bcart <- bamVarImp$importance %>%
    tbl_df %>%
    mutate(var = rownames(bamVarImp$importance))
  
  setwd(save.folder)
  capture.output(
    summary.fit.bcart ,
    file = paste(
      target.time,
      group.name,
      this.model,
      this.variable,
      "fit.txt",
      sep = "_"
    )
  )
  save(
    fit.bcart.sel,
    file = paste(
      target.time,
      group.name,
      this.model,
      this.variable,
      "model.rda",
      sep = "_"
    )
  )
  write.csv(
    imp.var.bcart,
    file = paste(
      target.time,
      group.name,
      this.model,
      this.variable,
      "varimp.csv",
      sep = "_"
    )
  )
  
  #' fit random forest
  #specify model
  this.model <- "rf"
  print(paste("Analyzing ", this.model, this.variable, group.name))
  
  #set seed for repeatability
  set.seed(235)
  
  #fit model
  fit.rf.sel <- caret::train(
    catch ~ .,
    data = this.train.b.sel,
    method = "rf",
    ntree = 100,
    importance = TRUE,
    tuneGrid = expand.grid(.mtry = (sqrt(
      ncol(this.train.b.sel)
    ))),
    trControl = fitControl
  )
  #summary model fit
  summary.fit.rf <- try(summary(fit.rf.sel))
  
  #extract variable importance
  rfVarImp  <- varImp(fit.rf.sel)
  imp.var.rf <- rfVarImp$importance %>%
    tbl_df %>%
    mutate(var = rownames(rfVarImp$importance))
  
  #save model outputs
  setwd(save.folder)
  capture.output(
    summary.fit.rf ,
    file = paste(
      target.time,
      group.name,
      this.model,
      this.variable,
      "fit.txt",
      sep = "_"
    )
  )
  save(
    fit.rf.sel,
    file = paste(
      target.time,
      group.name,
      this.model,
      this.variable,
      "model.rda",
      sep = "_"
    )
  )
  write.csv(
    imp.var.rf,
    file = paste(
      target.time,
      group.name,
      this.model,
      this.variable,
      "varimp.csv",
      sep = "_"
    )
  )
  
  #' fit boosted regression tree
  
  #specify model
  this.model <- "brt"
  print(paste("Analyzing ", this.model, this.variable, group.name))
  
  # set seed for repeatability
  set.seed(235)
  
  #get model fit
  fit.brt.sel <- caret::train(
    catch ~ .,
    data = this.train.b.sel,
    method = 'bstTree',
    tuneGrid = expand.grid(
      mstop = c(1, 3, 10, 15, 50),
      maxdepth = c(2, 4, 8, 12),
      nu = c(0.01, 0.1, 0.3)
    ),
    trControl = fitControl
  )
  
  #  preProc=c('center','scale'),
  #summary model fit
  summary.fit.brt <- try(summary(fit.brt.sel))
  
  #extract variable importance
  brtVarImp  <- varImp(fit.brt.sel)
  imp.var.brt <- brtVarImp$importance %>%
    tbl_df %>%
    mutate(var = rownames(brtVarImp$importance))
  
  #save model output
  setwd(save.folder)
  capture.output(
    summary.fit.brt ,
    file = paste(
      target.time,
      group.name,
      this.model,
      this.variable,
      "fit.txt",
      sep = "_"
    )
  )
  save(
    fit.brt.sel,
    file = paste(
      target.time,
      group.name,
      this.model,
      this.variable,
      "model.rda",
      sep = "_"
    )
  )
  write.csv(
    imp.var.brt,
    file = paste(
      target.time,
      group.name,
      this.model,
      this.variable,
      "varimp.csv",
      sep = "_"
    )
  )
  
  #'bagged MARS fitting
  
  #specify model
  this.model <- "bmars"
  print(paste("Analyzing ", this.model, this.variable, group.name))
  
  #set seed for repeatability
  set.seed(235)
  
  #fit model
  fit.bmars.sel <- caret::train(
    catch ~ .,
    data = this.train.b.sel,
    method = 'bagEarth',
    tuneGrid = expand.grid(degree = c(1, 2), nprune = seq(10, 90, 20)),
    trControl = fitControl,
    maximize = FALSE,
    keepX = FALSE,
    na.action = na.omit
  )
  
  #save summary model fit
  summary.fit.bmars <- try(summary(fit.bmars.sel))
  
  #extract variable importance
  bmarsVarImp  <- varImp(fit.bmars.sel)
  imp.var.bmars <- bmarsVarImp$importance %>%
    tbl_df %>%
    mutate(var = rownames(bmarsVarImp$importance))
  
  #save model ouputs
  setwd(save.folder)
  capture.output(
    summary.fit.bmars ,
    file = paste(
      target.time,
      group.name,
      this.model,
      this.variable,
      "fit.txt",
      sep = "_"
    )
  )
  save(
    fit.bmars.sel,
    file = paste(
      target.time,
      group.name,
      this.model,
      this.variable,
      "model.rda",
      sep = "_"
    )
  )
  write.csv(
    imp.var.bmars,
    file = paste(
      target.time,
      group.name,
      this.model,
      this.variable,
      "varimp.csv",
      sep = "_"
    )
  )
  
  # collect resamples
  results <-
    resamples(
      list(
        GeneralizedLinearModel = fit.glm.sel,
        BoostedGeneralizedAdditiveModel = fit.elm.sel,
        BaggedCART = fit.bcart.sel,
        RandomForest = fit.rf.sel,
        BoostedRegressionTrees = fit.brt.sel,
        BaggedMultivariateAdaptiveRegressionSplines = fit.bmars.sel
      )
    )
  #create dataframe with model fits
  results.frame <- results$values %>%
    tbl_df %>%
    mutate(group = group.name, time = target.time) %>%
    gather(model, values, -Resample,-group,-time) %>%
    separate(model, c("model_fit", "metric"), sep = "~") %>%
    left_join(species.names, by = "group")
  
  #save fit results
  setwd(save.folder)
  write.csv(
    results.frame,
    file = paste(target.time, group.name, "catch_fitresults.csv", sep = "_")
  )
  
  #predict using models on test data
  print(paste("Creating model ensemble ", this.variable, group.name))
  #model ensemble
  model.ensemble <-  list(
    glm = fit.glm.sel,
    elm = fit.elm.sel,
    bcart = fit.bcart.sel,
    rf = fit.rf.sel,
    brt = fit.brt.sel,
    bmars = fit.bmars.sel
  )
  #subset test data leaving only predictors
  this.test.data.sel <- this.test.data %>%
    dplyr::select(-Time,-run,-catch, -biomass)
  
  #' predict new values using test data
  print(paste("Predicting test data ", group.name))
  #predict on test data
  preds.log <-
    predict(model.ensemble, this.test.data.sel) %>%
    bind_rows() %>%
    tbl_df %>%
    mutate(obs = this.test.data$catch) %>%
    gather(model, pred, -obs) %>%
    mutate(
      group = group.name,
      residual = obs - pred,
      time = target.time,
      variable = this.variable
    ) %>%
    left_join(species.names, by = "group") %>%
    dplyr::select(time,
                  variable,
                  group,
                  group_name,
                  guild,
                  model,
                  obs,
                  pred,
                  residual)
  
  #write predictions
  write.csv(preds.log,
            file = paste(target.time, group.name, this.variable, "pred.csv", sep = "_"))
  
}

#' Analyze biomass data
#remove catch from data
this.variable <- "biomass"
this.train.data.v <- this.train.data %>%
  dplyr::select(-catch)

# fit generalized linear model using package caret

this.model <- "glm"

print(paste("Analyzing ", this.model, this.variable, group.name))

f.glm <- future({
  caret::train(biomass ~ .,
               data = this.train.data.v,
               method = 'glm')
}) %plan% multiprocess
fit.glm <- value(f.glm)

summary.fit <- try(summary(fit.glm))

#separates variables that are singularities, and are therefore not independent variables
coeff.frame <- summary.fit$coefficients %>%
  as.data.frame() %>% rownames() %>% .[-1]

this.train.b.sel <- this.train.data.v %>%
  dplyr::select(biomass, one_of(coeff.frame)) %>%
  as.data.frame() %>%
  mutate_all(funs(as.numeric))

# use K-fold cross-validation to prevent overfitting
# fit generalized linear model
# removes availability variables that are not randomized
fitControl <- trainControl(
  ## 10-fold CV
  method = "repeatedcv",
  number = 10,
  ## repeated ten times
  repeats = 10,
  savePredictions = "final"
)

#set seed to ensure repeatability
set.seed(235)

print(paste("Analyzing ", this.model, group.name))

#fit model

f.glm.sel <- future({
  caret::train(
    biomass ~ .,
    data = this.train.b.sel,
    method = 'glm',
    trControl = fitControl,
    tuneGrid = expand.grid(.parameter = sqrt(ncol(
      this.train.b.sel
    )))
  )
}) %plan% multiprocess
fit.glm.sel <- value(f.glm.sel)

#get summary of model fit
summary.fit.glm <- try(summary(fit.glm.sel))

#get variable importance
lmVarImp  <- varImp(fit.glm.sel)
imp.var.glm <- lmVarImp$importance %>%
  tbl_df %>%
  mutate(var = rownames(lmVarImp$importance))

#save model ouput
setwd(save.folder)
capture.output(
  summary.fit.glm ,
  file = paste(
    target.time,
    group.name,
    this.model,
    this.variable,
    "fit.txt",
    sep = "_"
  )
)
save(
  fit.glm.sel,
  file = paste(
    target.time,
    group.name,
    this.model,
    this.variable,
    "model.rda",
    sep = "_"
  )
)
write.csv(
  imp.var.glm,
  file = paste(
    target.time,
    group.name,
    this.model,
    this.variable,
    "varimp.csv",
    sep = "_"
  )
)


}

}  
    #set seed so results are repeatable
    
    #fit model
    f.bmars.sel <- future({
      set.seed(235)
      
      try(caret::train(
        biomass ~ .,
        data = this.train.b.sel,
        method = 'bagEarth',
        tuneGrid = expand.grid(degree = c(2), nprune = seq(10, 90, 20)),
        trControl = fitControl,
        maximize = FALSE,
        keepX = FALSE,
        na.action = na.omit,
        metric = "RMSE"
      ))
    }) %plan% multiprocess
    fit.bmars.sel  <- value(f.bmars.sel)
    
   while (!resolved(fit.bmars.sel)) {
      
      f.elm.sel <- future({
        try(caret::train(
          biomass ~ .,
          data = this.train.b.sel,
          method = "elm",
          tuneGrid = expand.grid(
            nhid = 5:50,
            actfun = c('sig', 'sin', 'purelin', 'radbas', 'poslin')
          ),
          trControl = fitControl,
          metric = "RMSE",
          type = 'Regression'
        ))
      }) %plan% multiprocess
      
    }
  
    
    while (!resolved(fit.bmars.sel)) {
      
    
    }
    
    #' fit Extreme Learning Machine
    #specify model
    this.model <- "elm"
    print(paste("Analyzing ", this.model, this.variable, group.name))
    

 
    
    
    summary.fit.elm <- try(summary(fit.elm.sel))
    #extract variable importance
    elmVarImp  <- varImp(fit.elm.sel)
    imp.var.elm <- elmVarImp$importance %>%
      tbl_df %>%
      mutate(var = rownames(elmVarImp$importance))
    
    #save model ouputs
    setwd(save.folder)
    capture.output(
      summary.fit.elm ,
      file = paste(
        target.time,
        group.name,
        this.model,
        this.variable,
        "fit.txt",
        sep = "_"
      )
    )
    save(
      fit.elm.sel,
      file = paste(
        target.time,
        group.name,
        this.model,
        this.variable,
        "model.rda",
        sep = "_"
      )
    )
    write.csv(
      imp.var.elm,
      file = paste(
        target.time,
        group.name,
        this.model,
        this.variable,
        "varimp.csv",
        sep = "_"
      )
    )
    
    #' Fit bagged CART
    #specify model
    this.model <- "bcart"
    print(paste("Analyzing ", this.model, this.variable, group.name))
    
    f.bcart.sel <- future({
      caret::train(
        biomass ~ .,
        data = this.train.b.sel,
        method = 'treebag',
        trControl = fitControl,
        metric = 'RMSE'
      )
    }) %plan% multiprocess
    fit.bcart.sel  <- value(f.bcart.sel)
    
    #get model summary
    summary.fit.bcart <- try(summary(fit.bcart.sel))
    #extract variable importance
    bamVarImp  <- varImp(fit.bcart.sel)
    imp.var.bcart <- bamVarImp$importance %>%
      tbl_df %>%
      mutate(var = rownames(bamVarImp$importance))
    
    #save output
    setwd(save.folder)
    capture.output(
      summary.fit.bcart ,
      file = paste(
        target.time,
        group.name,
        this.model,
        this.variable,
        "fit.txt",
        sep = "_"
      )
    )
    save(
      fit.bcart.sel,
      file = paste(
        target.time,
        group.name,
        this.model,
        this.variable,
        "model.rda",
        sep = "_"
      )
    )
    write.csv(
      imp.var.glm,
      file = paste(
        target.time,
        group.name,
        this.model,
        this.variable,
        "varimp.csv",
        sep = "_"
      )
    )
    
    #' fit random forest
    this.model <- "rf"
    print(paste("Analyzing ", this.model, group.name))
    
     
    #fit model
    f.rf.sel <- future({
      set.seed(235)
      
      caret::train(
          biomass ~ .,
          data = this.train.b.sel,
          method = "rf",
          ntree = 100,
          importance = TRUE,
          tuneGrid = expand.grid(.mtry = (sqrt(
            ncol(this.train.b.sel)
          ))),
          trControl = fitControl
        )
    }) %plan% multiprocess
    fit.rf.sel  <- value(f.rf.sel)
    
    #get summary of model
    summary.fit.rf <- try(summary(fit.rf.sel))
    #extract variable importance
    rfVarImp  <- varImp(fit.rf.sel)
    imp.var.rf <- rfVarImp$importance %>%
      tbl_df %>%
      mutate(var = rownames(rfVarImp$importance))
    
    #save outputs
    setwd(save.folder)
    capture.output(
      summary.fit.rf,
      file = paste(
        target.time,
        group.name,
        this.model,
        this.variable,
        "fit.txt",
        sep = "_"
      )
    )
    save(
      fit.rf.sel,
      file = paste(
        target.time,
        group.name,
        this.model,
        this.variable,
        "model.rda",
        sep = "_"
      )
    )
    write.csv(
      imp.var.rf,
      file = paste(
        target.time,
        group.name,
        this.model,
        this.variable,
        "varimp.csv",
        sep = "_"
      )
    )
    
    #' fit boosted regression tree
    this.model <- "brt"
    print(paste("Analyzing ", this.model, this.variable, group.name))
    
    #fit model
    
    
    f.brt.sel <- future({
      #set seed so results are repeatable
      set.seed(235)
      #create tuning grid, use getModelInfo() to see parameters
      tunegrid <- expand.grid(
        mstop = c(1, 3, 10, 15, 50),
        maxdepth = c(2, 4, 8, 12),
        nu = c(0.01, 0.1, 0.3)
      )
      
      caret::train(
        biomass ~ .,
        data = this.train.b.sel,
        method = 'bstTree',
        tuneGrid = tunegrid,
        trControl = fitControl
      )
    }) %plan% multiprocess
    fit.brt.sel  <- value(f.brt.sel)
    
   
    #  preProc=c('center','scale'),
    
    #get summary of data
    summary.fit.brt <- try(summary(fit.brt.sel))
    #extract variable importance
    brtVarImp  <- varImp(fit.brt.sel)
    imp.var.brt <- brtVarImp$importance %>%
      tbl_df %>%
      mutate(var = rownames(brtVarImp$importance))
    #save outputs
    setwd(save.folder)
    capture.output(
      summary.fit.brt ,
      file = paste(
        target.time,
        group.name,
        this.model,
        this.variable,
        "fit.txt",
        sep = "_"
      )
    )
    save(
      fit.brt.sel,
      file = paste(
        target.time,
        group.name,
        this.model,
        this.variable,
        "model.rda",
        sep = "_"
      )
    )
    write.csv(
      imp.var.brt,
      file = paste(
        target.time,
        group.name,
        this.model,
        this.variable,
        "varimp.csv",
        sep = "_"
      )
    )
    
    
    ## bagged MARS fitting
    
    #specify model
    this.model <- "bmars"
    print(paste("Analyzing ", this.model, this.variable, group.name))
    
    #save summary fit
    summary.fit.bmars <- try(summary(fit.bmars.sel))
    #extract variable importance
    bmarsVarImp  <- varImp(fit.bmars.sel)
    imp.var.bmars <- bmarsVarImp$importance %>%
      tbl_df %>%
      mutate(var = rownames(bmarsVarImp$importance))
    
    #write outputs
    setwd(save.folder)
    capture.output(
      summary.fit.bmars,
      file = paste(
        target.time,
        group.name,
        this.model,
        this.variable,
        "fit.txt",
        sep = "_"
      )
    )
    save(
      fit.bmars.sel,
      file = paste(
        target.time,
        group.name,
        this.model,
        this.variable,
        "model.rda",
        sep = "_"
      )
    )
    write.csv(
      imp.var.bmars,
      file = paste(
        target.time,
        group.name,
        this.model,
        this.variable,
        "varimp.csv",
        sep = "_"
      )
    )
    
    #collect resamples
    results <-
      resamples(
        list(
          GeneralizedLinearModel = fit.glm.sel,
          BoostedGeneralizedAdditiveModel = fit.elm.sel,
          BaggedCART = fit.bcart.sel,
          RandomForest = fit.rf.sel,
          BoostedRegressionTrees = fit.brt.sel,
          BaggedMultivariateAdaptiveRegressionSplines = fit.bmars.sel
        )
      )
    results.frame <- results$values %>%
      tbl_df %>%
      mutate(group = group.name, time = target.time) %>%
      gather(model, values, -Resample,-group,-time) %>%
      separate(model, c("model_fit", "metric"), sep = "~") %>%
      left_join(species.names, by = "group")
    
    setwd(save.folder)
    write.csv(
      results.frame,
      file = paste(target.time, group.name, "biomass_fitresults.csv", sep = "_")
    )
    
    #predict using models on test data
    print(paste("Creating model ensemble", group.name))
    
    #create model ensemble
    model.ensemble <-  list(
      glm = fit.glm.sel,
      elm = fit.elm.sel,
      bcart = fit.bcart.sel,
      rf = fit.rf.sel,
      brt = fit.brt.sel,
      bmars = fit.bmars.sel
    )
    #format test data, leave only predictors
    this.test.data.sel <- this.test.data %>%
      dplyr::select(-Time,-run,-catch, -biomass)
    
    #predict new values using test data
    print(paste("Predicting test data ", group.name))
    
    preds.log <-
      predict(model.ensemble, this.test.data.sel) %>%
      bind_rows() %>%
      tbl_df %>%
      mutate(obs = this.test.data$biomass) %>%
      gather(model, pred, -obs) %>%
      mutate(
        group = group.name,
        residual = obs - pred,
        time = target.time,
        variable = this.variable
      ) %>%
      left_join(species.names, by = "group") %>%
      dplyr::select(time,
                    variable,
                    group,
                    group_name,
                    guild,
                    model,
                    obs,
                    pred,
                    residual)
    
    #write predictions
    write.csv(preds.log,
              file = paste(target.time,
                           group.name,
                           this.variable,
                           "pred.csv",
                           sep = "_"))
    
    save(
      model.ensemble,
      file = paste(
        this.variable,
        target.time,
        group.name,
        "model_ensemble.rda",
        sep = "_"
      ))
    
    
  }
}

#close clusters
stopCluster(cl)
registerDoSEQ()

#copy results to Azure folder
for (eachfile in list.files(path = save.folder,
                            pattern = paste(target.time, group.name, sep = "_"))) {
  system(
    paste(
      "azure storage blob upload -a morzariadatasc -k goyzslE+pHjw9t9/r7JSVuP81sqZLKE0yJEKAojmw0IFas8iJvkk2bwRSxZpCKDhLjUtuhjqYWUU0cBFOv557A== ",
      save.folder,
      "/",
      eachfile,
      "modelres",
      " --quiet",
      sep = ""
    ),
    wait = TRUE
  )
  
}

# fit.rf.sel <- caret::train(
#   biomass ~ .,
#   data = this.train.data.v,
#   method = "rf",
#   ntree = 500,
#   tunelength = 10,
#   importance = TRUE,
#   # tuneGrid = expand.grid(.mtry = (sqrt(
#   #    ncol(this.train.b.sel)
#   # ))),
#   trControl = fitControl
# )

#' #' fit Boosted Generalized Additive Model
#' this.model <- "bgam"
#' print(paste("Analyzing ",this.model,group.name))
#' #' 
#' fit.bgam.sel <- caret::train(
#'   catch ~ .,
#'   data = this.train.b.sel,
#'   method = 'gamboost',
#'   trControl = fitControl,
#'   tuneGrid = expand.grid(.mstop = seq(100, 1000, 100), .prune = c(5)),
#'   metric = 'Rsquared',
#'   maximize = FALSE
#' )
#' 
#' summary.fit.bgam <- summary(fit.bgam.sel)
#' bgamVarImp  <- varImp(fit.bgam.sel)
#' imp.var.bgam <- bgamVarImp$importance %>%
#'   tbl_df %>%
#'   mutate(var = rownames(bgamVarImp$importance))
#' setwd(save.folder)
#' capture.output(
#'   summary.fit.bgam ,
#'   file = paste(target.time, group.name, this.model, "catch_fit.txt", sep = "_")
#' )
#' save(
#'   fit.bgam.sel,
#'   file = paste(
#'     target.time,
#'     group.name,
#'     this.model,
#'     "catch_model.rda",
#'     sep = "_"
#'   )
#' )
#' write.csv(
#'   imp.var.bgam,
#'   file = paste(
#'     target.time,
#'     group.name,
#'     this.model,
#'     "catch_varimp.csv",
#'     sep = "_"
#'   )
#' )
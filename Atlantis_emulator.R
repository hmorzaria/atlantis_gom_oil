#' Code to build emulator based on Atlantis runs for
#' GOM sensitivity analysis
#' @author Hem Nalini Morzaria Luna
#' @date June 2017

# List of packages for session
.packages = c("readr","stringi","data.table",
              "tidyverse","stringr","R.utils","magrittr","future",
              "parallel","doSNOW","dismo","gbm","dplyr","sparklyr")

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

new.names <- c("run","group","catch","biomass","benthicFeedingSharks","bioerodingFish", "bivalves",
               "blackDrum", "blacktipShark","blueCrab", "blueMarlin", "bluefinTuna", "deepSerranidae",
               "deepWaterFish","filterFeedingSharks", "flatfish",
               "gagGrouper","greaterAmberjack","jacks","kingMackerel",
               "ladyfish","largePelagicFish","largeReefFish","largeSharks",
               "littleTunny","lutjanidae", "mediumPelagicFish", "menhaden", 
               "mullets", "otherBillfish","otherDemersalFish","otherShrimp", 
               "otherTuna","pinfish", "pompano", "redDrum","redGrouper",
               "redSnapper", "scamp","sciaenidae", "seatrout","shallowSerranidae",
               "sheepshead","skatesAndRays", "smallDemersalFish","smallPelagicFish",
               "smallReef fish","smallSharks","snook","spanishMackerel","spanishSardine","swordfish", 
               "vermilionSnapper", "yellowfinTuna")

#will use last 5 years of biomass
#target.time <- c(5474.5)
target.time <- c(5370.0,5400.0,5430.0,5460.0,5474.5)

train.files <- 1:1000

#change this for test or train sets
index.files <- train.files # test.files
pathfile <- "emulator_data_m.csv" # 


for(eachindex in train.files){
  
  #get biomass file
  this.biom <- paste("GOM_OUTBiomIndx_",eachindex,".txt",sep="")
  
  print(paste(this.biom, eachindex)) 
  
  #read biomass
  setwd(result.folders)
  run.biom <- try(read_delim(this.biom, delim=" ",na=c("NA","N/A",""," "),col_names=TRUE)) %>% 
    as.data.frame() %>% .[,1:92] %>% tbl_df %>% 
    filter(Time %in% target.time) %>% 
    dplyr::select(-Time) %>% 
    summarise_all(funs(mean)) %>% 
    mutate(run=eachindex) %>% 
    dplyr::select(run,everything()) %>% 
    gather(group,biomass,-run)
  
  #get catch file
  this.catch <- paste("GOM_OUTCatch_",eachindex,".txt",sep="")
  
  print(this.catch) 
  
  #read biomass
  setwd(result.folders)
  run.catch <- fread(this.catch, header=TRUE) %>% 
    dplyr::select(Time:SQU) %>% 
    filter(Time %in% target.time) %>% 
    dplyr::select(-Time) %>% 
    summarise_all(funs(mean)) %>% 
    mutate(run=eachindex) %>% 
    dplyr::select(run,everything()) %>% 
    gather(group,catch,-run)
  
  
  all.data <- run.catch %>% 
    full_join(run.biom) %>% 
    mutate(catch=(ifelse(is.na(catch),0,catch)))
  
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
  
  #replicate availability matrix
  prey.data.brt <- do.call("rbind", replicate(91, all.prey.data, simplify = FALSE))
  
  #combine biomass and availability matrix
 brt.frame <- all.data %>% 
   bind_cols(prey.data.brt) 
 
 setwd(work.folder)
 names(brt.frame) %>% 
   write.csv("availability_names.csv")
  
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

# will use polygons most affected by oil, same as Ainsworth et al. 2017
target.cell <- c(20, 21, 23, 31, 38, 39, 43, 52, 53, 54, 56, 63)
target.time <- c(5370.0,5400.0,5430.0,5460.0,5474.5)

setwd(work.folder)

species.names <- fread("Group_Names.csv", header=TRUE) %>% tbl_df %>% 
  mutate(`predator name`=tolower(`Long Name`)) %>% 
  rename(group=predator)

test.files <- setwd(result.folders) %>% 
  list.files(pattern="BoxBiomass_")

index.files <- 1:length(test.files)
pathfile <- "cell_biomass.csv" 

for(eachindex in index.files){
  
  this.biom <- test.files[eachindex]
  print(paste(this.biom, eachindex, ":",length(test.files)))
        
  #read biomass
  setwd(result.folders)
  run.biom <- try(read_delim(this.biom, delim=" ",na=c("NA","N/A",""," "),col_names=TRUE)) %>% 
    as.data.frame() %>% 
    tbl_df %>% 
    filter(Box %in% target.cell) %>% 
    filter(Time %in% target.time) %>% 
    dplyr::select(-Time, -Box) %>% 
    summarise_all(funs(mean)) %>% 
    mutate(run=eachindex) %>% 
    dplyr::select(run,everything()) %>% 
    gather(group,cell_biomass,-run) %>% 
    left_join(species.names,by="group")

  if(eachindex == 1){
    
    setwd(work.folder)
    write_delim(run.biom, path=pathfile, delim=",", append = FALSE, col_names = TRUE)
  } else{
    
    setwd(work.folder)
    write_csv(run.biom, path=pathfile,append = TRUE, col_names = FALSE)
    
  }
  
}


# Following commands used a Azure HDInsight cluster
#load file into Spark

# specifying yarn-client as our master, yarn is Hadoop NextGen
# http://spark.apache.org/docs/latest/running-on-yarn.html

work.folder <- "~/GOMcode"
setwd(work.folder)

config <- spark_config() 
config$spark.network.timeout <- "10d"

sc <- spark_connect (master = "yarn-client",  config = config)

#sc <- spark_connect (master = "yarn-client")

#check spark settings
print(sc)


#spark reads data set
#data should be stored in default container in Azure storage account.
#syntax is wasb[s]://<BlobStorageContainerName>@<StorageAccountName>.blob.core.windows.net/<path>
# unencrypted access (with the wasb: prefix) and SSL encrypted access (with wasbs). 
# https://docs.microsoft.com/en-us/azure/hdinsight/hdinsight-hadoop-use-blob-storage

origins <- file.path("wasbs://morzaria-2017-07-06t18-40-30-552z@morzariadatasc.blob.core.windows.net")

#data was partitioned above into train (75% simulations) and test (25%) sets

emulator_tbl <- spark_read_csv(sc,
                               path = origins,
                               name = 'emulator_data_y',
                               header = TRUE, overwrite = TRUE)

mean_emulator_tbl <- spark_read_csv(sc,
                                    path = origins,
                                    name = 'emulator_data_m',
                                    header = TRUE,
                                    overwrite = TRUE)

train.set <- 1:750
test.set <- 751:1000

this.data <- mean_emulator_tbl %>% 
  
  
  runs <- this.data %>% distinct(run)
test.run <- collect(runs)

train_tbl_b <- this.data %>% 
  filter(run %in% train.set)

test_tbl_b <- this.data %>% 
  filter(run %in% test.set)

train_tbl_c <- this.data %>% 
  filter(run %in% train.set) %>% 
  filter(catch!=0)

test_tbl_c <- this.data %>% 
  filter(run %in% test.set) %>% 
  filter(catch!=0)


# spark_write_parquet(train_tbl_b, path="wasbs://meantrainbiomass@morzariadatasc.blob.core.windows.net")
#                     
# spark_write_parquet(test_tbl_b, path=origins)
# 
# spark_write_parquet(train_tbl_c, path=origins)
# 
# spark_write_parquet(test_tbl_c, path=origins)
# 

#list spark tables
src_tbls(sc)

feature.cols <- this.data %>% 
  dplyr::select(gag_grouper:refractory_detritus_sediment) %>% 
  colnames()


# Bundle the models into a single list object
ml_models <- function(data, variable){
  
  ## Gradient boosted trees
  ml_gbt <- data %>% 
    ml_gradient_boosted_trees(variable, feature.cols, type="regression")  
  
  ## Decision Tree
  ml_dt <- data %>% 
    ml_decision_tree(variable, feature.cols, type="regression")
  
  ## Random Forest
  ml_rf <- data %>% 
    ml_random_forest(variable, feature.cols, type="regression")
  
  
  return(list(
    "Decision Tree" = ml_dt,
    "Random Forest" = ml_rf,
    "Gradient Boosted Trees" = ml_gbt
  ))
}


# run all the models
bio.emulator <- ml_models(train_tbl_b,"biomass")
catch.emulator <- ml_models(train_tbl_c,"catch")

# Create a function for scoring biomass
score_test_data <- function(model, data=test_tbl_b){
  pred <- sdf_predict(model, data)
  dplyr::select(pred, biomass, prediction) 
}

# Score all the models
ml_score <- lapply(bio.emulator, score_test_data)


# Lift function
calculate_lift <- function(ml_score) {
  ml_score %>%
    mutate(bin = ntile(desc(prediction), 10)) %>% 
    group_by(bin) %>% 
    summarize(count = sum(biomass)) %>% 
    mutate(prop = count / sum(count)) %>% 
    arrange(bin) %>% 
    mutate(prop = cumsum(prop)) %>% 
    dplyr::select(-count) %>% 
    collect() %>% 
    as.data.frame()
}

# Initialize results
ml_gains <- data.frame(bin = 1:10, prop = seq(0, 1, len = 10), model = "Base")

# Calculate lift
for(i in names(ml_score)){
  ml_gains <- ml_score[[i]] %>%
    calculate_lift %>%
    mutate(model = i) %>%
    rbind(ml_gains, .)
}

# Plot results
ggplot(ml_gains, aes(x = bin, y = prop, colour = model)) +
  geom_point() + geom_line() +
  ggtitle("Lift chart for predicting biomass - Test Data Set") + 
  xlab("") + ylab("")

# Calculate AUC and accuracy
perf_metrics <- data.frame(
  model = names(ml_score),
  AUC = 100 * sapply(ml_score, ml_classification_eval, "biomass", "prediction"),
  Accuracy = 100 *  sapply(ml_score, ml_classification_eval, "biomass", "prediction",metric="accuracy"),
  Precision = 100 *  sapply(ml_score, ml_classification_eval, "biomass", "prediction",metric="precision"),
  row.names = NULL, stringsAsFactors = FALSE)

# Plot results
gather(perf_metrics, metric, value, AUC, Accuracy, Precision) %>%
  ggplot(aes(reorder(model, value), value, fill = metric)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  coord_flip() +
  xlab("") +
  ylab("Percent") +
  ggtitle("Performance Metrics")

# Initialize results
feature_importance <- data.frame()

# Calculate feature importance
for(i in c("Decision Tree", "Random Forest", "Gradient Boosted Trees")){
  feature_importance <- ml_tree_feature_importance(sc, ml_models[[i]]) %>%
    mutate(Model = i) %>%
    mutate(importance = as.numeric(levels(importance))[importance]) %>%
    mutate(feature = as.character(feature)) %>%
    rbind(feature_importance, .)
}

# Plot results
feature_importance %>%
  ggplot(aes(reorder(feature, importance), importance, fill = Model)) + 
  facet_wrap(~Model) +
  geom_bar(stat = "identity") + 
  coord_flip() +
  xlab("") +
  ggtitle("Feature Importance")


# Run models for average of last 5 years of data
#data was partitioned above into train (75% simulations) and test (25%) sets

train_tbl_a <- spark_read_csv(sc,
                              path = origins,
                              name = 'emulator_data_train',
                              header = TRUE)

test_tbl_a <- spark_read_csv(sc,
                             path = origins,
                             name = 'emulator_data_test',
                             header = TRUE)



#list spark tables
src_tbls(sc)

#em_data %>%
#  filter(biomass==3.734023, gag_grouper==0, red_grouper==0, scamp == 0)

feature.cols_a <- train_tbl_a %>% 
  dplyr::select(gag_grouper:refractory_detritus_sediment) %>% 
  colnames()

## Generalized linear regression # NOT WORKING because some response variables were all 0s
ml_glr_a <- train_tbl_a %>% 
  ml_generalized_linear_regression("biomass", feature.cols_a)

## Decision Tree
ml_dt_a <- train_tbl_a %>% 
  ml_decision_tree("biomass", feature.cols, type="regression")

## Random Forest
ml_rf_a <- train_tbl_a %>% 
  ml_random_forest("biomass", feature.cols, type="regression")

## Neural Network NOT SURE WHY NOT WORKING
#ml_nn <- train_tbl %>%
#  ml_multilayer_perceptron("biomass", feature.cols, layers = c(94,50,20,10))

## Gradient boosted trees
ml_gbt_a <- train_tbl_a %>% 
  ml_gradient_boosted_trees("biomass", feature.cols, type="regression")  


# Bundle the modelss into a single list object
ml_models_a <- list(
  "Decision Tree" = ml_dt_a,
  "Random Forest" = ml_rf_a,
  "Gradient Boosted Trees" = ml_gbt_a
)

# Create a function for scoring
score_test_data_a <- function(model, data=test_tbl_a){
  pred <- sdf_predict(model, data)
  dplyr::select(pred, biomass, prediction)
}

# Score all the models
ml_score_a <- lapply(ml_models_a, score_test_data_a)


# Lift function
calculate_lift_a <- function(ml_score_a) {
  ml_score_a %>%
    mutate(bin = ntile(desc(prediction), 10)) %>% 
    group_by(bin) %>% 
    summarize(count = sum(biomass)) %>% 
    mutate(prop = count / sum(count)) %>% 
    arrange(bin) %>% 
    mutate(prop = cumsum(prop)) %>% 
    dplyr::select(-count) %>% 
    collect() %>% 
    as.data.frame()
}

# Initialize results
ml_gains_a <- data.frame(bin = 1:10, prop = seq(0, 1, len = 10), model = "Base")

# Calculate lift
for(i in names(ml_score_a)){
  ml_gains_a <- ml_score_a[[i]] %>%
    calculate_lift_a %>%
    mutate(model = i) %>%
    rbind(ml_gains, .)
}

# Plot results
ggplot(ml_gains_a, aes(x = bin, y = prop, colour = model)) +
  geom_point() + geom_line() +
  ggtitle("Lift chart for predicting biomass - Test Data Set") + 
  xlab("") + ylab("")

# Calculate AUC and accuracy
perf_metrics_a <- data.frame(
  model = names(ml_score_a),
  AUC = 100 * sapply(ml_score_a, ml_classification_eval, "biomass", "prediction"),
  Accuracy = 100 *  sapply(ml_score_a, ml_classification_eval, "biomass", "prediction",metric="accuracy"),
  Precision = 100 *  sapply(ml_score_a, ml_classification_eval, "biomass", "prediction",metric="precision"),
  row.names = NULL, stringsAsFactors = FALSE)

# Plot results
gather(perf_metrics_a, metric, value, AUC, Accuracy, Precision) %>%
  ggplot(aes(reorder(model, value), value, fill = metric)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  coord_flip() +
  xlab("") +
  ylab("Percent") +
  ggtitle("Performance Metrics")

# Initialize results
feature_importance_a <- data.frame()

# Calculate feature importance
for(i in c("Decision Tree", "Random Forest", "Gradient Boosted Trees")){
  feature_importance_a <- ml_tree_feature_importance(sc, ml_models_a[[i]]) %>%
    mutate(Model = i) %>%
    mutate(importance = as.numeric(levels(importance))[importance]) %>%
    mutate(feature = as.character(feature)) %>%
    rbind(feature_importance, .)
}

# Plot results
feature_importance_a %>%
  ggplot(aes(reorder(feature, importance), importance, fill = Model)) + 
  facet_wrap(~Model) +
  geom_bar(stat = "identity") + 
  coord_flip() +
  xlab("") +
  ggtitle("Feature Importance")

spark_disconnect(sc)

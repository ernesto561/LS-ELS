library(RSAGA)
library(exactextractr)
library(tidyverse)
library(sf)
library(terra)
library(pROC)
library(DescTools)
library(caret)
library(units)
library(tmap)
library(xgboost)
library(Matrix)
library(parallel)
library(doParallel)

Sys.setenv(TZ='America/El_Salvador')
Sys.setlocale("LC_ALL", "en_US.UTF-8")

n <- 1
cluster <- makeCluster(detectCores() - n) #  n is the number of remaining cores; convention to leave 1 core for OS
registerDoParallel(cluster)

#path to SAGA executable
#At this time, RSAGA only works with SAGA no greater than 8.4.1
env <- rsaga.env(r'(C:\Users\mreyes.AMBIENTE\saga-8.4.1_x64)')
#env <- rsaga.env(r'(C:\Users\ernes\saga-8.4.1_x64)')

##Read data

#folder name
area <- "els_alos"

#El Salvador outline (for maps)
els_lim <- read_sf("input/limits_els/els_outline.shp")

#Slope units
su <- read_sf("input/su/SU_ELS_limits.shp") 
#Landslides
ls_point <- read_sf("input/landslides/landslides_test.shp") 

#Checks if at least one landslide is inside a slope unit
su$frane <- lengths(st_intersects(su, ls_point)) > 0
su$frane <- as.integer(ifelse(su$frane == "TRUE", 1, 0))

#Plot of su and landslides
ggplot()+geom_sf(data = su, aes(fill=factor(frane)),linewidth=0.08)+
  scale_fill_discrete(name = "Independent variable", labels = c("no landslides", "landslides"))+
  theme_bw(14) 

ggsave("output/ls_su.png", width = 12, height = 6)

# #DEM
# dem <- rast(paste0("input/continuous/", area, "/els_alos_30m.tif"))
#for some reason, when the file is converted from WGS84 to UTM in QGIS the dimensions are slightly different
# 
# ####Variables preprocessing####
# terra::writeRaster(dem, paste0("input/continuous/", area, "/dem.sdat"), overwrite=TRUE)
#terra::writeRaster fails with large rasters. Converted to sdat with 30 m in both x and y resolution in QGIS.

#The input is a DEM in sdat format called "dem.sdat"

rsaga.slope.asp.curv(in.dem = paste0("input/continuous/", area, "/dem.sdat"),
                     out.slope = paste0("input/continuous/", area, "/slope"), unit.slope = 1,
                     out.cprof = paste0("input/continuous/", area, "/cprof"), out.cplan = paste0("input/continuous/", area, "/cplan"),
                     out.aspect = paste0("input/discrete/", area, "/asp.dat"), unit.aspect = 1,
                     method = "poly2zevenbergen",
                     env = env)

asp <- rast(paste0("input/discrete/", area, "/asp.sdat"))
m <- c(0, 22.5, 1,
       22.5, 67.5, 2,
       67.5, 112.5, 3,
       112.5, 157.5, 4,
       157.5, 202.5, 5,
       202.5, 247.5, 6,
       247.5, 292.5, 7,
       292.5, 337.5, 8,
       337.5, 360, 1)
rcl <- matrix(m, ncol=3, byrow=TRUE)
asp <- terra::classify(asp, rcl)
writeRaster(asp, paste0("input/discrete/", area, "/asp.sdat"), overwrite=TRUE)
 
################################
## Slope unit analysis
################################

crs_els <- st_crs(su)

read_r <- function(x){r <- rast(x);crs(r)<-crs_els$wkt;return(r)}

#Continuous variables
vars_cont <- map(list.files(paste0("input/continuous/", area, "/"), pattern="*.sdat$", full.names = T), read_r)
su_vars_cont <- exact_extract(rast(vars_cont), su, c('median', 'stdev'), force_df =T, full_colnames = T, max_cells_in_memory = 12e+07)

#Discrete variables
vars_disc <- map(list.files(paste0("input/discrete/", area, "/"), pattern="*.sdat$", full.names = T), read_r)
su_vars_disc <- exact_extract(rast(vars_disc), su, 'majority', force_df =T, full_colnames = T, max_cells_in_memory = 12e+07)

su_model <- data.frame(cbind(su, su_vars_cont, su_vars_disc))

all<-data.frame(su_model)


#Removes NoData Values
all=na.exclude(all)
rownames(all) <- NULL

all$frane <- as.factor(all$frane)
all$majority.asp <- as.factor(all$majority.asp)

#Random selection of negative slope units
all_random <- all
#Selection of negative slope units using median slope < 5°
all_slo5 <- all %>% dplyr::filter(frane == 1 | (frane == 0 & median.slope <=5))

##Selection of negative slope units using Mahalanobis Distance°
#Categorical variables are excluded. 

md_vars_ls <- all %>% dplyr::filter(frane==1) %>%
  dplyr::select(-c(DN, geometry, frane, majority.asp)) 

all_md <- all %>% 
  mutate(md = mahalanobis(all %>% 
                            dplyr::select(-c(DN, geometry, frane, majority.asp)), colMeans(md_vars_ls), cov(md_vars_ls)))

cutoff <- qchisq(p = 0.95 , df = ncol(md_vars_ls)-1)

all_md_cut <- all_md %>% dplyr::filter(frane == 1 | (frane == 0 & md>cutoff)) %>% 
  dplyr::select(-c(md))

# #Mahalanobis distance plot
# ggplot()+geom_sf(data = st_as_sf(all_md), aes(fill=md),color=NA)+scale_fill_viridis_c()+theme_bw(14)+
#   labs(title = "Malahalobis distance for slope units")+theme(plot.title.position = "plot")
# 
# ggsave("output/md.png", width = 12, height = 6)

ggplot(st_as_sf(all_md)) +
  geom_sf(aes(fill=md), data = ~ subset(., md < 1000), color=NA) +
  geom_sf(color = "orange", data = ~ subset(., md >= 1000), size=20) +
  scale_fill_gradient(low="blue", high="red")+theme_bw(14)+
  labs(title = "Mahalanobis distance for slope units")+theme(plot.title.position = "plot")

ggsave("output/md.png", width = 12, height = 6)




##########################
######Variable maps######
##########################

map_vars <- function(sf, var){
  name <- "var"
  titulo <- name
  tm <- tm_shape(sf)+
    tm_fill(var)+
    tm_shape(els_lim)+
    tm_borders()
}

su_model_sf <- st_as_sf(su_model)

vars_maps <- names(dplyr::select(su_model, -c(DN, frane, geometry)))
map_list <- pmap(list(list(su_model_sf), vars_maps), map_vars)

varmap1 <- tmap_arrange(map_list[1:4])
tmap_save(varmap1, filename = "output/varmap1.png")
varmap2 <- tmap_arrange(map_list[5:9])
tmap_save(varmap2, filename = "output/varmap2.png")


################################
#Variable plots
################################

vars_to_plot <- dplyr::select(all, -c(DN, geometry))

#ggpairs(vars_to_plot)

################################
#Calibration and validation data
################################

#For this function, df is a dataframe with all the variables and per is
#the percentage of calibration data

#Returns a list with a df with the variables with the training data, a list with the Y training values labels,
#a df with the samples used for training the model (for ploting the samples), a df with the test data, 
#a list with the test labels and a df with all the values for prediction

cal_val <- function(df, per){

  X_data <- df
  y_data <- df[["frane"]]
  
  #Split the data set
  set.seed(50)
   
  # 70/30 split
  split_indices <- createDataPartition(y_data, p = per, list = FALSE)
  X_tr <- X_data[split_indices, ]
  #Downsampling in order to create a balanced dataset
  X_tr_bal <- downSample(X_tr%>%dplyr::select(-c(frane)), X_tr$frane, yname = "frane")
  #Dummy contrast coding for categorical variables 
  #(https://cran.r-project.org/web/packages/xgboost/vignettes/discoverYourData.html)
  X_train <- sparse.model.matrix(frane ~ ., data = X_tr_bal%>%dplyr::select(-c(DN, geometry)))[,-1]
  Y_train <- as.numeric(as.character(X_tr_bal$frane))
  X_tst <- X_data[-split_indices, ] %>% relocate(frane, .after = last_col())
  X_test <- sparse.model.matrix(frane ~ ., data = X_tst%>%dplyr::select(-c(DN, geometry)))[,-1] 
  Y_test <- as.numeric(as.character(X_tst$frane))
  
  X_data_pred <- sparse.model.matrix(frane ~ ., data = X_data %>%dplyr::select(-c(DN, geometry)))[,-1]
  
  return(list(X_train, Y_train, X_tr_bal, X_test, Y_test, X_data_pred))

}

calval_random <- cal_val(all_random, 0.7)
calval_slo5 <- cal_val(all_slo5, 0.7)
calval_md <- cal_val(all_md_cut, 0.7)

#######################
##Plots samples########
#######################

###Maps###
su_map <- function(sf, var, n, style, palette, title){
  name <- "var"
  titulo <- name
  tm <- tm_shape(sf)+
    tm_fill(var, n=n, style = style, palette = palette )+
    tm_shape(els_lim)+
    tm_borders(lwd = 2)+
    tm_layout(
      title = title,
      title.position = c("left", "top"),
      title.size = 1.1,
      bg.color = "#fcfcfc",
      inner.margins = c(0.06, 0.01, 0.09, 0.01),
      outer.margins = 0,
      frame.lwd = 0.2
    )
  #return(tm)
  tmap_save(tm=tm, filename=paste0("output/", title, ".png"))
}

su_map(st_as_sf(calval_random[[3]]), "frane", 2,"cat", "-RdYlGn", "XGBoost - random negative samples" )
su_map(st_as_sf(calval_slo5[[3]]), "frane", 2,"cat", "-RdYlGn", "XGBoost - slope less than 5" )
su_map(st_as_sf(calval_md[[3]]), "frane", 2,"cat", "-RdYlGn", "Mahalanobis distance" )


#####################
####Modeling#########
#####################


#############################
####Hyperparameter tuning####
#############################

hyperparam_grid <- expand.grid(
  nrounds = seq(from = 100, to = 300, by = 100),
  eta = c(0.025, 0.05, 0.1, 0.3),
  max_depth = c(4, 5, 6),
  gamma = c(0, 1, 2),
  colsample_bytree = c(0.5, 0.75, 1.0),
  min_child_weight = c(1, 3, 5),
  subsample = 1
)

tune_control <- caret::trainControl(
  method = "cv", # cross-validation
  number = 10, # with n folds
  verboseIter = FALSE, # no training log
  allowParallel = TRUE,
  savePredictions = T,
  classProbs = T,
  summaryFunction = twoClassSummary
)

tuned_model <- function(X_train, Y_train){
  y_tr = as.factor(Y_train)
  levels(y_tr)=c("No","Yes")
  bst <- caret::train(
    x = X_train,
    y = y_tr,
    trControl = tune_control,
    tuneGrid = hyperparam_grid,
    method = "xgbTree", #  to use XGB
    verbose = FALSE,
    verbosity = 0
)
  return(bst)

}

model_cal_rdm <- tuned_model(calval_random[[1]], calval_random[[2]])
model_cal_slo5 <- tuned_model(calval_slo5[[1]], calval_slo5[[2]])
model_cal_md <- tuned_model(calval_md[[1]], calval_md[[2]])


stopCluster(cluster)

registerDoSEQ()

unregister <- function() {
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
}

# run xgboost with optimal values
# The input is the the output of calval function and the output of tuned model
model_opt <- function(calval, model_cal){
xgboost(
  booster= "gbtree",#model used
  objective = "binary:logistic",
  data = calval[[1]],#data
  label = calval[[2]],#binary outcome
  eta = model_cal$bestTune$eta,#learning rate, determines shrinkage of the iterations
  verbose = FALSE,
  nrounds = model_cal$bestTune$nrounds,
  #optimized parameters
  max_depth = model_cal$bestTune$max_depth,#depth of tree higher values may lead to overfitting
  min_child_weight = model_cal$bestTune$min_child_weight,
  subsample = model_cal$bestTune$subsample,
  gamma = model_cal$bestTune$gamma, #Minimum loss reduction required to make a further partition on a leaf node of the tree. The larger gamma is, the more conservative the algorithm will be.
  colsample_bytree = model_cal$bestTune$colsample_bytree,
  eval_metric = 'auc',
)
}

model_opt_random <- model_opt(calval_random, model_cal_rdm)
model_opt_slo5 <- model_opt(calval_slo5, model_cal_slo5)
model_opt_md <- model_opt(calval_md, model_cal_md)


#Predictions with test data
model_random_pred_tst = predict(model_opt_random, calval_random[[4]], type = "response")
model_slo5_pred_tst = predict(model_opt_slo5, calval_slo5[[4]], type = "response")
model_md_pred_tst = predict(model_opt_md, calval_md[[4]], type = "response")


#AUC test data
auc_test_random <- roc(response = calval_random[[5]], predictor = model_random_pred_tst)
auc_test_slo5 <- roc(response = calval_slo5[[5]], predictor = model_slo5_pred_tst)
auc_test_md <- roc(response = calval_md[[5]], predictor = model_md_pred_tst)


auc_test_random_value <-auc(auc_test_random)
auc_test_slo5_value <-auc(auc_test_slo5)
auc_test_md_value <-auc(auc_test_md)


auc_test_values <- data.frame(model = c("random", "slope < 5", "md"), AUC = c(auc_test_random_value, auc_test_slo5_value, auc_test_md_value))


#AUC for every fold with training data
#https://stackoverflow.com/a/69261452/4268720

auc_folds <- function(model_cal){
  sapply(X = unique(model_cal$pred$Resample),
                 FUN = function(x) {
                   r <- model_cal$pred[model_cal$pred$Resample == x,]
                   R <- auc(response = r$obs, predictor = r$Yes)
                   return(R)
                 }, simplify = T) %>%
  enframe() 
}

auc_folds_rdm = auc_folds(model_cal_rdm) %>% mutate(model = "random")
auc_folds_slo5 = auc_folds(model_cal_slo5) %>% mutate(model = "slope < 5")
auc_folds_md = auc_folds(model_cal_md) %>% mutate(model = "md")


auc_folds_df <- bind_rows(auc_folds_rdm, auc_folds_slo5, auc_folds_md)

p <- ggplot(auc_folds_df, aes(x=model, y=value)) + 
  geom_boxplot()+geom_point(data = auc_test_values, aes(x=model, y=AUC))+labs(y="AUC")+theme_bw()


#Brier score for every fold with training data
bs_folds <- function(model_cal){
  sapply(X = unique(model_cal$pred$Resample),
         FUN = function(x) {
           r <- model_cal$pred[model_cal$pred$Resample == x,]
           R <- BrierScore(x = r$obs, pred = r$Yes)
           return(R)
         }, simplify = T) %>%
    enframe() 
}

bs_folds_rdm = bs_folds(model_cal_rdm) %>% mutate(model = "random")
bs_folds_slo5 = bs_folds(model_cal_slo5) %>% mutate(model = "slope < 5")
bs_folds_md = bs_folds(model_cal_md) %>% mutate(model = "md")


bs_folds_df <- bind_rows(bs_folds_rdm, bs_folds_slo5, bs_folds_md)

bs_p <- ggplot(bs_folds_df, aes(x=model, y=value)) + 
  geom_boxplot()+geom_point(data = auc_test_values, aes(x=model, y=AUC))+labs(y="AUC")+theme_bw()

###################################
####Predictions with all data######
###################################

#For the prediction data are use the calval_random data, because is the total of data
model_random_resp_all = predict(model_opt_random, calval_random[[6]], type = "response")
model_slo5_resp_all = predict(model_opt_slo5, calval_random[[6]], type = "response")
model_md_resp_all = predict(model_opt_md, calval_random[[6]], type = "response")


su_random <- bind_cols(all, data.frame(prob=model_random_resp_all))
su_slo5 <- bind_cols(all, data.frame(prob=model_slo5_resp_all))
su_md <- bind_cols(all, data.frame(prob=model_md_resp_all))


###################
####Maps########
###################

su_map(st_as_sf(su_random), "prob", 4, "jenks", "-RdYlGn", "XGBoost - random negative samples")
su_map(st_as_sf(su_slo5), "prob", 4, "jenks", "-RdYlGn", "XGBoost - samples slo<5")
su_map(st_as_sf(su_md), "prob", 4, "jenks", "-RdYlGn", "XGBoost - Mahalanobis distance")





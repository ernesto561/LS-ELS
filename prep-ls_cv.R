library(RSAGA)
library(exactextractr)
library(tidyverse)
library(sf)
library(terra)
library(pROC)
library(caret)
library(ROCR)
library(units)
library(tmap)
library(xgboost)
library(Matrix)
library(parallel)
library(doParallel)

<<<<<<< HEAD
n <- 4
cluster <- makeCluster(detectCores() - n) #  n is the number of cores remaining; convention to leave 1 core for OS
=======
cluster <- makeCluster(detectCores() - 1) # convention to leave 1 core for OS
>>>>>>> 8461949c4ea3d3b165cbe8352d564c3a0aea6fb9
registerDoParallel(cluster)

#path to SAGA executable
#At this time, RSAGA only works with SAGA 8.4.1
#env <- rsaga.env(r'(C:\Users\mreyes.AMBIENTE\saga-8.4.1_x64)')
env <- rsaga.env(r'(C:\Users\ernes\saga-8.4.1_x64)')

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
ggplot()+geom_sf(data = su, aes(fill=factor(frane)),linewidth=0.02)+
  scale_fill_discrete(labels=c('No landslides', 'landslides'), name="Landslide inventory")+
  theme_bw(16) 

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

all_random <- all
all_slo5 <- all %>% dplyr::filter(frane == 1 | (frane == 0 & median.slope <=5))

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
#Calibration and validation data
################################

#For this function, df is a dataframe with all the variables and per is
#the percentage of calibration data

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
  #Dummy contrast coding for categorical variables (https://cran.r-project.org/web/packages/xgboost/vignettes/discoverYourData.html)
  X_train <- sparse.model.matrix(frane ~ ., data = X_tr_bal%>%dplyr::select(-c(DN, geometry)))[,-1]
  Y_train <- as.numeric(as.character(X_tr_bal$frane))
  X_tst <- X_data[-split_indices, ] %>% relocate(frane, .after = last_col())
  X_test <- sparse.model.matrix(frane ~ ., data = X_tst%>%dplyr::select(-c(DN, geometry)))[,-1] 
  
  X_data_pred <- sparse.model.matrix(frane ~ ., data = X_data %>%dplyr::select(-c(DN, geometry)))[,-1]
  
  return(list(X_train, Y_train, X_tr_bal, X_test, X_data_pred))

}

calval_random <- cal_val(all_random, 0.7)
calval_slo5 <- cal_val(all_slo5, 0.7)

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
  return(tm)
}

su_map(st_as_sf(calval_random[[3]]), "frane", 2,"cat", "-RdYlGn", "XGBoost - random negative samples" )
su_map(st_as_sf(calval_slo5[[3]]), "frane", 2,"cat", "-RdYlGn", "XGBoost - slope < 5" )


#####################
####Modeling#########
#####################
model <- function(X_train, Y_train){
   default_model <- xgboost(data = X_train,
                           label = Y_train,
                           booster = "gbtree",
                           objective = "binary:logistic",
                           nrounds = 100,
                           verbose = 0) 
return(default_model)

}

model_random <- model(calval_random[[1]], calval_random[[2]])
model_slo5 <- model(calval_slo5[[1]], calval_slo5[[2]])


##Importance matrix and plot

importance <- xgb.importance(feature_names = colnames(calval_random[[1]]), model = model_random)
xgb.ggplot.importance(importance_matrix = importance)+theme_bw()+theme(plot.title.position = "plot")

p_random <- predict(model_random, calval_random[[5]], type = "response")
pred_random <- bind_cols(all, data.frame(p_random))
su_map_random <- left_join(su_model, pred_random, by = c("DN"), keep = FALSE)

p_slo5 <- predict(model_slo5, calval_random[[5]], type = "response")
pred_slo5 <- bind_cols(all, data.frame(p_slo5))
su_map_slo5 <- left_join(su_model, pred_slo5, by = c("DN"), keep = FALSE)

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

tuned_model <- function(X_train, Y_train, hyperparam_grid){
  y_val = as.factor(Y_train)
  levels(y_val)=c("No","Yes")
  bst <- caret::train(
    x = X_train,
    y = y_val,
    trControl = tune_control,
    tuneGrid = hyperparam_grid,
    method = "xgbTree", #  to use XGB
    verbose = FALSE,
    verbosity = 0
)
  return(bst)

}


model_cal <- tuned_model(calval_random[[1]], calval_random[[2]], hyperparam_grid)
model_random_final <- xgb.cv(data = calval_random[[1]],
                             label = calval_random[[2]],
                             booster = "gbtree",
                             objective = "binary:logistic",
                             nfold = 10,
                             prediction = TRUE,
                             eval_metric = "auc",
                             nrounds=model_cal$bestTune$nrounds,
                             eta=model_cal$bestTune$eta,
                             max_depth=model_cal$bestTune$max_depth,
                             gamma=model_cal$bestTune$gamma,
                             colsample_bytree=model_cal$bestTune$colsample_bytree,
                             min_child_weight=model_cal$bestTune$min_child_weight,
                             subsample=model_cal$bestTune$subsample) 

z <- lapply(model_random_final$folds, function(x){
  pred <- model_random_final$pred[x]
  true <- (calval_random[[2]])[x]
  index <- x
  out <- data.frame(pred, true, index)
  auc = data.frame(auc = print(auc(out$true, out$pred)))
  return(auc)
}) %>% bind_rows()

names(z) <- paste("fold", 1:10, sep = "_")

stopCluster(cluster)

pred <- predict(model_random_final, calval_random[[4]])


registerDoSEQ()

unregister <- function() {
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
}


###################
####all_map########
###################

su_map(st_as_sf(su_map_random), "p_random", 4, "jenks", "-RdYlGn", "XGBoost - random negative samples")
su_map(st_as_sf(su_map_slo5), "p_slo5", 4, "jenks", "-RdYlGn", "XGBoost - samples slo<5")


dd.roc <- sapply(X = unique(model_random_final$pred$Resample),
                 FUN = function(x) {
                   r <- model_random_final$pred[model_random_final$pred$Resample == x,]
                   return(r)
                   R <- auc(response = r$obs, predictor = r$Yes)
                   #data.frame(auc=R[[auc]])
                 }, simplify = F) %>%
  bind_rows(.id = "Resample") %>%
  as_tibble()


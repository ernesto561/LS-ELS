library(RSAGA)
library(exactextractr)
library(tidyverse)
library(sf)
library(terra)
library(pROC)
library(usdm) 
library(splitstackshape) 
library(caret)
library(ROCR)
library(sperrorest)
library(earth)
library(units)
library(dismo)
library(MIAmaxent)
library(ENMeval)
library(tmap)
library(caret)
library(GGally)
library(xgboost)

memory.limit(800000)

#path to SAGA executable
#At this time, RSAGA only works with SAGA 8.4.1
env <- rsaga.env(r'(C:\Users\mreyes.AMBIENTE\saga-8.4.1_x64)')
#env <- rsaga.env(r'(C:\Users\ernes\saga-8.4.1_x64)')

#number of replications
n <- 10

##Read data

#folder name
area <- "els_alos"

#El Salvador outline (for maps)
els_lim <- read_sf("input/limits_els/els_outline.shp")

#Slope units
su <- read_sf("input/su/SU_ELS.shp") 
#Landslides
ls_point <- read_sf("input/landslides/landslides_test.shp") 

#Checks if at least one landslide is inside a slope unit
su$frane <- lengths(st_intersects(su, ls_point)) > 0
su$frane <- as.integer(ifelse(su$frane == "TRUE", 1, 0))

#Plot of su and landslides
ggplot()+geom_sf(data = su, aes(fill=factor(frane)),linewidth=0.02)+theme_bw() 

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

# rsaga.geoprocessor(lib = "ta_morphometry",
#                    module = "TPI Based Landform Classification",
#                    param = list(DEM = paste0("input/continuous/", area, "/dem.sdat"),
#                                 LANDFORMS = paste0("input/discrete/", area, "/lcl.sdat"),
#                                 RADIUS_A_MIN = 0,
#                                 RADIUS_A_MAX = 100,
#                                 RADIUS_B_MIN = 0,
#                                 RADIUS_B_MAX = 1000),
#                                 env = env)

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
su_vars_cont <- exact_extract(rast(vars_cont), su, c('median', 'stdev'), force_df =T, full_colnames = T)

#Discrete variables
vars_disc <- map(list.files(paste0("input/discrete/", area, "/"), pattern="*.sdat$", full.names = T), read_r)
su_vars_disc <- exact_extract(rast(vars_disc), su, 'majority', force_df =T, full_colnames = T)

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

# map_vars <- function(sf, var){
#   name <- "var"
#   titulo <- name
#   tm <- tm_shape(sf)+
#     tm_fill(var)+
#     tm_shape(els_lim)+
#     tm_borders()
# }
# 
# su_model_sf <- st_as_sf(su_model)
# 
# vars_maps <- names(dplyr::select(su_model, -c(DN, frane, geometry)))
# map_list <- pmap(list(list(su_model_sf), vars_maps), map_vars)
# 
# varmap1 <- tmap_arrange(map_list[1:4])
# tmap_save(varmap1, filename = "output/varmap1.png")
# varmap2 <- tmap_arrange(map_list[5:9])
# tmap_save(varmap2, filename = "output/varmap2.png")
# 

################################
#Variable plots
################################

vars_to_plot <- dplyr::select(all, -c(DN, geometry))

#ggpairs(vars_to_plot)

################################
#Calibration and validation data
################################

cal_val <- function(df){

  X_data <- df
  y_data <- df[["frane"]]
  
  X_data_pred <- X_data %>% dplyr::select(-c(DN, geometry))
  X_data_pred_m <- sparse.model.matrix(frane ~ ., data = X_data_pred)[,-1]
  
  #Split the data set
  set.seed(50)
  
  # 70/30 split
  split_indices <- createDataPartition(y_data, p = 0.7, list = FALSE)
  X_tr <- X_data[split_indices, ]
  #Downsampling in order to create a balanced dataset
  X_tr_bal <- downSample(X_tr%>%dplyr::select(-c(frane, DN, geometry)), X_tr$frane, yname = "frane") 
  #Dummy contrast coding for categorical variables (https://cran.r-project.org/web/packages/xgboost/vignettes/discoverYourData.html)
  X_train <- sparse_matrix <- sparse.model.matrix(frane ~ ., data = X_tr_bal)[,-1]
  Y_train <- as.numeric(as.character(X_tr_bal$frane)) 
  X_tst <- X_data[-split_indices, ]
  X_test <- X_tst %>% relocate(frane, .after = last_col())
  
  return(list(X_train, Y_train, X_test, X_data_pred_m))

}

calval_random <- cal_val(all_random)
calval_slo5 <- cal_val(all_slo5)

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

##Importance matrix and plot

importance <- xgb.importance(feature_names = colnames(calval_random[[1]]), model = model_random)

xgb.ggplot.importance(importance_matrix = importance)+theme_bw()+theme(plot.title.position = "plot")


p <- predict(model_random, calval_random[[4]], type = "response")

pred <- bind_cols(all, data.frame(p))

su_map_random <- left_join(su_model, pred, by = c("DN"),keep = FALSE)



###################
####all_map########
###################

#all is the input data, model is the output of

all_map <- function(all, model, name){

all_map<-matrix(nrow=nrow(all), ncol=n)

for(i  in 1:n){
  all_map[,i]<-predict(model[[i]], all, type=c("response"))
}

all_map_avg <-apply(all_map,1, mean)
all_map_avg<- round(all_map_avg, digits=4)

all_map_data<-data.frame(DN=all$DN,score=all_map_avg)

su_map = left_join(su, all_map_data)
sf::st_write(su_map, paste0("output/", name, ".shp"), append = FALSE)

return(su_map)

}

els_random <- all_map(all_random, mars_random[[1]], "els_random")
els_slo5 <- all_map(all_slo5, mars_slo5[[1]], "els_slo5")



###Final map###
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

su_map(els_random, "score", 4, "jenks", "-RdYlGn", "Landslide susceptibility - random negative samples")
su_map(els_slo5, "score", 4, "jenks", "-RdYlGn", "Landslide susceptibility - slopes <=5")

su_map(st_as_sf(su_map_random), "p", 4, "jenks", "-RdYlGn", "XGBoost")


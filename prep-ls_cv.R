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
#env <- rsaga.env(r'(C:\Users\mreyes.AMBIENTE\saga-8.4.1_x64)')
env <- rsaga.env(r'(C:\Users\ernes\saga-8.4.1_x64)')

#number of replications
n <- 10

##Read data

#folder name
area <- "els_alos"

#El Salvador outline (for maps)
#els_lim <- read_sf("input/limits_els/els_outline.shp")

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
su_vars_cont <- exact_extract(rast(vars_cont), su, c('median', 'stdev'))

#Discrete variables
vars_disc <- map(list.files(paste0("input/discrete/", area, "/"), pattern="*.sdat$", full.names = T), read_r)
su_vars_disc <- exact_extract(rast(vars_disc), su, 'majority')

su_model <- data.frame(cbind(su, su_vars_cont, su_vars_disc))

all<-data.frame(su_model)


#Removes NoData Values
all=na.exclude(all)
rownames(all) <- NULL

all$frane <- as.factor(all$frane)

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

cal_val <- function(all){

  X_data <- all
  y_data <- all[["frane"]]
  
  #Split the data set
  set.seed(50)
  
  # 70/30 split
  split_indices <- createDataPartition(y_data, p = 0.7, list = FALSE)
  X_tr <- X_data[split_indices, ]
  X_train <- downSample(X_tr%>%dplyr::select(-c(frane)), X_tr$frane, yname = "frane")
  X_tst <- X_data[-split_indices, ]
  X_test <- X_tst %>% relocate(frane, .after = last_col())
  
  return(list(X_train, X_test))

}

calval_random <- cal_val(all_random)
calval_slo5 <- cal_val(all_slo5)

#####################
####Modeling#########
#####################
model <- function(X_train){
  X_train_sel <- X_train %>% dplyr::select(-c(DN, geometry))
  default_model <- xgboost(data = model.matrix(~.+0,data = X_train_sel[,-c("frane")]),
                           label = X_train_sel$frane,
                           booster = "gbtree",
                           objective = "binary:logistic",
                           nrounds = 100,
                           verbose = 0) 
return(default_model)

}

model_random <- model(calval_random[[1]])

xgbpred <- function(model, data, ...) {
  predict(model, newdata=as.matrix(data), ...)
}

p <- predict(default_model, as.matrix(X_data), type = "response")

pred <- bind_cols(all, data.frame(p))

############################
####ROC models function#####
############################

roc_model <- function(vec, all_predict){
  
  #ROC curve 
  casi_all<-NULL;score_all<-NULL; xrocall<-NULL;roc_all<-NULL
  
  for(i  in 1:n){
    casi_all[[i]]<- data.frame(vec[[i]]$frane)
  }
  
  for(i  in 1:n){
    casi_all[[i]]<-as.vector(casi_all[i])
    score_all[[i]]<-as.vector(all_predict[,i])
  }

  for(i in 1:n){
    xrocall[[i]]<-data.frame(score= score_all[[i]],response=casi_all[[i]])
  }
  
  for(i in 1:n){
    colnames(xrocall[[i]])<- c("score", "response")
  }
  
  for(i in 1:n){
    roc_all[[i]]<-roc(xrocall[[i]]$response ~ xrocall[[i]]$score, xrocall)}
  
  return(list(roc_all))

}

#Validation data
roc_random <- roc_model(calval_random[[2]], mars_random[[2]])
roc_slo5 <- roc_model(calval_slo5[[2]], mars_slo5[[2]])


#################################
####Mean ROC models function#####
#################################

#Input is vec (validation or calibration data) and all_predict (output of mars_model)

roc_model_m <- function(vec, all_predict){
casiallc<-matrix(nrow = nrow(vec[[1]]), ncol=n)

for(i  in 1:n){
  casiallc[,i]<-as.vector(vec[[i]]$frane)
}

roc_all_m<-NULL
casesall<-as.vector(casiallc)
scoreall<-as.vector(all_predict)
dataall<-data.frame(score=scoreall,response=casesall)
roc_all_m<-roc(dataall$response ~ dataall$score, dataall)
youdenall_m<- min(coords(roc_all_m, "b", ret="t", best.method="youden"))
youdenall_m<- round(youdenall_m, digits=3)
auc_all_m<-round(roc_all_m$auc, digits = 3)

return(list(roc_all_m, youdenall_m, auc_all_m))

}

roc_m_random <- roc_model_m(calval_random[[2]], mars_random[[2]])
roc_m_slo5 <- roc_model_m(calval_slo5[[2]], mars_slo5[[2]])

#################################
####Plot ROC models function#####
#################################

#Inputs are roc_all (output of roc_model) and roc_all_m (output of roc_model_m)  

plot_roc <- function(roc_all, roc_all_m){

mean_roc <- data.frame(x=1-roc_all_m[[1]]$specificities, y=roc_all_m[[1]]$sensitivities)

p<- ggroc(roc_all[[1]], legacy.axes = T)+geom_line(color="gray")+
  geom_line(data = mean_roc, aes(x, y), color="red", inherit.aes = FALSE)+
  annotate("text", x=0.75, y=0.75, label= paste0("AUC = ", roc_m_random[[3]]))+
  theme_bw()+theme(legend.position="none")+coord_equal()

return(p)

}

plot_roc(roc_random, roc_m_random)
plot_roc(roc_slo5, roc_m_slo5)

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

su_map(pp, "p", 4, "jenks", "-RdYlGn", "XGBoost")


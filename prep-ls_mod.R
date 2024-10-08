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

memory.limit(800000)

#path to SAGA executable
#At this time, RSAGA only works with SAGA 8.4.1
env <- rsaga.env(r'(C:\Users\mreyes.AMBIENTE\saga-8.4.1_x64)')
#env <- rsaga.env(r'(C:\Users\ernes\saga-8.4.1_x64)')

#number of replications
n <- 5

##Read data

#folder name
area <- "els_alos"

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
#terra::writeRaster fails with large rasters. Converted to sdat with 30 m resolution in QGIS.

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

##########################
######Variable plots######
##########################

map <- function(var){
  name <- "var"
  titulo <- name
  tm <- tm_shape(dem)+
    tm_raster(style="cont")+
    tm_layout(main.title = titulo, 
              main.title.size = 1.5, 
              main.title.position = c("left", "top"), 
              legend.title.size = 1.2,
              legend.text.size = 1,
              inner.margins = c(0.2, 0.2, 0.2, 0.2))+
    tm_scale_bar(position = c("left", "bottom"), text.size = 0.7)+
    tm_compass(position = c("right", "top"))#+
    # #tm_add_legend(type = "fill", 
    #               labels = c("Depósito"),
    #               col = NA,
    #               border.lwd = 1.5,
    #               border.col = "red")
  #tmap_save(tm, paste0("../imagenes/", zona_comp, name,".png"), asp=get_asp_ratio(tm))
  return(tm)
}



################################
#Calibration and validation data
################################

cal_val <- function(all){

all$frane<-as.factor(all$frane)

pos_all<-all[which(all$frane=="1"),]

####sampling all area####
all_cases<-NULL
all_sampling14<-NULL
all_cal14<-NULL
all_val14<-NULL

##usa questo ciclo for per prendere i casi in modo casuale(solo i negativi ovviamente) e tutti i positivi
## Traducci?n: se usa este bucle para tomar los casos al azar (solo los negativos, por supuesto) y todos los positivos

for(i  in 1:n){
  all_cases[[i]]<-stratified(all, "frane", nrow(pos_all))
  all_cases[[i]]$frane <- as.factor(all_cases[[i]]$frane)}

##Aqui se muestrean solo un porcentaje del total para los deslizamientos (p. ej. 0.7=70%, 3er parametro de la funcion stratified de splitstackshape)
for(i  in 1:n){
  all_sampling14[[i]]<-stratified(all_cases[[i]], "frane", .7, select = list(frane=c("0", "1")), bothSets = TRUE)
}

for(i  in 1:n){
  all_cal14[[i]]<-all_sampling14[[i]]$SAMP1
  all_cal14[[i]]<-data.frame(all_cal14[[i]])
}

for(i  in 1:n){
  all_val14[[i]]<-all_sampling14[[i]]$SAMP2
  all_val14[[i]]<-data.frame(all_val14[[i]])
}

return(list(all_cal14, all_val14))

}

calval_random <- cal_val(all)

#####################
####MARS modeling####
#####################
mars_model <- function(cal, val){

mars_all<-NULL
all_predict<-matrix(nrow=nrow(val[[1]]), ncol = n)
AUC_all<-matrix(nrow=1, ncol=n)

####Regression analysis
####Cambiar esta linea si se cambian variables####

for(i  in 1:n){
  mars_all[[i]]<- earth (reformulate(names(dplyr::select(all, -c(DN, frane, geometry))), "frane"), data = cal[[i]], trace = 1, degree=1,  glm=list(family=binomial), Scale.y=FALSE)}

for(i  in 1:n){
  all_predict[,i]<-predict(mars_all[[i]], val[[i]], type=c("response"))
}

summary_earth <- lapply(1:n, function(x) capture.output(summary(mars_all[[x]]), file = "salida/summary_earth.csv", append = TRUE))

return(all_predict)

}

mars_random <- mars_model(calval_random[[1]], calval_random[[2]])

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
roc_random <- roc_model(calval_random[[2]], mars_random)

#################################
####Mean ROC models function#####
#################################

#Input is vec(validation or calibration data) and all_predict (output of mars_model)

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

roc_m_random <- roc_model(calval_random[[2]], mars_random)

#################################
####Plot ROC models function#####
#################################

#Inputs are roc_all (output of roc_model) and roc_all_m (output of roc_model_m)  

plot_roc <- function(roc_all, roc_all_m){

mean_roc <- data.frame(x=1-roc_all_m$specificities, y=roc_all_m$sensitivities)

p<- ggroc(roc_all, legacy.axes = T)+geom_line(color="gray")+
  geom_line(data = mean_roc, aes(x, y), color="red", inherit.aes = FALSE)+theme_bw()+theme(legend.position="none")

return(p)

}

plot_roc(roc_random, roc_m_random)

###################
####all_map####
###################
all_map<-matrix(nrow=nrow(all), ncol=n)

for(i  in 1:n){
  all_map[,i]<-predict(mars_all[[i]], all, type=c("response"))
}

all_map_avg <-apply(all_map,1, mean)
all_map_avg<- round(all_map_avg, digits=4)

all_map_data<-data.frame(Id=all$Id,score=all_map_avg)

su_map = left_join(su, all_map_data)
sf::st_write(su_map, "salida/su_apa_mean_sd.shp", append = FALSE)
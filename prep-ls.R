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

memory.limit(800000)

#path to SAGA executable
env <- rsaga.env(r'(C:\Users\mreyes.AMBIENTE\saga-8.4.1_x64)')
#env <- rsaga.env(r'(C:\Users\ernes\saga-8.4.1_x64)')

#number of replications
n <- 10

##Read data

#folder name
area <- "els"

#Slope units
su <- read_sf("input/su/SU_ELS_LCC.shp") 
#Landslides
# ls_pol <- read_sf("input/landslides/san_vicente.shp") %>%
#   dplyr::filter(Evento != "Sismo") %>% st_transform(crs(su))
ls_point <- read_sf("input/landslides/landslides_test_lcc.shp") 

# #Plot of su and landslides
# ggplot()+geom_sf(data = su, fill=NA)+geom_sf(data=ls_pol, fill="red", color=NA)+theme_bw()

#DEM
#dem <- rast(paste0("input/continuous/", area, "/dem.tif"))

#Checks if at least one landslide is inside a slope unit
su$frane <- lengths(st_intersects(su, ls_point)) > 0
su$frane <- as.integer(ifelse(su$frane == "TRUE", 1, 0))

# #Landslides as polygons
# #Checks if polygons area is at least 0.1% of a slope unit as used by Goddard et al (2023) based on the 
# #results by Jacobs et al. (2023)
# su <- su %>% dplyr::mutate(area_su = st_area(su))
# 
# #https://gis.stackexchange.com/a/397962/36679
# su_ls <- st_intersection(ls_pol,su)
# 
# su_ls <- su_ls %>%
#   dplyr::mutate(area_ls = st_area(su_ls)) %>% 
#   st_drop_geometry() %>%
#   group_by(DN) %>% 
#   summarise(area_ls = sum(area_ls))
#   
# su_ls <- left_join(su, su_ls) %>%
#   st_drop_geometry() %>%
#   dplyr::mutate(ls_area_su = area_ls/area_su*100) %>%
#   drop_units() %>%
#   dplyr::mutate(su_ls = ifelse(ls_area_su >= 0.1 & !is.na(ls_area_su), 1,0)) %>%
#   dplyr::select(DN, su_ls)
#   
# su_ls <- left_join(su, su_ls)
# 
#Plot of su and landslides
#ggplot()+geom_sf(data = su_ls, aes(fill=factor(su_ls)))+geom_sf(data=ls_pol, fill="red", color=NA)+theme_bw()
ggplot()+geom_sf(data = su, aes(fill=factor(frane)),linewidth=0.05)+theme_bw() 
 
# writeRaster(dem, paste0("input/continuous/", area, "/dem.sdat"), overwrite=TRUE)
# 
# rsaga.slope.asp.curv(in.dem = paste0("input/continuous/", area, "/dem.sdat"), 
#                      out.slope = paste0("input/continuous/", area, "/slope"), unit.slope = 1,
#                      out.cprof = paste0("input/continuous/", area, "/cprof"), out.cplan = paste0("input/continuous/", area, "/cplan"),
#                      out.aspect = paste0("input/discrete/", area, "/asp.dat"), unit.aspect = 1,
#                      method = "poly2zevenbergen",
#                      env = env)
# 
# rsaga.geoprocessor(lib = "ta_morphometry", 
#                    module = "TPI Based Landform Classification",
#                    param = list(DEM = paste0("input/continuous/", area, "/dem.sdat"), 
#                                 LANDFORMS = paste0("input/discrete/", area, "/lcl.sdat"),
#                                 RADIUS_A_MIN = 0, 
#                                 RADIUS_A_MAX = 100,
#                                 RADIUS_B_MIN = 0,
#                                 RADIUS_B_MAX = 1000),
#                                 env = env)
# 
# asp <- rast(paste0("input/discrete/", area, "/asp.sdat"))
# m <- c(0, 22.5, 1,
#        22.5, 67.5, 2,
#        67.5, 112.5, 3,
#        112.5, 157.5, 4,
#        157.5, 202.5, 5,
#        202.5, 247.5, 6,
#        247.5, 292.5, 7,
#        292.5, 337.5, 8,
#        337.5, 360, 1)
# rcl <- matrix(m, ncol=3, byrow=TRUE)
# asp <- terra::classify(asp, rcl)
# writeRaster(asp, paste0("input/discrete/", area, "/asp.sdat"), overwrite=TRUE)

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

# vif_test <- subset( all, select = -c(Id) )
# sink("salida/vif.txt")  #Escribe la salida de VIF a vif.txt
# vif(vif_test)
# sink() 

# ##########################
# #Maxent modeling
# ##########################
# all_maxent <- all %>%
#   dplyr::rename(RV=frane) %>%
#   dplyr::select(-c(DN, geometry)) %>%
#   mutate_at(c('RV'), ~na_if(., 0))
# 
# all_maxentDVs <- deriveVars(all_maxent, 
#                            transformtype = c("L","M","D","B"))
# 
# all_maxentDVselect <- selectDVforEV(all_maxentDVs$dvdata, alpha = 0.001, quiet = TRUE)
# all_maxentEVselect <- selectEV(all_maxentDVselect$dvdata, alpha = 0.001, interaction = TRUE)
# 
# all_maxent_model <- chooseModel(all_maxentDVselect$dvdata, 
#                                 formula(reformulate(names(all_maxentDVselect$dvdata)[-1], "RV")))
# 
# all_m <- all
# 
# all_maxent_modelPreds <- projectModel(model = all_maxent_model,
#                                transformations = all_maxentDVs$transformations,
#                                data = all_m)
# 
# all_maxentfinal <- left_join(su, all_maxent_modelPreds$output)
# sf::st_write(all_maxentfinal, "salida/su_els_maxent.shp", append = FALSE)

##########################
#MARS modeling
##########################

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

####Cargar datos de calibraci?n y validaci?n previos####

#Usar las l?neas a continuaci?n s?lo se usa si se quieren usar datos previos de calibraci?n y validaci?n 
#para comparar variables usando los mismos set de datos

#all_cal14 <- readRDS("salida/all_cal14.rds")
#all_val14 <- readRDS("salida/all_val14.rds")

####model_all####
mars_all<-NULL
all_predict<-matrix(nrow=nrow(all_val14[[i]]), ncol = n)
AUC_all<-matrix(nrow=1, ncol=n)

####Regression analysis
####Cambiar esta linea si se cambian variables####

for(i  in 1:n){
  mars_all[[i]]<- earth (reformulate(names(dplyr::select(all, -c(DN, frane, geometry))), "frane"), data = all_cal14[[i]], trace = 1, degree=1,  glm=list(family=binomial), Scale.y=FALSE)}

for(i  in 1:n){
  all_predict[,i]<-predict(mars_all[[i]], all_val14[[i]], type=c("response"))
}

summary_earth <- lapply(1:n, function(x) capture.output(summary(mars_all[[x]]), file = "salida/summary_earth.csv", append = TRUE))

#################################
####ROC curve validation data####
#################################

#ROC curve for validation data
for(i  in 1:n){
  all_val14[[i]]$frane<-as.factor(all_val14[[i]]$frane)
}

for(i  in 1:n){
  AUC_all[,i]<-auc(all_val14[[i]]$frane,all_predict[,i])
}

write.csv(AUC_all, "salida/AUC_all_val.csv")

casi_all<-NULL

for(i  in 1:n){
  casi_all[[i]]<- data.frame(all_val14[[i]]$frane)
}

score_all<-NULL

for(i  in 1:n){
  casi_all[[i]]<-as.vector(casi_all[i])
  score_all[[i]]<-as.vector(all_predict[,i])
}


xrocall<-NULL

for(i in 1:n){
  xrocall[[i]]<-data.frame(score= score_all[[i]],response=casi_all[[i]])
}

for(i in 1:n){
  colnames(xrocall[[i]])<- c("score", "response")
}

roc_all<-NULL
for(i in 1:n){
  roc_all[[i]]<-roc(xrocall[[i]]$response ~ xrocall[[i]]$score, xrocall)}

youdenall<-NULL
for(i in 1:n){youdenall[[i]]<- min(coords(roc_all[[i]], "b", ret="t", best.method="youden"))}

youdenall_avg<-mean(youdenall)
youdenall_avg<- round(youdenall_avg, digits=3)
AUC_all_avg<- mean(AUC_all)
AUC_all_avg <- round(AUC_all_avg, digits = 3)

#################################
####AUC function####
#################################

auc <- function(vec, name){
  
  #ROC curve 
  vec %>% map_depth(., 1,~.x %>% mutate(across(c("frane"), as.factor)))
  
  for(i  in 1:n){
    AUC_all[,i]<-auc(vec[[i]][["frane"]],all_predict[,i])
  }
  
  write.csv(AUC_all, paste0("salida/AUC_", name))
  AUC_all_avg<- round(mean(AUC_all), digits = 3)
  return(AUC_all_avg)
}

auc(all_val14, "mars_random_val.csv")

auc <- function(vec, name){
  
  #ROC curve 
  all_cal14 <- vec %>% 
    map_depth(., 1,~.x %>% mutate(across(c("frane"), as.factor)))
}

######confusion#####
#COnfusion matrix for validation data
casiallc<-matrix(nrow = nrow(all_val14[[1]]), ncol=n)

for(i  in 1:n){
  
  casiallc[,i]<-as.vector(all_val14[[i]]$frane)
}

roc_all_m<-NULL
casesall<-as.vector(casiallc)
scoreall<-as.vector(all_predict)
dataall<-data.frame(score=scoreall,response=casesall)
roc_all_m<-roc(dataall$response ~ dataall$score, dataall)
youdenall_m<- min(coords(roc_all_m, "b", ret="t", best.method="youden"))
youdenall_m<- round(youdenall_m, digits=3)
auc_all_m<-round(roc_all_m$auc, digits = 3)
plotrocall<-plot.roc(roc_all[[1]], xlab="1-Specificity",ylab="Sensibility", main="ROC curve validation data", col= "grey", lwd = 0.1,asp = NA)
text(x=0.8, y= 1, paste("Cut-off =", youdenall_m), cex=1, col= "dark green")
text(x=0.2, y= 0.1, paste("AUC =",auc_all_m), cex=1, col= "red")

for (i in 2:n) {
  plot.roc(roc_all[[i]], add = TRUE,  col ='grey', lwd = 0.1)
}

plot.roc(roc_all_m, legacy.axes = TRUE, col= "red", lwd = 0.1, add=TRUE)

dev.print(pdf, 'salida/roc_curve_val.pdf')
dev.print(png, 'salida/roc_curve_val.png',width = 528, height = 553)

dataall$response<-as.factor( ifelse(dataall$response== "0", "0", "1") )
respall <-as.factor( ifelse(dataall$score< youdenall_m, "0", "1") )
conf_all <- confusionMatrix(respall, dataall$response, positive="1")

sink("salida/conf_matrix_val.txt")  #Escribe la matriz de confusion a conf_matrix.txt
printconfall<-print(conf_all, digits=max(3, getOption("digits") - 3 ), printStats= TRUE)
sink() 


#################################
####ROC curve calibration data####
#################################

#ROC curve for calibration data
all_predict<-matrix(nrow=nrow(all_cal14[[i]]), ncol = n)
AUC_all<-matrix(nrow=1, ncol=n)

for(i  in 1:n){
  all_predict[,i]<-predict(mars_all[[i]], all_cal14[[i]], type=c("response"))
}

for(i  in 1:n){
  all_cal14[[i]]$frane<-as.factor(all_cal14[[i]]$frane)
}

for(i  in 1:n){
  AUC_all[,i]<-auc(all_cal14[[i]]$frane,all_predict[,i])
}

write.csv(AUC_all, "salida/AUC_all_cal.csv")


casi_all<-NULL

for(i  in 1:n){
  casi_all[[i]]<- data.frame(all_cal14[[i]]$frane)
}


score_all<-NULL

for(i  in 1:n){
  casi_all[[i]]<-as.vector(casi_all[i])
  score_all[[i]]<-as.vector(all_predict[,i])
}


xrocall<-NULL

for(i in 1:n){
  xrocall[[i]]<-data.frame(score= score_all[[i]],response=casi_all[[i]])
}

for(i in 1:n){
  colnames(xrocall[[i]])<- c("score", "response")
}

roc_all<-NULL
for(i in 1:n){
  roc_all[[i]]<-roc(xrocall[[i]]$response ~ xrocall[[i]]$score, xrocall)}

youdenall<-NULL
for(i in 1:n){youdenall[[i]]<- min(coords(roc_all[[i]], "b", ret="t", best.method="youden"))}

youdenall_avg<-mean(youdenall)
youdenall_avg<- round(youdenall_avg, digits=3)
AUC_all_avg<- mean(AUC_all)
AUC_all_avg <- round(AUC_all_avg, digits = 3)

######confusion#####
#COnfusion matrix for calibration data

casiallc<-matrix(nrow = nrow(all_cal14[[1]]), ncol=n)

for(i  in 1:n){
  
  casiallc[,i]<-as.vector(all_cal14[[i]]$frane)
}

roc_all_m<-NULL
casesall<-as.vector(casiallc)
scoreall<-as.vector(all_predict)
dataall<-data.frame(score=scoreall,response=casesall)
roc_all_m<-roc(dataall$response ~ dataall$score, dataall)
youdenall_m<- min(coords(roc_all_m, "b", ret="t", best.method="youden"))
youdenall_m<- round(youdenall_m, digits=3)
auc_all_m<-round(roc_all_m$auc, digits = 3)
plotrocall<-plot.roc(roc_all[[1]], xlab="1-Specificity",ylab="Sensibility", main="ROC curve calibration data", col= "grey", lwd = 0.1,asp = NA)
text(x=0.8, y= 1, paste("Cut-off =", youdenall_m), cex=1, col= "dark green")
text(x=0.2, y= 0.1, paste("AUC =",auc_all_m), cex=1, col= "red")

for (i in 2:n) {
  plot.roc(roc_all[[i]], add = TRUE,  col ='grey', lwd = 0.1)
}

plot.roc(roc_all_m, legacy.axes = TRUE, col= "red", lwd = 0.1, add=TRUE)

dev.print(pdf, 'salida/roc_curve_cal.pdf')
dev.print(png, 'salida/roc_curve_cal.png',width = 528, height = 553)

dataall$response<-as.factor( ifelse(dataall$response== "0", "0", "1") )
respall <-as.factor( ifelse(dataall$score< youdenall_m, "0", "1") )
conf_all <- confusionMatrix(respall, dataall$response, positive="1")

sink("salida/conf_matrix_cal.txt")  #Escribe la matriz de confusion a conf_matrix.txt
printconfall<-print(conf_all, digits=max(3, getOption("digits") - 3 ), printStats= TRUE)
sink() 

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
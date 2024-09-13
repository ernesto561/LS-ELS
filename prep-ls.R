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

memory.limit(800000)

#path to SAGA executable
env <- rsaga.env(r'(C:\Users\mreyes.AMBIENTE\saga-8.4.1_x64)')
#env <- rsaga.env(r'(C:\Users\ernes\saga-8.4.1_x64)')

#number of replications
n <- 10

#Slope units
su <- read_sf("input/su/su_san_vicente.shp") 
#Landslides
ls <- read_sf("input/landslides/san_vicente.shp") %>%
  dplyr::filter(Evento != "Sismo")

#DEM
dem <- rast("input/continuous/dem.tif")

#Checks if at least one landslide is inside a slope unit
su$frane <- lengths(st_intersects(su, ls)) > 0
su$frane <- as.integer(su$frane == "TRUE")

writeRaster(dem, "input/continuous/dem.sdat", overwrite=TRUE)

rsaga.slope.asp.curv(in.dem = "input/continuous/dem.sdat", 
                     out.slope = "input/continuous/slope", unit.slope = 1,
                     out.cprof = "input/continuous/cprof", out.cplan = "input/continuous/cplan",
                     out.aspect = "input/discrete/asp.dat", unit.aspect = 1,
                     method = "poly2zevenbergen",
                     env = env)

asp <- rast("input/discrete/asp.dat")
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

#Continuous variables
vars_cont <- map(list.files("input/continuous/", pattern="*.sdat$", full.names = T), rast)
su_vars_cont <- exact_extract(rast(vars_cont), su, c('median', 'stdev'))

#Discrete variables
vars_disc <- map(list.files("input/discrete/", pattern="*.sdat$", full.names = T), rast)
su_vars_disc <- exact_extract(rast(vars_disc), su, 'majority')

su_model <- data.frame(cbind(su, su_vars_cont)) |> dplyr::select(-c(geometry, gridcode)) 

all<-data.frame(su_model)

#Removes NoData Values
all=na.exclude(all)

rownames(all) <- NULL

vif_test <- subset( all, select = -c(Id) )
sink("salida/vif.txt")  #Escribe la salida de VIF a vif.txt
vif(vif_test)
sink() 

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

####An?lisis de regresi?n####
####Cambiar esta linea si se cambian variables####

for(i  in 1:n){
  mars_all[[i]]<- earth (reformulate(names(dplyr::select(all, -c(Id, frane))), "frane"), data = all_cal14[[i]], trace = 1, degree=1,  glm=list(family=binomial), Scale.y=FALSE)}

for(i  in 1:n){
  all_predict[,i]<-predict(mars_all[[i]], all_val14[[i]], type=c("response"))
}

for(i  in 1:n){
  all_val14[[i]]$frane<-as.factor(all_val14[[i]]$frane)
}

for(i  in 1:n){
  AUC_all[,i]<-auc(all_val14[[i]]$frane,all_predict[,i])
}

write.csv(AUC_all, "salida/AUC_all_val.csv")
summary_earth <- lapply(1:n, function(x) capture.output(summary(mars_all[[x]]), file = "salida/summary_earth.csv", append = TRUE))


#################
####ROC curve####
#################

#Curva ROC creada para los datos de validaci?n

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


######confusion#####
#Matriz de confusi?n creada para los datos de validaci?n
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


#################
####ROC curve####
#################

#Curva ROC creada para los datos de calibraci?n
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
#Matriz de confusi?n creada para los datos de calibraci?n

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

su_map = left_join(su_els, all_map_data)
sf::st_write(su_map, "salida/su_apa_mean_sd.shp", append = FALSE)
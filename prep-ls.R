library(RSAGA)
library(exactextractr)
library(tidyverse)
library(sf)
library(terra)


su_els <- read_sf("input/su_els_apaneca.shp") 
ls <- read_sf("input/landslides_test.shp")

su_els$contains_point <- lengths(st_intersects(su_els, ls)) > 0

plot(su_els["contains_point"])
plot(ls, add=TRUE)

su_ls <- sum(su_els$contains_point == 1)
su_nols <- length(su_els$contains_point)-su_ls

#path to SAGA executable
#env <- rsaga.env(r'(C:\Users\mreyes.AMBIENTE\saga-8.4.1_x64)')
env <- rsaga.env(r'(C:\Users\ernes\saga-8.4.1_x64)')

dem <- rast("input/dem.tif") 

writeRaster(dem, "input/dem.sdat", overwrite=TRUE)

rsaga.slope.asp.curv(in.dem = "input/dem.sdat", out.slope = "input/slope",
                     out.cprof = "input/cprof", out.cplan = "input/cplan",
                     method = "poly2zevenbergen",
                     env = env)



vars <- map(list.files("input/test", pattern="*.sdat$", full.names = T), rast)

su_vars <- exact_extract(rast(vars), su_els, 'mean')

brazil <- cbind(su_els, exact_extract(dem, su_els, 'mean'))

mean_var(su_els, vars[4][[1]])

mean_vars <- pmap(list(vars, su_els), terra::extract)

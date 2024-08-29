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

 
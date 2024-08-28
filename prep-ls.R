library(tidyverse)
library(sf)

su_els <- read_sf("input/su_els_apaneca.shp")
ls <- read_sf("input/landslides_test.shp")

su_els$contains_point <- lengths(st_intersects(su_els, ls)) > 0


## Script takes two shapefiles with polygons of different levels of granularity and aggregates from finer to coarser
## specifically, this script aggregates from census tract to nyc neighborhoods

library(data.table)
library(dplyr)
library(rgdal)
library(ggplot2)
library(maptools)
library(sp)
library(rgeos)
library(geojsonio)
library(PBSmapping)

## load functions
source("/Users/billbachrach/Dropbox (hodgeswardelliott)/Team_NYC/Bill Bachrach/useful functions/HWE_FUNCTIONS.r")
source("/Users/billbachrach/Dropbox (hodgeswardelliott)/Team_NYC/Bill Bachrach/useful functions/-.R")
source("/Users/billbachrach/Dropbox (hodgeswardelliott)/Team_NYC/Bill Bachrach/useful functions/useful minor functions.R")

## names of the polygons in both the fine (small) and coarse (large) shapefiles
region.large <- "neighborhood"
region.small <- "BoroCT2010"

## read in larger shapefiles
largeshape.url <- "http://data.beta.nyc//dataset/0ff93d2d-90ba-457c-9f7e-39e47bf2ac5f/resource/35dd04fb-81b3-479b-a074-a27a37888ce7/download/d085e2f8d0b54d4590b1e7d1f35594c1largecitiesnycneighborhoods.geojson"
large.map <- geojson_read(as.location(largeshape.url),
                          method="local",
                          what="sp")
large.map <- large.map %>% spTransform(CRS("+proj=longlat +datum=WGS84"))

## set region of interest to more generic name
names(large.map) <- gsub(region.large,"region.large",names(large.map))

large.f <- large.map %>% 
  fortify(region="region.large")
nyc.large <- merge(large.f
                   ,large.map@data
                   ,by.x="id"
                   ,by.y="region.large")


## read in smaller shapefiles 
small.map <- readOGR("/Users/billbachrach/Dropbox (hodgeswardelliott)/Team_NYC/Bill Bachrach/Data Sources/Census/Shapefiles/nyct2010_16d/nyct2010.shp"
                     , layer="nyct2010")
small.map <- small.map %>% spTransform(CRS("+proj=longlat +datum=WGS84"))
## set region of interest to more generic name
names(small.map) <- gsub(region.small,"region.small",names(small.map))
small.f <- small.map %>% 
  fortify(region="region.small")
nyc.small <- merge(small.f
                   ,small.map@data
                   ,by.x="id"
                   ,by.y="region.small")


## Determine which smaller polygons fit entirely into the larger polygons
dat <- nyc.small[,c("long","lat","id")]
colnames(dat) <- c("Longitude","Latitude","id")
coordinates(dat) <- ~ Longitude + Latitude
proj4string(dat) <- proj4string(large.map)
assign.vec <- over(dat, large.map)[,"region.large"]
nyc.small[,"region.large"] <- assign.vec


## moving on to smaller polygons which do not have an exact match 
na.obs <- which(is.na(nyc.small[,"region.large"]))
id.vec <- unique(nyc.small[,"id"])

## for each smaller polygon, count how larger polygons its points lie in
## number of region.larges each minor polygon could belong to
rl_assign.list <- lapply(id.vec, function(x){
  out <- nyc.small %>% filter(id==x & !is.na(region.large)) %>% 
    group_by(region.large) %>% 
    summarize(count=n()) %>% 
    mutate(prop=count/sum(count))
  return(out)
}
)

# multi_rl <- which(unlist(lapply(rl_assign.list, function(x) nrow(x)))>1)

## creating a vector which identifies region.large if either there is a single region.large or 90% of points lie within a region.large
assign.vec <- unlist(lapply(rl_assign.list, function(x){
  out <- NA
  if(nrow(x)==1){
    out <- as.character(
      as.data.frame(x)[1,"region.large"]
    )
  } else{
    major <- (x[,"prop"]) >= 0.9
    if(sum(major)==1){
      out <- as.character(
        as.data.frame(x)[major,"region.large"]
      )
    }
  }
  return(out)
}
)
)


## which of those are still NA
na.rls <- which(is.na(assign.vec))

## put the data into a dataframe key 
region.large.key <- as.data.frame(cbind(id.vec,
                                        assign.vec),
                                  stringsAsFactors=F)
colnames(region.large.key) <- c("id","region.large")


## following for loop tests intersesmallions of the smalls which have not been assigned and possible match region.larges
## the region.large with the most overlapping area gets assigned to the small

max_rl.out <- list()

for(i in 1:length(na.rls)){
  itemno <- which(small.map$region.small == region.large.key[na.rls[i],1])
  nbrhd.tmp <- as.character(as.data.frame(rl_assign.list[[na.rls[i]]])[,"region.large"])
  large.no <- which(large.map$region.large%in% nbrhd.tmp)
  large.rls <- as.character(large.map$region.large[large.no])
  
  area.out <- list()
  for(j in 1:length(large.no)){
    # assign("last.warning",NULL,envir=baseenv())
    x <- large.no[j]
    ps.1 <- as.PolySet(
      SpatialPolygons2PolySet(large.map[x,])
    )
    ps.2 <- as.PolySet(
      SpatialPolygons2PolySet(small.map[itemno,])
    )
    
    ps.3 <- joinPolys(ps.1,ps.2)
    
    ps.3 <- try(PolySet2SpatialPolygons(joinPolys(ps.1,ps.2)),silent=T)
    if(class(ps.3)!="try-error"){
      area <- ps.3@polygons[[1]]@area
      
      area.out[[j]] <- area } else {
        area.out[[j]] <- NA
        cat("try-error at iteration",i,"\n")
      }
  }
  
  if(sum(is.na(area.out))>0){
    drops <- which(is.na(area.out))
    area.out[drops] <- NULL
    large.rls <- large.rls[-drops]
  }
  
  if(length(area.out)>0){
    # cat("pulling area from results\n")
    small.area <- small.map@polygons[[itemno]]@area
    area.out <- unlist(area.out)
    max.rl <- large.rls[which.max(area.out)]
    max_rl.out[[i]] <- max.rl
    itertell.fun(i,25)
  } else {
    max_rl.out[[i]] <- NA
    cat("no intersecting area in iteration",i,"\n")
  }
}

region.large.key[na.rls,"region.large"] <- unlist(max_rl.out)

## put back in the original region names
colnames(region.large.key) <- c(region.small,region.large)

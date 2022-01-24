install.packages("xlsx")
install.packages("readxl")
install.packages("tidyr")
install.packages("maptools")
install.packages("geosphere")
install.packages("rgl")
library(tidyr)
library(xts)
library(sp)
library(spacetime)
library(xlsx)
library(readxl)
library(gstat)
library(sp)
library(spacetime)
library(raster)
library(rgdal)
library(rgeos) 
library(maptools)
library(lattice)
library(geosphere)
library(rgl)
setwd(choose.dir())


# read first sheet
my_data <- read_excel("WholeAlberta3.xlsx", sheet = 2)
sorteddata<-gather(my_data, na.rm = FALSE, convert = FALSE)
valuedata<-sorteddata$value
valuedata
class(valuedata)
windpwva<-data.frame(valuedata)

windflocations <- read_excel("metafin.xlsx", sheet = 1)
#windflocations <- read_excel("meta Alberta3.xlsx", sheet = 1)
summary(windflocations)
coordinates(windflocations)=~X+Y

crs3TM83 <- CRS("+proj=tmerc +lat_0=0 +lon_0=-114 +k=0.9999 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs", doCheckCRSArgs = FALSE)
proj4string(windflocations) <- crs3TM83

plot(windflocations)
plot(NewsiteSp,add=T,col="red")

Windsp <- SpatialPoints(windflocations@coords,crs3TM83)
Windsp
plot(Windsp)

timedf <- read_excel("time.xlsx", sheet = 2)
windfarmtime<-timedf$time
class(windfarmtime)
timewindf <-strptime(timedf$time,"%m/%d/%Y %H:%M:%S")
windtimeStamp<-as.POSIXlt(timewindf,tz="MST")
windtimeStamp[2]
timefinal<-as.POSIXct(timewindf,tz="MST")

windDF <- STFDF(Windsp,timefinal,data=windpwva) 

xs1 <- as(windDF, "Spatial") 
class(xs1)
xs1
attr(xs1, "time")
scales <- list(x = list(rot = 0))
stplot(windDF, mode = "xt", scales = scales, xlab ="windfarm id" )
##variogramSurface(model, dist_grid, covariance = FALSE)


var<-variogramST(valuedata~1,data = windDF,assumeRegular = F,na.omit=T)
plot(var,wireframe=T)
plot(var,map=F)
plot(var,map=T) 

pars.l <- c(sill.s = 0, range.s = 0, nugget.s = 0,sill.t = 0, range.t = 1, nugget.t = 0,sill.st = 0, 
                                                             
            range.st = 1, nugget.st = 0, anis = 1)

pars.u <- c(sill.s = 200, range.s = 80000, nugget.s = 100,sill.t = 200, range.t = 24, nugget.t = 100,
                                  
            sill.st = 200, range.st = 1000, nugget.st = 100,anis = 500) 


separable <- vgmST("separable", space=vgm(0.2,"Exp", 80000,0.1),time =vgm( 0.7,"Exp",5, 0.1), sill=28)

plot(var,separable,map=T)
plot(var,separable,map=F)

vvvgm<- fit.StVariogram(var, separable, fit.method=0)
Separvgm<-vvvgm
attr(vvvgm, "MSE")

Alberta <- readOGR(dsn = "C:/Users/msala/Desktop/wind canada", layer = "Alberta_CD2016")

Alb<- spTransform(Alberta,CRS("+init=epsg:3395"))

plot(Alb) 

tm.grid <- seq(as.POSIXct('2008-2-1 00:00 CET') ,as.POSIXct('2008-10-1 24:00 CET'),length.out=6) 

tm.grid

proj4string(sp.grid) <- crs3TM83

grid.ST <- STF(Windsp,tm.grid) 
grid.ST

pred <- krigeST(valuedata~1, data=windDF, modelList=vvvgm, newdata=grid.ST)

stplot(pred)

my_data2 <- read_excel("WholeAlberta3.xlsx", sheet = 3)
sorteddata2<-gather(my_data2, na.rm = FALSE, convert = FALSE)
valuedata2<-sorteddata2$value
valuedata2
class(valuedata2)
windpwva2<-data.frame(valuedata2)
windDF2 <- STFDF(Windsp,tm.grid,data=windpwva2) 
tm.grid

stplot(windDF2)
xs1 <- as(windDF, "Spatial") 
class(xs1)
xs1
attr(xs1, "time")
scales <- list(x = list(rot = 0))
stplot(windDF, mode = "xt", scales = scales, xlab ="windfarm id" )
predictkrig<-pred$var1.pred
predvalues<-data.frame(predictkrig)
observedvalues<-data.matrix(windpwva2, rownames.force = NA)

plot(observedvalues, predictkrig, main = NULL,
     xlab = "Observations", ylab = "predictions",
     pch = 19, frame = FALSE)
abline(lm(predictkrig ~ observedvalues ), col = "blue")
observedvalues<-data.matrix(windpwva2, rownames.force = NA)
error<-observedvalues-predictkrig
errordata<-data.frame(error)
errordf <- STFDF(Windsp,tm.grid,data=errordata)
stplot(errordf)
samplewind<-spsample(Windsp, 1000, type="regular", bb = bbox(Windsp))
class(samplewind)
plot(spsample(Windsp, 1000, type="regular", bb = bbox(Windsp)))
gridsample.ST <- STF(samplewind,tm.grid)
writeOGR(spdatarandom, dsn=getwd(),layer="Samp",driver="ESRI Shapefile")
predsample <- krigeST(valuedata~1, data=windDF, modelList=vvvgm, newdata=gridsample.ST)
predsample
stplot(predsample)
class(predsample)
predictkrig2<-predsample$var1.pred
ranwind<-data.frame(samplewind)
samplewind
subpred2<-predictkrig2[1:1008]
datapred2<-data.frame(subpred2)
spsamplewind<-SpatialPointsDataFrame(samplewind,datapred2)
spdatarandom
Albnew<-spTransform(Alb,crs3TM83)

plot(Albnew)
plot(samplewind, add=T, col="red", pch=20, cex=2)

write.csv(spsamplewind, file = "spsamplewind.csv")

prodSumModel <- vgmST("productSum",space = vgm(2, "Exp", 30000, 0.01), time = vgm(35, "Exp", 20, 0.8), k=0.05) 

plot(var,prodSumModel,map=F)
productsum2 <- fit.StVariogram(var, prodSumModel, fit.method=0)
attr(productsum2,"MSE")

#Automatic fit
prodSumModel_Vgm <- fit.StVariogram(var, prodSumModel,method = "L-BFGS-B",lower=pars.l)

plot(var, prodSumModel_Vgm,all=T,wireframe=T, zlim=c(0,30),  zlab=list("gamma", rot=45),  xlab=list("distance (m)", rot=40),  ylab=list("time lag (houres)", rot=-40),  scales=list(arrows=F, z = list(distance = 5)))

attr(prodSumModel_Vgm, "MSE")

SimplesumMetric <- vgmST("simpleSumMetric",space = vgm(2,"Exp", 30000, 0.01),
                         
      time = vgm(35,"Exp", 20, 0.8), joint = vgm(1,"Lin", 500, 0), nugget=1, stAni=500)
time.grid

SimplesumMetric_Vgm <- fit.StVariogram(var, SimplesumMetric,method = "L-BFGS-B",lower=pars.l)
SimplesumMetric2 <- fit.StVariogram(var, SimplesumMetric, fit.method=0)
plot(var,SimplesumMetric,map=F)
plot(var, SimplesumMetric_Vgm,all=T,wireframe=T, zlim=c(0,30),  zlab=list("gamma", rot=45),  xlab=list("distance (m)", rot=40),  ylab=list("time lag (houres)", rot=-40),  scales=list(arrows=F, z = list(distance = 5)))
attr(simplesumMetric_Vgm, "MSE")
attr(SimplesumMetric2, "MSE")

sumMetric <- vgmST("sumMetric", space = vgm(psill=7,"Exp", range=40000, nugget=1.5),
                   
                 time = vgm(psill=22,"Exp", range=10, nugget=0.9), 
                
                  joint = vgm(psill=1,"Lin", range=50, nugget=0.05), stAni=300)

sumMetric_Vgm <- fit.StVariogram(var, sumMetric, method="L-BFGS-B",lower=pars.l,upper=pars.u,tunit="hours")                   
sumMetric2<-fit.StVariogram(var, sumMetric, fit.method=0)
plot(var,sumMetric,map=F)
attr(sumMetric_Vgm, "MSE")
attr(sumMetric2, "MSE")
plot(var, sumMetric_Vgm,all=T,wireframe=T, zlim=c(0,30),  zlab=list("gamma", rot=45),  xlab=list("distance (m)", rot=40),  ylab=list("time lag (houres)", rot=-40),  scales=list(arrows=F, z = list(distance = 5)))

separable <- vgmST("separable",  space=vgm(0.2, "Exp", 80000, 0.1),  time =vgm(   0.7, "Exp",   5, 0.1), sill=28)

plot(var,separable,map=F)

 
metric <- vgmST("metric", joint = vgm(23,"Exp", 130000, 2), stAni=16600)

plot(var,metric,map=F)
plot(var,metric,map=T)
metricvgm<- fit.StVariogram(var, metric, fit.method=0)
attr(metricvgm, "MSE")
#autofit for metric
metric_Vgm2 <- fit.StVariogram(var, metric, method="L-BFGS-B",lower=pars.l)
attr(metric_Vgm2, "MSE")
plot(var,metric_Vgm2,map=F)
plot(var,metric_Vgm2,map=T)
plot(var, metric_Vgm2,all=T,wireframe=T, zlim=c(0,30),  zlab=list("gamma", rot=45),  xlab=list("distance (m)", rot=40),  ylab=list("time lag (houres)", rot=-40),  scales=list(arrows=F, z = list(distance = 5)))
#autofit for separable
separable_Vgm2 <- fit.StVariogram(var, separable, fit.method=11,method="L-BFGS-B", stAni=5, lower=pars.l,upper=pars.u)
attr(separable_Vgm2, "MSE")
#plot all
plot(var,list(metric_Vgm2, sumMetric_Vgm, prodSumModel_Vgm,Separvgm,SimplesumMetric2),all=T,wireframe=T, zlim=c(0,30),  zlab=list("gamma", rot=45),  xlab=list("distance (m)", rot=40),  ylab=list("time lag (houres)"
        , rot=-40),  scales=list(arrows=F, z = list(distance = 5))) 
plot(var,list(metric_Vgm2, sumMetric_Vgm, prodSumModel_Vgm,Separvgm,SimplesumMetric2),all=T)
plot(var,list(metric_Vgm2, sumMetric_Vgm, prodSumModel_Vgm,Separvgm,SimplesumMetric2),map=F,all=T)

#separable_Vgm2 <- fit.StVariogram(var, separable, fit.method=11,method="L-BFGS-B", stAni=5, lower=pars.l,upper=pars.u)
#attr(separable_Vgm2, "MSE")
pred <- krigeST(PPB~1, data=timeDF, modelList=sumMetric_Vgm, newdata=grid.ST)
summetvgm <- fit.StVariogram(var, prodSumModel, fit.method=0)
attr(summetvgm, "MSE")
plot(var,summetvgm,map=F)
predmetr <- krigeST(valuedata~1, data=windDF, modelList=metricvgm, newdata=grid.ST)
pred <- krigeST(valuedata~1, data=windDF, modelList=vvvgm, newdata=grid.ST)
pred3 <- krigeST(valuedata~1, data=windDF, modelList=sumMetric_Vgm, newdata=grid.ST)

predmetric<-krigeST(valuedata~1, data = windDF, newdata = grid.ST, metric_Vgm2, nmax = 20)

predmetric
NewsiteSp<-read_excel("newsite.xlsx", sheet = 1)
coordinates(NewsiteSp)=~X+Y
proj4string(NewsiteSp) <- crs3TM83
newwfarm<- SpatialPoints(NewsiteSp@coords,crs3TM83)
stgridfornewsite<-STF(newwfarm,timefinal)
stgridfornewsite
predmetnew<-krigeST(valuedata~1, data = windDF, newdata = stgridfornewsite, metric_Vgm2, nmax = 20)
predmetnew[1]
prenewkrig<-predmetnew$var1.pred
prednewvalues<-data.frame(prenewkrig)
####observedvalues<-data.matrix(windpwva2, rownames.force = NA)
observed1889data<-read_excel("Windfarm1889.xlsx", sheet = 1)
valueobserved1889<-observed1889data$value
##observed1889<-data.matrix(observed1889data, rownames.force = NA)
class(observed1889data)
subpred1889<-valueobserved1889[1:8784]
plot(subpred1889, prenewkrig, main = "Cross validation plot ",
     xlab = "Observations", ylab = "predictions",
     pch = 19, frame = FALSE)
summary(prenewkrig)


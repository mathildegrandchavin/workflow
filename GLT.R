
#data temperature


getwd()
setwd("D:/quercus/Utilisateurs/mgrandchavin/climate data")

#install packages
install.packages("rgdal")

#whivh packages are alreay installed.

a<-installed.packages()
packages<-a[,1] 
is.element(c("sp","rgdal","raster","dismo","maptools"), packages)

install.packages("dismo")

#load packages
library(raster)
library(rgdal)
library(sp)
library(dismo)
library(maptools)
##upload GLT T max data

GLT_Tmax_1<-raster("Tmax/Tmax_0320_1.tif")

GLT_Tmax_1<-raster("Tmax/Tmax_0320_1.tif")
GLT_Tmax_2<-raster("Tmax/Tmax_0320_2.tif")
GLT_Tmax_3<-raster("Tmax/Tmax_0320_3.tif")
GLT_Tmax_4<-raster("Tmax/Tmax_0320_4.tif")
GLT_Tmax_5<-raster("Tmax/Tmax_0320_5.tif")
GLT_Tmax_6<-raster("Tmax/Tmax_0320_6.tif")
GLT_Tmax_7<-raster("Tmax/Tmax_0320_7.tif")
GLT_Tmax_8<-raster("Tmax/Tmax_0320_8.tif")
GLT_Tmax_9<-raster("Tmax/Tmax_0320_9.tif")
GLT_Tmax_10<-raster("Tmax/Tmax_0320_10.tif")
GLT_Tmax_11<-raster("Tmax/Tmax_0320_11.tif")
GLT_Tmax_12<-raster("Tmax/Tmax_0320_12.tif")
GLT_Tmax_13<-raster("Tmax/Tmax_0320_13.tif")

extent(GLT_Tmax_1)

#upload GLT T min data

setwd

GLT_Tmin_1<-raster("Tmin/Tmin_0320_1.tif")
GLT_Tmin_2<-raster("Tmin/Tmin_0320_2.tif")
GLT_Tmin_3<-raster("Tmin/Tmin_0320_3.tif")
GLT_Tmin_4<-raster("Tmin/Tmin_0320_4.tif")
GLT_Tmin_5<-raster("Tmin/Tmin_0320_5.tif")
GLT_Tmin_6<-raster("Tmin/Tmin_0320_6.tif")
GLT_Tmin_7<-raster("Tmin/Tmin_0320_7.tif")
GLT_Tmin_8<-raster("Tmin/Tmin_0320_8.tif")
GLT_Tmin_9<-raster("Tmin/Tmin_0320_9.tif")
GLT_Tmin_10<-raster("Tmin/Tmin_0320_10.tif")
GLT_Tmin_11<-raster("Tmin/Tmin_0320_11.tif")
GLT_Tmin_12<-raster("Tmin/Tmin_0320_12.tif")
GLT_Tmin_13<-raster("Tmin/Tmin_0320_13.tif")

#upload mesoTair_p0320_vOK

mesoTair_Tmin_1<-raster("mesoTair_p0320_vOK/Tmin/Tmin_0320_1.tif")

mesoTair_Tmin_2<-raster("mesoTair_p0320_vOK/Tmin/Tmin_0320_2.tif")
mesoTair_Tmin_3<-raster("mesoTair_p0320_vOK/Tmin/Tmin_0320_3.tif")
mesoTair_Tmin_4<-raster("mesoTair_p0320_vOK/Tmin/Tmin_0320_4.tif")
mesoTair_Tmin_5<-raster("mesoTair_p0320_vOK/Tmin/Tmin_0320_5.tif")
mesoTair_Tmin_6<-raster("mesoTair_p0320_vOK/Tmin/Tmin_0320_6.tif")
mesoTair_Tmin_7<-raster("mesoTair_p0320_vOK/Tmin/Tmin_0320_7.tif")
mesoTair_Tmin_8<-raster("mesoTair_p0320_vOK/Tmin/Tmin_0320_8.tif")
mesoTair_Tmin_9<-raster("mesoTair_p0320_vOK/Tmin/Tmin_0320_9.tif")
mesoTair_Tmin_10<-raster("mesoTair_p0320_vOK/Tmin/Tmin_0320_10.tif")
mesoTair_Tmin_11<-raster("mesoTair_p0320_vOK/Tmin/Tmin_0320_11.tif")
mesoTair_Tmin_12<-raster("mesoTair_p0320_vOK/Tmin/Tmin_0320_12.tif")
mesoTair_Tmin_13<-raster("mesoTair_p0320_vOK/Tmin/Tmin_0320_13.tif")

stack_mesoTair<-stack(mesoTair_Tmin_1,mesoTair_Tmin_2,mesoTair_Tmin_3,mesoTair_Tmin_4,mesoTair_Tmin_5,mesoTair_Tmin_6,mesoTair_Tmin_7,mesoTair_Tmin_8,mesoTair_Tmin_9,mesoTair_Tmin_10,mesoTair_Tmin_11,mesoTair_Tmin_12,mesoTair_Tmin_13)
summary(stack_mesoTair)

memory.limit(size=5000)

matrixw<-as.matrix(stack_mesoTair)
summary(matrixw)

pts.mesoTair_Tmin_1<-getValues(mesoTair_Tmin_1)
pts.mesoTair_Tmin_2<-getValues(mesoTair_Tmin_1)
pts.mesoTair_Tmin_3<-getValues(mesoTair_Tmin_1)
pts.mesoTair_Tmin_4<-getValues(mesoTair_Tmin_1)
pts.mesoTair_Tmin_5<-getValues(mesoTair_Tmin_1)
pts.mesoTair_Tmin_6<-getValues(mesoTair_Tmin_1)
pts.mesoTair_Tmin_7<-getValues(mesoTair_Tmin_1)
pts.mesoTair_Tmin_8<-getValues(mesoTair_Tmin_1)
pts.mesoTair_Tmin_9<-getValues(mesoTair_Tmin_1)
pts.mesoTair_Tmin_10<-getValues(mesoTair_Tmin_1)
pts.mesoTair_Tmin_11<-getValues(mesoTair_Tmin_1)
pts.mesoTair_Tmin_12<-getValues(mesoTair_Tmin_1)
pts.mesoTair_Tmin_13<-getValues(mesoTair_Tmin_1)
?

summary(pts.mesoTair_Tmin_13)
df<-rbind(pts.mesoTair_Tmin_1,pts.mesoTair_Tmin_2,pts.mesoTair_Tmin_3,pts.mesoTair_Tmin_4,pts.mesoTair_Tmin_5,pts.mesoTair_Tmin_6,pts.mesoTair_Tmin_7,pts.mesoTair_Tmin_8,pts.mesoTair_Tmin_9,pts.mesoTair_Tmin_10)


list.mesoTair_Tmin<-c(mesoTair_Tmin_1,mesoTair_Tmin_2,mesoTair_Tmin_3,mesoTair_Tmin_4,mesoTair_Tmin_5,mesoTair_Tmin_6,mesoTair_Tmin_7,mesoTair_Tmin_8,mesoTair_Tmin_9,mesoTair_Tmin_10,mesoTair_Tmin_11,mesoTair_Tmin_12,mesoTair_Tmin_13)
is.list(list.mesoTair_Tmin)


library(tibble)
list.mesoTair_Tmin<-as_tibble(list.mesoTair_Tmin)
library(foreach)

xx<-lapply(list.mesoTair_Tmin,getValues)

xx<-foreach (list.mesoTair_Tmin) %do% {
	
getValues(list.mesoTair_Tmin)}




summary(pts.mesoTair_Tmin_1)



#get a data.frame with : T min, T max, month

#get info
bbox(GLT_Tmin_1)
ncol(GLT_Tmin_2)
ncol(GLT_Tmax_1)
res(GLT_Tmax_1)

#plot the data
plot(GLT_Tmin_1)
str(GLT_Tmin_1)
projection(GLT_Tmin_1)
extent(GLT_Tmin_1)





##pairwise correlations between variables

par(mfrow=c(2,2))
plot(GLT_Tmin_1,GLT_Tmax_1)
plot(GLT_Tmin_1,GLT_Tmin_2)
plot(GLT_Tmin_2,GLT_Tmax_2)
plot(GLT_Tmax_2,GLT_Tmax_3)

#stack




empilement_Tmin<-stack(GLT_Tmin_1,GLT_Tmin_2,GLT_Tmin_3,GLT_Tmin_4,GLT_Tmin_5,GLT_Tmin_6,GLT_Tmin_7,GLT_Tmin_8,GLT_Tmin_9,GLT_Tmin_10,GLT_Tmin_11,GLT_Tmin_12,GLT_Tmin_13)
empilement_Tmax<-stack(GLT_Tmax_1,GLT_Tmax_2,GLT_Tmax_3,GLT_Tmax_4,GLT_Tmax_5,GLT_Tmax_6,GLT_Tmax_7,GLT_Tmax_8,GLT_Tmax_9,GLT_Tmax_10,GLT_Tmax_11,GLT_Tmax_12,GLT_Tmax_13)


#error: not the same extent
lapply(c(GLT_Tmax_1,GLT_Tmax_2,GLT_Tmax_3,GLT_Tmax_4,GLT_Tmax_5,GLT_Tmax_6,GLT_Tmax_7,GLT_Tmax_8,GLT_Tmax_9,GLT_Tmax_10,GLT_Tmax_11,GLT_Tmax_12,GLT_Tmax_13),extent)

##13 n'a pas le même extent



GLT_Tmax_13<-resample(GLT_Tmax_13,GLT_Tmax_11,method="ngb")
empilement_Tmax<-stack(GLT_Tmax_1,GLT_Tmax_2,GLT_Tmax_3,GLT_Tmax_4,GLT_Tmax_5,GLT_Tmax_6,GLT_Tmax_7,GLT_Tmax_8,GLT_Tmax_9,GLT_Tmax_10,GLT_Tmax_11,GLT_Tmax_12,GLT_Tmax_13)

extent(GLT_Tmax_12)

empilement_Tmin
empilement_Tmax

empilement_GLT<-stack(GLT_Tmin_1,GLT_Tmin_2,GLT_Tmin_3,GLT_Tmin_4,GLT_Tmin_5,GLT_Tmin_6,GLT_Tmin_7,GLT_Tmin_8,GLT_Tmin_9,GLT_Tmin_10,GLT_Tmin_11,GLT_Tmin_12,GLT_Tmin_13,GLT_Tmax_1,GLT_Tmax_2,GLT_Tmax_3,GLT_Tmax_4,GLT_Tmax_5,GLT_Tmax_6,GLT_Tmax_7,GLT_Tmax_8,GLT_Tmax_9,GLT_Tmax_10,GLT_Tmax_11,GLT_Tmax_12,GLT_Tmax_13)

plot(empilement_GLT)

##occurences

install.packages("data.table")
install.packages("bit64")
library(bit64)
library(data.table)
getwd()
setwd("S:/Utilisateurs/mgrandchavin/data_occ/occurences")
occurences<-fread("occurences.csv")
summary(occurences)
class(occurences)



SpOccurences<-occurences
summary(SpOccurences)

coordinates(SpOccurences)<-c("decimalLongitude","decimalLatitude")
summary(SpOccurences)
coordinates(SpOccurences)
class(SpOccurences)

str(SpOccurences)

##set a coordinate reference system
##we need to set the right CRS, do we need to ensure that the angles or areas are correctly represented?
projection(SpOccurences)<-"+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs"

plot(GLT_Tmin_1)
plot(SpOccurences)

##CLEAN DATA

#STEP 1: DATA WITH NO SPATIAL DATA (H)

output<-foreach { sp=PlantList
if data has no spatial data

return data frame (value = H)

#STEP 2: DUPLICATE DATA (G)

#1) TRUE REPLICATE (SAME COORDONATES)
#2) GEOGRAPHICAL REPLICATES

#STEP 3: PERFORM A ENVIRONMENTAL CONGRUENCE ANALYSIS: IS DATA IN SEA OR LAKES?

##To avoid such spurious effects we use the function nearestcell in package biogeo in R (Robertson et al. 2016), 
##which checks whether neighboring environmental cells have congruent (e.g. terrestrial) environments.
## If that is the case, the function assigns new coordinates to the record and creates a new field with the names ‘x.orginal’ and ‘y.original’
## to keep track of the transformation. Records with transformed coordinates are evaluated again for duplicates, as outlined previously in step 2.

#STEP 4: IDENTIFICATION OF POTENTIAL GEOLOCATION ISSUES 
#USE OF THE GTS DATABASE (BEECH ET AL 2017)
#CONSIDER COUNTRIES WHERE THE SP IS CONSIDERED INTRODUCED AND NATURALIZED THROUGH A GLOBAL INVASIVE SPECIES DATABASE (GISD)
#IF NOT IN THE COMBINED LIST OF COUNTRIES IN WHICH A SP IS KNOWN TO OCCUR: F

#STEP 5 : 

#MAT AND MET: WE USED THE DATA CLEANING WORKFLOW OF SERRA DIAZ ET AL (2018)

##clean the data

#install coordinateCleaner

install.packages("devtools")
library(devtools)

install_github("ropensci/CoordinateCleaner")

#set up the libraries and data

install.packages(c("countrycode","ggplot2"))
library(countrycode)
library(CoordinateCleaner)
library(dplyr)
library(ggplot2)
library(rgbif)
library(sp)


#data



setwd("C:/Users/grandchavin/desktop/data_occ")
occurences<-data.table::fread("occurences/occurences.csv")

#do that on the species directories

x<-foreach (sp="outputGBIFdir/...")




Hely<-occurences[4]
Pinus<-occurences[species=="Pinus uncinata"]
View(Pinus)
summary(Hely)
summary(Pinus)


dat<-Hely %>%
dplyr::select(species,decimalLongitude,decimalLatitude,countryCode,individualCount,gbifID,family,taxonRank,coordinateUncertaintyInMeters,year,basisOfRecord,institutionCode)

datPinus<-Pinus %>%
dplyr::select(species,decimalLongitude,decimalLatitude,countryCode,individualCount,gbifID,family,taxonRank,coordinateUncertaintyInMeters,year,basisOfRecord,institutionCode)

summary(datPinus)
summary(dat)
#remove records without coordinates

dat<-dat%>%
	filter(!is.na(decimalLongitude))%>%
	filter(!is.na(decimalLatitude))

datPinus<-datPinus%>%
	filter(!is.na(decimalLongitude))%>%
	filter(!is.na(decimalLatitude))


summary(dat)

#visualize the data on a map to get an overview

wm<-borders("world",colour="gray50",fill="gray50")
ggplot()+coord_fixed()+wm+
	geom_point(data=datPinus,aes(x=decimalLongitude,y=decimalLatitude),colour="darkred",size=0.5)+
theme_bw()

#convert country code from ISO2c to ISO3c
datPinus$countryCode<-countrycode(datPinus$countryCode,origin="iso2c",destination="iso3c")

#flag problems
#see coordinates, zero coordinates, coordinate-country mismatches, coordinates assignes to country and province centrooids, coordinates within city areas, outlier coordinates and coordinates asigned to biodiversity institutions
 install.packages("rnaturalearthdata")
library(rnaturalearthdata)
library(CoordinateCleaner)
dfPinus<-data.frame(datPinus)
flags<-clean_coordinates(x=dfPinus,lon="decimalLongitude",lat="decimalLatitude",countries="countryCode",species="species",tests=c("capitals","centroids","equal","gbif","sinstituions","zeros","countries"))
summary(flags)
plot(flags,lon="decimalLongitude",lat="decimalLatitude")

summary(dat)
summary(flags)

#eclude problematic records

dat_clean<-dfPinus[flags$.summary,]

dfPinus[flags$.summary,]



#flagged records

dfPinus[!flags$.summary,]

dat.flag<-dfPinus[!flags$.summary,]

##or using magrittr pipe, results directly in a data.frame comprising only cleaned records

names(dat)[2:3]<-c("decimallongitude","decimallatitude")


clean <-dfPinus>%
	cc_val()%>%
	cc_equ()%>%
	cc_cap()%>%
	cc_cen()%>%
	cc_coun(iso3="countryCode")%>%
	cc_gbif()%>%
	cc_inst()%>%
	cc_sea()%>%
	cc_zero()%>%
	cc_outl()%>%
	cc_dupl()

datPinus
clean

##rajouter les colonnes des données "flag" dans datPinus
	
datPinus%>%
	as_tibble()%>%
	mutate(val=cc_val(.,value="flagged"),equ=cc_equ(.,value="flagged"),cap=cc_cap(.,value="flagged"),cen=cc_cen(.,value="flagged"),gbif=cc_gbif(.,value="flagged"),inst=cc_inst(.,value="flagged"),sea=cc_sea(.,value="flagged"),zero=cc_zero(.,value="flagged"),outl=cc_outl(.,value="flagged"),dupl=cc_dupl(.,value="flagged"))
summary(datPinus)

#records with low coordinate precision

summary(dat_clean)

hist(dat_clean$coordinateUncertaintyInMeters / 1000, breaks=20)

#remove records with a precision below 10 km
#do I need to reduce that to 1 km?
dat_clean<-dat_clean%>%
	filter(coordinateUncertaintyInMeters/1000<=10|is.na(coordinateUncertaintyInMeters))

hist(dat_clean$coordinateUncertaintyInMeters /1000, breaks=20)

table(datPinus$basisOfRecord)
table(dat_clean$year)
summary(dat_clean)

#remove records from before the WW2


dat_clean<-dat_clean%>%
	filter(year>1945)

#if we want to clean the taxon names check the taxize R package
table(dat_clean$family)

##checking plant species range :
use species ranges from external sources as reference and flag all records falling outside these ranges
for plants check the botanical countries of the Wold Checklist of selected plant families
#we need to check :
#which species fall out of their ranges
#which species are from the pyrenees
#which species are from mountain range
#create a SpatialPolygonDataFrame and put it into the fncion cc_iucn
#uwe the GTS database to identify countries where the target species is native
#identify countries where the species is introduced and naturalized ; GLOBAL INVASIVE SPECIES DATAASE (GISD)
#GLOBAL REGISTER OF INTRODUCED AND INVASIVE SPECIES
#country spatial layer obtained from the global administrative database

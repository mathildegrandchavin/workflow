#### SPECIESS SORT AND COORDINATE CLEAN
#### AUTHOR

### 0 Set vars  ====
library (data.table)
library (tidyverse)
library(foreach)
library(dplyr)
library(raster)
library(countrycode)
library(rnaturalearthdata)
library(CoordinateCleaner)

setwd ("D:/quercus/Utilisateurs/mgrandchavin")

pasteSpName<-function (x) paste(strsplit (x,' ')[[1]],collapse="_")

### 1 Read in all occurrences and invasive species data ====
allOcc <- data.table::fread('data_occ/occurences/occurences.csv')
allOcc <- as_tibble(allOcc)

### a enlever :)
#allOcc = allOcc[1:100,]

summary(allOcc)
allOcc$taxonRank

gisd <- read.csv2("data_occ/gisd.csv")
summary(gisd)


### Read in environmental data
mesoTair_Tmax_list<-list.files(path="climate data/mesoTair_p0320_vOK/Tmax",all.files=TRUE)
mesoTair_Tmax_list = mesoTair_Tmax_list[-c(1,2,3)]
mesoTair_Tmax_list <- lapply(mesoTair_Tmax_list, function(x){paste("climate data/mesoTair_p0320_vOK/Tmax",x,sep="/")})
mesoTair_Tmax<-lapply(mesoTair_Tmax_list,raster)

mesoTair_Tmin_list<-list.files(path="climate data/mesoTair_p0320_vOK/Tmin",all.files=TRUE)
mesoTair_Tmin_list = mesoTair_Tmin_list[-c(1,2)]
mesoTair_Tmin_list <- lapply(mesoTair_Tmin_list, function(x){paste("climate data/mesoTair_p0320_vOK/Tmin",x,sep="/")})
mesoTair_Tmin<-lapply(mesoTair_Tmin_list,raster)
mesoTair_Tmin_list

microTair_Tmax_list<-list.files(path="climate data/µTair_p0320/Tmax_p0320",all.files=TRUE)
microTair_Tmax_list = microTair_Tmax_list[-c(1,2)]
microTair_Tmax_list <- lapply(microTair_Tmax_list, function(x){paste("climate data/µTair_p0320/Tmax_p0320",x,sep="/")})
microTair_Tmax<-lapply(microTair_Tmax_list,raster)

microTair_Tmin_list<-list.files(path="climate data/µTair_p0320/Tmin_p0320",all.files=TRUE)
microTair_Tmin_list = microTair_Tmin_list[-c(1,2)]
microTair_Tmin_list <- lapply(microTair_Tmin_list, function(x){paste("climate data/µTair_p0320/Tmin_p0320",x,sep="/")})
microTair_Tmin<-lapply(microTair_Tmin_list,raster)

envsmeso<-raster::stack(mesoTair_Tmax,mesoTair_Tmin)
envsmicro<-raster::stack(microTair_Tmax,microTair_Tmin)


### 2 Select taxa to species level
allOcc = allOcc %>% 
  filter (taxonRank =='SPECIES')

### 3 Create species list and species folders ====

#species list with non duplicated rows
allSp =  pull(allOcc,species) %>% unique

#create a directory for all species
pathOutput = paste0(getwd(),'/species')
dir.create (pathOutput)




#foreach (sp=allSp) %dopar% ###this tells how many times the loop should be iterated and makes the link with the species list
sp = "Galium verum"
#{
	#get the occurence data of species from the general df 
	sp_raw_occ = allOcc %>% filter (species == sp)

	#create a folder for each species
	sp_dir = paste0(pathOutput,'/',pasteSpName(sp))
	dir.create (sp_dir)
	readr::write_rds (x = sp_raw_occ,
                  file = paste0(sp_dir,'/occRaw_',pasteSpName(sp),'.rds') )
	###clean data
	dat_sp<-sp_raw_occ %>%
	#select only the columns we need for the cleaning process
	dplyr::select(species,decimalLongitude,decimalLatitude,countryCode,individualCount,gbifID,family,taxonRank,coordinateUncertaintyInMeters,year,basisOfRecord,institutionCode)	


	#step 1: remove records without coordinates
	dat_sp<-dat_sp %>%
		filter(!is.na(decimalLongitude))%>%
		filter(!is.na(decimalLatitude))

	#Convert country code from ISOc to ISO3c
	dat_sp$countryCode<-countrycode(dat_sp$countryCode,origin="iso2c",destination="iso3c")
	
	#step 2: flag problems
	#duplicated data
	
	dat_sp<-data.frame(dat_sp)
	flags<-CoordinateCleaner::clean_coordinates(dat_sp,lon="decimalLongitude",lat="decimalLatitude",countries="countryCode",tests=c("validity","duplicates","seas","institutions","centroids","outliers","capitals","urban","zeros","countries"))
	
	#data frame with problematic records
	dat_sp_flag<-dat_sp[!flags$.summary,]
	#data frame with clean records
	dat_sp_clean<-dat_sp[flags$.summary,]

	#is the species invasive?
	invasivematch<- dat_sp_clean$species %in% gisd$Species
	dat_sp_clean<-dat_sp_clean[! invasivematch ,]

	#keep only precise coordinates (precision >0,05 km)
	dat_sp_clean<-dat_sp_clean%>%
	filter(coordinateUncertaintyInMeters/1000<=0.05|is.na(coordinateUncertaintyInMeters))

	#write in the directory the clean data
	readr::write_rds (x = dat_sp_clean,
                  file = paste0(sp_dir,'/occClean_',pasteSpName(sp),'.rds') )

#}

#foreach (sp=allSp) %dopar%

#{	
	#read the rds file of clean data
	dat_sp_clean<-readr::read_rds(file = paste0(sp_dir,'/occClean_',pasteSpName(sp),'.rds'))


	#fit a SDM
	
	#make sure latitude and longitude are numeric
	dat_sp_clean$decimalLatitude<-as.numeric(dat_sp_clean$decimalLatitude)
	dat_sp_clean$decimalLongitude<-as.numeric(dat_sp_clean$decimalLongitude)
	
	#spatial thinning of 1 km
	
	output<-spThin::thin(dat_sp_clean,lat.col="decimalLatitude",long.col="decimalLongitude",spec.col="species",thin.par=1,reps=100,locs.thinned.list.return=TRUE,write.files=FALSE,verbose=FALSE)
	
	#maximise the number of localities we proceed with
	
	#find the iteration that returns the max number of occurrences
	maxThin<-which(sapply(output, nrow)==max(sapply(output,nrow)))
	
	#if there is more than one max, pick the first one
	
	maxThin<-output[[ifelse(length(maxThin)>1,maxThin[1],maxThin)]]
	
	#subset occs to match only thinned occs
	dat_sp_clean<-dat_sp_clean[as.numeric(rownames(maxThin)),]
	
	#put the species occurence in the same crs as the raster
	#first set the original crs
	dat_sp_clean_sf<-sf::st_as_sf(dat_sp_clean,coords=c("decimalLongitude","decimalLatitude"),crs=sf::st_crs(4326))


	#extract the crs of the raster
	crs_destination<-sf::st_crs(raster::crs(envsmeso))
	
	#transform the occurence data crs into the raster crs
	dat_sp_clean_sf<-sf::st_transform(dat_sp_clean_sf,crs_destination)

	#select the zone in which to take occurence data
	#our extent here is the pyrenees
	
	StudyAreaCoord<-rbind(c(3299875,2074100),c(3762200,2074100),c(3762200,2328325),c(3299875,2328325),c(3299875,2074100))
	
	p<-sf::st_sf(sf::st_sfc(sf::st_polygon(list(StudyAreaCoord))),crs=crs_destination)

	#we crop so that the environmental raster has the same extent as the study area
	envsmeso<-crop(envsmeso,StudyAreaCoord)
	occs.geom<-dat_sp_clean_sf[,"geometry"]
	
	
	#we intersect the occ points with the study area to select only those inside
	intersection<-sf::st_intersection(p,occs.geom)
	
	intersect.rowNums<-as.numeric(which(!(is.na(intersection))))
	occs.geom<-occs.geom[intersect.rowNums,]
	

	#SDM with MESO DATA
	#extract environmental values at occ grid cells
	
	loc.vals<-raster::extract(envsmeso,occs.geom[,"geometry"])

	#remove occs without environmental values
	dat_sp_clean_sf<-dat_sp_clean_sf[!is.na(loc.vals),]

	#selection of background #for now, we choose the same background extent as the study area
	#crop the environmental rasters by the background extent shape
	envsmesoCrop<-raster::crop(envsmeso,StudyAreaPoly)

	#mask environmental variables and take a random sample of background values 
	
	envsmesoMask<-raster::mask(envsmesoCrop,p)
	#envsmesoMask<-raster::mask(envsmesoCrop,StudyAreaPoly)

	bg.xy<-dismo::randomPoints(envsmesoMask,n=50000, tryf=10)

	#convert matrix output to data frame
	
	bg.xy<-as.data.frame(bg.xy)
	

	#partition occurrence data
	
	occs.xy<-as.data.frame((sf::st_coordinates(dat_sp_clean_sf)))

	
	group.data<-ENMeval::get.checkerboard2(occ=occs.xy,env=envsmesoMask,bg=bg.xy,aggregation.factor=2)

	#pull out the occurences and background partition group numbers from the list

	
	occs.grp<-group.data[[1]]

	bg.grp<-group.data[[2]]

	
	#build MaxEnt Model

	#define the vector of regularization mutipliers to test

	rms<-seq(1,2,1)


	#try something to fix the pb of variables of not the same length
	
	occs.grp <- rbind(occs.xy[1, ] , occs.grp)
	occs.grp <- occs.grp[-1,]
	
	bg.grp<-rbind(bg.xy[1,],bg.grp)
	bg.grp<-bg.grp[-1,]
	
	#iterate model building over all chosen paramter settings
	
	colnames(bg.xy) <- c("x", "y")
	colnames(occs.xy) <- c("x", "y")
	e<-ENMeval::ENMevaluate(occ=occs.xy,env=envsmesoMask,bg.coords=bg.xy,RMvalues=rms,fc="L",method="user",occ.grp=occs.grp,bg.grp=bg.grp,clamp=TRUE,algorithm="maxnet")
  f<-dismo::maxent(envsmeso, occs.xy,bg.xy)
  
  
	#unpack the results data frame the list of models, and the rasterstack of raw predictions
	evaltbl<-e@results
	evalMods<-e@models
	names(evalMods)<-e@tune.settings$tune.args
	evalPreds<-e@predictions

	#view ENMeval results

	ENMeval::evalplot.stats(e,stats="auc.val","rm","fc")

	#select model from the list
	mod<-evalMods[["rm.1_fc.L"]]

	#generate cloglog prediction
	pred<-predictMaxnet(mod,envsmesoMask,type="cloglog",clamp=TRUE)
	
	#write model and model prediction
	readr::write_rds (x = mod,
                  file = paste0(sp_dir,'/model_',pasteSpName(sp),'.rds') )

	readr::write_rds (x = pred,
                  file = paste0(sp_dir,'/modelPred',pasteSpName(sp),'.rds') )

#}

#plot
plot(pred)
mapview::mapview()

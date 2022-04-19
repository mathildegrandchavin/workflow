###########
# ECOSILVA
# MATHILDE GRAND-CHAVIN
# SDM WORKFLOW TEMPLATE

###1 -Occurence data
#entrée : presences/absences.csv
#sortie : presence_absence_background.csv




###2 - Predictors
#entrée: pred.tif
#sortie: presence_absence_backgroud_pred.csv

#VIF
model<-lm(xy~variable1+v2)
summary(model)
library(car)
vif(model)
#variable selection
#study area


###3 - Model fitting

#entrée: presence_absence_background_pred.csv
#sortie: presence_absence_background_k.csv
#sortie: model objects
###4 - Assessment

#entrée: model object.rds
#sortie: stats.csv



###5- Model selection

#entrée: stats.csv
#sortie: selected/ averaged model

###6- Thresholding

#entrée: averaged model
#sortie: raw output
#sortie : cloglog transf
#sortie : binary pred


###7- Mapping

#raw
#transformed
#binary

###8- Transfer

#entrée: newpred.tif

###9 - Reporting


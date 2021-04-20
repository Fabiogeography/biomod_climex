##################################################################
##### Ensemble biomod2 and CLIMEX SDM projections #################
##### Written by: Regan Early ####################################
##### Written on: 08/11/2018 #####################################
##### Modified on:  20/04/2021  ##################################
##################################################################

library(raster)
library(REAT) ## Calculate weighted coefficient of variance. If you have trouble installing REAT scroll to the bottom

##### Set up parameters #####
id <- "TA" ## A unique identifier for the model run (e.g. indicating species name, environmental variables used, etc...). Should be short. In  
eval.stat <- "ROC" ## Choice of evaluation statistic. Support is provided here for ROC, TSS, and KAPPA, though others are possible. Only AUC can be calculated without a choice of threshold. See ?BIOMOD_Modeling for all options that could be applied to biomod results and modEvAmethods("multModEv") for all options that could be applied to climex results using this code.
incl.climex <- T ## Include CLIMEX in ensemble even if its evaluation statistic falls below the threshold
thresh <- 0.8 ## The value above which models should be included in the ensemble.
wd <- "E:/NON_PROJECT/TUTA_ABSOLUTA/DISTRIBUTION_MODELS/R/TO_SHARE_GITHUB/"
wd.out <- paste0(wd,id,"/")

##### Load in rasters of biomod2 model projections and save the ones to be averaged as tifs #####

### Note the filepath given below is based on the folders automatically created by biomod2
### It will be correct if you used the code above to set the working directory: wd.out <- paste0(wd,id,"/")
rs.biomod <- stack(paste0(wd.out, id, "/proj_",id,"_Full/proj_",id, "_Full_",id,".gri"))
rs.biomod <- rs.biomod / 1000 ## convert to 0-1 numerical scale

### Identify the biomod2 models to be ensembled

## Obtain evaluation statistics from cross-validation for biomod models
eval <- read.csv(paste0(wd.out, eval.stat,"_eval_",id,".csv"))
rownames(eval) <- eval$model

## Select models whose evaluation statistic falls above the threshold
eval.biomod <- as.matrix(eval[(rownames(eval)!="climex"), grep("PA", colnames(eval))]) ## Make a table without the CLIMEX stats
eval.biomod[eval.biomod<thresh] <- NA ## Remove the names of models whose evaluation statistic falls below the chosen threshold

nms <- matrix(names(rs.biomod), ncol=ncol(eval.biomod), nrow=(length(names(rs.biomod))/ncol(eval.biomod)))
nms[is.na(eval.biomod)] <- NA
nms <- na.omit(as.vector(nms)) ## A vector of names of models whose evaluation statistic falls above the threshold

rs.biomod <- rs.biomod[[nms]] # Only the rasters that fall above the chosen threshold

##### Load in rasters of CLIMEX model projection/s and save the one/s to be averaged as tifs ###
climex <- raster(paste0(wd, "san_ei"))
climex <- climex/100 ## convert to 0-1 numeric scale
climex <- projectRaster(climex, crs=crs(rs.biomod)) ## transform projection to match biomod models
climex <- resample(climex, rs.biomod[[1]]) ## resample to resolution of biomod models
names(climex) <- "climex"

##### Calculate ensemble #####

### If CLIMEX should be included in the ensemble regardless of it's evaluation statistic, set it's evaluation statistic to be 1.
### Note this also sets the weight for CLIMEX to 1. 
if(isTRUE(incl.climex)) {
  eval.climex <- 1
} else {
  ### If CLIMEX should only be included in the ensemble if its mean evaluation statistic is above the threshold, obtain the mean evaluation statistic from the multiple validations of the full model for CLIMEX
  eval.climex <- eval[(rownames(eval)=="climex"), "mean"] 
}

### Stack the raster projections from the models to be ensembled.
### Note if a single CLIMEX model is supplied, this approach includes CLIMEX only once in the ensemble.
### If many biomod2 models have been created, e.g. if several pseudo-absence datasets were selected, this may mean CLIMEX has little influence on the ensemble.
### This could be altered, but the authors recommend that inclusion in an ensemble and weightings are based on clearly documented, quantitative evaluation statistics.
if(eval.climex > thresh) { ## If the eval stat for climex falls at or above the chosen threshold:
  rs <- stack(rs.biomod, climex) 
  wts <- na.omit(as.vector(c(eval.biomod, eval.climex))) ## calculate the weights for the ensembles.
} else {  ## If the eval stat for climex falls below the chosen threshold:
  rs <- rs.biomod
  wts <- na.omit(as.vector(eval.biomod)) ## calculate the weights for the ensembles.
}

### Calculate and save the ensemble (the mean value of the projections from all selected models, weighted by the evaluation statistics of those models)
wmean <- weighted.mean(rs, wts) 
writeRaster(wmean, paste0(wd.out, "/wmean_", eval.stat, "_", id,".tiff"))

### Calculate and save the uncertainty associated with the ensemble (the coefficient of variance of the projections from all selected models, weighted by the evaluation statistics of those models)
wcv <- calc(rs, fun=REAT::cv, weighting=wts, is.sample=F, wmean=T) ## see documentation of REAT::cv for how weights are calculated
writeRaster(wcv, paste0(wd.out, "/wcv_", eval.stat, "_", id, ".tiff"))

##### If you have trouble installing REAT because the package version does not match your version of R try this: #####
library(rlang, lib.loc=...)
library(devtools, lib.loc=...) 

install_version("REAT", version = "3.0.2",
                repos = "http://cran.us.r-project.org",
                lib=...)

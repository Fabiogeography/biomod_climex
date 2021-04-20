##################################################################
##### Evaluate biomod and CLIMEX SDM projections ######
##### Written by: Regan Early ####################################
##### Written on: 08/11/2018 #####################################
##### Modified on:  20/04/2021  ##################################
##################################################################

library(raster) ## To load and manipulate rasters
library(biomod2) ## To run and evaluate biomod species distribution models
library(dplyr) ## For function written get_PAtab
library(modEvA) ## Calculate evaluation satistics

wd <- "E:/NON_PROJECT/TUTA_ABSOLUTA/DISTRIBUTION_MODELS/R/TO_SHARE_GITHUB/"

##### Set up run #####
nreps <- 4 ## Number of pseudo-absence datasets created and used
techs <- c('GLM','GAM', 'SRE','RF','ANN','FDA','MAXENT.Phillips','GBM','CTA','MARS') ## biomod techniques to use
id <- "TA" ## A unique identifier for the model run (e.g. indicating species name, environmental variables used, etc...). Should be short. In  
eval.stat <- "KAPPA" ## Choice of evaluation statistic. Support is provided here for ROC, TSS, and KAPPA, though others are possible. Only AUC can be calculated without a choice of threshold. See ?BIOMOD_Modeling for all options that could be applied to biomod results and modEvAmethods("multModEv") for all options that could be applied to climex results using this code.
generate.SDM <- T ## Whether this code will make biomod predictions (T) or load predictions made previously (F)
generate.PA <- F ## Whether biomod should generate and save its own pseudo-absence sets (T) or if the user will supply their own (F).

### Prepare output directory
wd.out <- paste0(wd,id,"/")
dir.create(wd.out)
setwd(wd.out)

### Function for extracting presence-absence table. By D Georges, https://rpubs.com/dgeorges/416446
get_PAtab <- function(bfd){
  dplyr::bind_cols(
    x = bfd@coord[, 1],
    y = bfd@coord[, 2],
    status = bfd@data.species,
    bfd@PA
  )}

### Load environmental data
env.bkg <- stack(paste0(wd, "env_bkg")) ### Data covering the region where the species is recorded present/absent and which is used to construct the SDM
env.p <- stack(paste0(wd, "env_p")) ### Data covering the region where the species' range will be projected

### If the Maxent.Phillips method of constructing an SDM is to be used, give the pathway to the folder where the software is saved
myBiomodOption <- BIOMOD_ModelingOptions(
  MAXENT.Phillips = list( path_to_maxent.jar = "E:/NON_PROJECT/TUTA_ABSOLUTA/DISTRIBUTION_MODELS/") ## Ensure there are no spaces or non-alphanumeric characters in this pathway.
)

if (isTRUE(generate.SDM)) {
  ##### Prepare data for biomod SDMs #####
  if (isTRUE(generate.PA)) { ## biomod places the pseudo-absences
    
    ### Load species presence data
    pres.sp <- read.csv(paste0(wd, "/pres.csv")) ## A data frame with columns 'sp' (containing 1, indicating the species is present), 'x' and 'y' (the coordinates of the presences)
    coordinates(pres.sp) <- ~x+y ## Convert the data frame to a points spatial data frame
    
    ### Enter presence data and use biomod to place pseudo-absences 
    ### Here pseudo-absences are placed randomly throughout the background. The number of pseudo-absences is made to be the same as the number of presences. 
    myBiomodData <- BIOMOD_FormatingData(resp.var = pres.sp, ## The response variable (i.e. presences)
                                         expl.var = env.bkg, ## The explanatory variables (i.e. background environmental data)
                                         resp.name = id,  ## The unique identifier for the model created
                                         PA.nb.rep = nreps, ## The number of pseudo-absence datasets created and used
                                         PA.nb.absences = nrow(pres.sp), ## Number of pseudo-absences
                                         PA.strategy = 'random', ## How to place the pseudo-absences
                                         na.rm = TRUE) ## Remove any presence / pseudo-absence points with missing environmental data from the analysis
    
    ### Extract and save the pseudo-absence data created above. A tibble with columns:
    ### 'x' and 'y' (coordinates), 
    ### 'status' (containing 1 or NA, indicating the species is present or absent respectively),
    ### 'PA1', 'PA2' ... 'PAn' where n is nreps, the number of pseudo-absence datasets created. Each column contains 'TRUE' or 'FALSE', indicating whether the corresponding point location is used in the named presence / pseudo-absence dataset. 
    pa <- get_PAtab(myBiomodData)
    write.csv(pa, paste0(wd.out,"/pseudoabs.csv"), row.names=F)
    
  } else { ## User uploads their own table of presences / pseudo-absences. A tibble with columns:
    ### 'x' and 'y' (coordinates), 
    ### 'status' (containing 1 or NA, indicating the species is present or absent respectively),
    ### 'PA1', 'PA2' ... 'PAn' where n is nreps, the number of pseudo-absence datasets created. Each column contains 'TRUE' or 'FALSE', indicating whether the corresponding point location is used in the named presence / pseudo-absence dataset. 
    pa.tab <- read.csv(paste0(wd,"/presence_pseudoabs.csv"))
    
    resp.var <- pa.tab[, c("x","y","status")] ## The presences
    coordinates(resp.var) <- ~x+y ## Convert the data frame to a points spatial data frame
    
    pa <- pa.tab[,grep("PA",colnames(pa.tab))] ## The PA.table indicating the pseudo-absences selection
    
    ### Create biomod data with the manual data
    myBiomodData <- BIOMOD_FormatingData(resp.var = resp.var, ## The response variable (i.e. presences)
                                         expl.var = env.bkg, ## The explanatory variables (i.e. background environmental data)
                                         resp.name = id, ## The unique identifier for the model created
                                         PA.nb.rep = nreps, ## The number of pseudo-absence datasets created and used
                                         PA.strategy = 'user.defined',
                                         PA.table = pa, ## Defined above.
                                         na.rm = TRUE) ## Remove any presence / pseudo-absence points with missing environmental data from the analysis
  }
  
  ##### Calculate and project the biomod SDMs #####
  myBiomodModelOut <- BIOMOD_Modeling(
    myBiomodData, ## Defined above.
    models = techs, ## SDM techniques. Defined above. 
    models.options = myBiomodOption, ## Defined above.
    NbRunEval=1, ## 1 cross-validation runs per pseudo-absence set. 
    DataSplit=70, ## Proportion of data randomly assigned to the calibration data. Models made with calibration data have '_RUN1' added to name (or RUN2 etc... if NBRunEval >1). Others have '_Full'
    VarImport=3, ## Number of permutations with which to evaluate variable importance
    models.eval.meth = eval.stat,
    SaveObj = TRUE, ## save data to hard drive
    rescal.all.models = F, ## Scales all model predictions from 0-1 with binomial GLM
    do.full.models = TRUE, ## In addition to models made with the calibration data for cross-validation. Have 'Full' added to name
    modeling.id = id)
  
  ##### Project the full models made with each presence and pseudo-absence dataset #####
  ### Name the models to be projected.
  ### Can project either models made with all presence / pseudo-absence data in a dataset ('_Full')
  ### or with the calibration data only ('_RUN1', '_RUN2' ... '_RUNn' where n is NBRunEval) 
  
  mods2proj <- vector()
  for(a in paste0("PA",c(1:nreps))) { mods2proj <- c(mods2proj, paste0(id,"_",a,"_Full_", techs)) }
  
  ### Project models
  myBiomodProj <- BIOMOD_Projection(
    modeling.output = myBiomodModelOut, ## The object containing the SDMs. Defined above.
    new.env = env.p, ## Environmental data covering the region where the species' range will be projected
    proj.name = paste0(id, '_Full'), ## Can project either models made with all presence / pseudo-absence data in a dataset ('_FULL') or with the calibration data only ('_RUN1', '_RUN2' etc...) 
    selected.models = mods2proj, ## Defined above.
    binary.meth = NULL, ## Whether the continuous projections are thresholded into suitably/unsuitable
    filtered.meth = NULL, ## Thresholding technique
    compress = 'xz', ## Compression format of objects stored on hard drive
    clamping.mask = F, ## A mask would identify locations where predictions are uncertain because the values of the environmental variables are outside the range used for calibrating the models
    output.format = '.grd' ## Format in which to save the outputted rasters
  )
  
  ##### Tabulate evaluation statistics for individual biomod models (from validation data #####
  myBiomodModelEval <- get_evaluations(myBiomodModelOut)
  eval <- as.data.frame(myBiomodModelEval[eval.stat,"Testing.data",,"RUN1",]) ## Get the cross-validation evaluation statistic for all models and pseudo-absence sets. Need to change "RUN1" if do multiple cross-evaluations.
  
} else {
  
  ### Load the evaluation statistics saved for all the models made in a previous biomod2 run.
  ### The value for models.eval.meth supplied in BIOMOD_modeling must be the same as the value given for eval.stat, above.
  ### Note the filepath given below is based on the folders automatically created by biomod2
  ### It will be correct if you used the code above to set the working directory: wd.out <- paste0(wd,id,"/")
  load(paste0(wd.out,id,"/.BIOMOD_DATA/",id,"/models.evaluation"))
  eval <- as.data.frame(models.evaluation[eval.stat,"Testing.data",,"RUN1",]) ## Get the cross-validation evaluation statistic for all models and pseudo-absence sets. Need to change "RUN1" if do multiple cross-evaluations.
  
  ## User uploads their own table of presences / pseudo-absences for climex evaluation. A tibble with columns:
  ### 'x' and 'y' (coordinates), 
  ### 'status' (containing 1 or NA, indicating the species is present or absent respectively),
  ### 'PA1', 'PA2' ... 'PAn' where n is nreps, the number of pseudo-absence datasets created. Each column contains 'TRUE' or 'FALSE', indicating whether the corresponding point location is used in the named presence / pseudo-absence dataset. 
  pa.tab <- read.csv(paste0(wd,"/presence_pseudoabs.csv"))
}

##### Load and evaluate CLIMEX model #####
### Obtain the correct name for the chosen evaluation statistic - as used by multModEv
if(eval.stat=="ROC") {eval.stat.climex <- "AUC"} 
if(eval.stat=="KAPPA") {eval.stat.climex <- "kappa"} 
if(eval.stat=="TSS") {eval.stat.climex <- eval.stat}

### Load CLIMEX model
climex <- raster(paste0(wd, "san_ei"))
climex <- climex/100 ## convert to 0-1 numeric scale
climex <- projectRaster(climex, crs=crs(env.p)) ## transform projection to match biomod models
climex <- resample(climex, env.p[[1]]) ## resample to resolution of biomod models

### Evaluate the CLIMEX model 
evaluate.climex <- function (x) {
  ### Obtain vectors of the predictions (pred) and the presence / pseudo-absences (occ)
  pa.rep <- as.data.frame(pa.tab[,c("x","y","status", paste0("PA",x))]) ## Select the appropriate pseudo-absence column
  colnames(pa.rep) <- c("x","y","status", "PA")
  pa.rep <- as.data.frame(pa.rep[pa.rep$PA==T,]) ## Select only the pseudo-absences used in the relevent pseudo-absence set
  occ <- pa.rep[pa.rep$PA==T,"status"] ## Create a vector of the presence / pseudo-absences
  occ[is.na(occ)] <- 0 ## Give pseudo-absences a value of 0.
  
  pred <- extract(climex, pa.rep[,c("x","y")]) ## Extract the climex predictions at the PA locations

  ### Obtain evaluation statistic
  out <- multModEv(obs.data = as.data.frame(occ), 
                   pred.data = as.data.frame(pred), 
                   measures = eval.stat.climex) 
                   # thresh = "preval") ## Note the threshold must be set appropriately if kappa or TSS are required.
  
  out <- out[,-1] ## Retain only the evaluation statistic
  
}

climex.eval <- lapply(1:nreps, evaluate.climex)

##### Unite evaluation statistics for biomod2 and CLIMEX models
eval["climex", ] <- climex.eval
eval$mean  <- apply(eval, MARGIN=1, FUN="mean")
eval$sd <- apply(eval[,names(eval) !="mean"], MARGIN=1, FUN="sd")
eval$model <- rownames(eval)
write.csv(eval, paste0(wd.out, eval.stat,"_eval_",id,".csv"), row.names = F)



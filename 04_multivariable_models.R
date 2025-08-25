################################################################################
# Ecological analysis of the association between tuberculosis incidence and
# ZIP code-level socioeconomic variables, 2008--2019
################################################################################

# Univariable models associating TB notification rates with socioeconomic variables.
# Spatial random effects are modeled with a BYM2 distribution.
# Models use the inla function on the pre-treated dataset obtained with the R code in file 01_data_pretreatment.R.

library(tidyverse)
library(magrittr)
library(sf)
library(Cairo)
library(INLA)
library(ggpubr)


#--------- Read the useful function bundle ------------
source( "include/functions.R", encoding = "UTF-8" )

#----- Read the pre-treated dataset ----------------
map_PMSI21_2 <- readRDS("pretreated_dataset_sf.rds" %>% respath())
n_areas <- nrow(map_PMSI21_2) / 2

#---------- Read the neighboring matrix of PMSI21 codes ---------
# Neighbour list object
tb.adj <- readRDS("neighboring_structure_PMSI21.rds" %>% datapath)

# Create the adjacency matrix
W <- spdep::nb2mat( neighbours = tb.adj, style = "B", zero.policy = TRUE )

#------ Read shapefiles of French "départements" (districts) and regions -------
map_dep <- readRDS(datapath("France_shapefiles/map_PMSI21_dep.rds"))  
map_reg <- readRDS(datapath("France_shapefiles/map_PMSI21_reg.rds"))


# Priors for the BYM2 model
prior.prec.bym2 <- list(
  prec = list(
    prior = "pc.prec",
    param = c(0.5 / 0.31, 0.01)
  ),
  phi = list(
    prior = "pc",
    param = c(0.5, 1 / 3) 
  )
)


#------------------ Multivariable model ----------------------------------
listvar <- c("unemploy_log_g", "crowded_house_log_g", "house_income_log_g") 

form0 <-'
  n ~ 
   f(ID_PMSI21, 
     model="bym2",
     graph = W,
     scale.model = TRUE,
     constr = TRUE,
     hyper = prior.prec.bym2
  ) + dens_2010 +
'

form1 <- 'f({listvar}, 
     model = "rw2", 
     scale.model = TRUE, 
     constr = TRUE
  ) 
  ' %>% str_glue %>% paste( collapse=" + ") 

form <-  paste( form0, form1 ) %>% formula
form

model_multiv_bym2 <- inla( 
  formula = form,
  family = "poisson",
  data = map_PMSI21_2,
  E = Ei_knn,
  control.predictor = list( compute=TRUE ),
  control.compute = list( config=TRUE, waic=TRUE, dic=TRUE, cpo=TRUE, openmp.strategy="huge", return.marginals.predictor = TRUE ),
  verbose=T
)


model_multiv_bym2$summary.hyperpar

#saveRDS(model_multiv_bym2, "model_multiv_bym2.rds" %>% respath)
#model_multiv_bym2 <- readRDS( "model_multiv_bym2.rds" %>% respath)


#-------Sensitivity analyses -------------------------- 
# Using the expected values obtained after imputation with the "proportional" method (see article for details)
model_multiv_bym2_prop <- inla( 
  formula = form,
  family = "poisson",
  data = map_PMSI21_2,
  E = Ei_prop,
  control.predictor = list(compute=TRUE),
  control.compute = list(config=TRUE, waic=TRUE, dic=TRUE, cpo=TRUE, openmp.strategy="huge", return.marginals.predictor = TRUE),
  verbose=T
)
#saveRDS(model_multiv_bym2_prop, "model_multiv_bym2_prop.rds" %>% respath)
#model_multiv_bym2_prop <- readRDS( "model_multiv_bym2_prop.rds" %>% respath)


# Using the expected values obtained after imputation with the random forest method (see article for details)
model_multiv_bym2_RF <- inla( 
  formula = form,
  family = "poisson",
  data = map_PMSI21_2,
  E = Ei_RF,
  control.predictor = list(compute=TRUE),
  control.compute = list(config=TRUE, waic=TRUE, dic=TRUE, cpo=TRUE, openmp.strategy="huge", return.marginals.predictor = TRUE),
  verbose=T
)
#saveRDS(model_multiv_bym2_RF, "model_multiv_bym2_RF.rds" %>% respath)
#model_multiv_bym2_RF <- readRDS("model_multiv_bym2_RF.rds" %>% respath)

model_bym2_chom_suroc_rev_dens$waic$waic
model_bym2_chom_suroc_rev_dens_ventil$waic$waic
model_bym2_chom_suroc_rev_dens_RF$waic$waic

model_bym2_chom_suroc_rev_dens$dic$dic
model_bym2_chom_suroc_rev_dens_ventil$dic$dic
model_bym2_chom_suroc_rev_dens_RF$dic$dic

model_bym2_chom_suroc_rev_dens$summary.fixed
model_bym2_chom_suroc_rev_dens_ventil$summary.fixed
model_bym2_chom_suroc_rev_dens_RF$summary.fixed

model_bym2_chom_suroc_rev_dens$summary.hyper
model_bym2_chom_suroc_rev_dens_ventil$summary.hyper
model_bym2_chom_suroc_rev_dens_RF$summary.hyper
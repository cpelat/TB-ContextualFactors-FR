################################################################################
# Ecological analysis of the association between tuberculosis incidence and
# ZIP code-level socioeconomic variables, 2008--2019
################################################################################

# Spatial models with INLA using the pre-treated dataset obtained with the R
# code in file 01_data_pretreatment.R.
# Central analysis : 
# BYM2 distribution for spatial random effects
# Imputation with the kNN algorithm for missing values in case characteristics 
# (to compute the expected case counts)

library(tidyverse)
library(magrittr)
library(sf)
library(Cairo)
library(INLA)
library(ggpubr)
library(flextable)
#library(RColorBrewer)

set_flextable_defaults(  font.size = 10 )
setEPS() # Option to generate EPS file with the postscript function



#--------- Read the useful function bundle ------------
source( "include/include_functions.R", encoding = "UTF-8" )

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



#----------------- Modeling with INLA -----------------------
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

# "Empty" model: solely a spatial BYM2 random effect at the PMSI21 (ZIP) code level 
form_bym2_empty <- n ~ f(ID_PMSI21, 
                        model = "bym2",
                        graph = W,
                        scale.model = TRUE,
                        constr = TRUE,
                        hyper = prior.prec.bym2
) 

model_bym2_empty <- inla(form_bym2_empty,
                         family="poisson",
                         data = map_PMSI21_2,
                         E = Ei_knn,
                         control.compute = list( config=TRUE, 
                                                 waic=TRUE, dic=TRUE, 
                                                 cpo=TRUE, 
                                                 openmp.strategy="huge"
                         ),
                         verbose=TRUE
)


#saveRDS(model_bym2_empty, file="model_bym2_empty.rds" %>% respath)
#model_bym2_empty <- readRDS("model_bym2_empty.rds" %>% respath)

model_bym2_empty$summary.fixed
model_bym2_empty$summary.hyperpar
var_bym2_empty <- model_bym2_empty %>% get_marginal_var_sample() %>% get_marginal_mean()
var_bym2_empty

model_bym2_empty$dic$dic
model_bym2_empty$waic$waic




model_bym2_empty$summary.fixed
model_bym2_empty$summary.hyperpar
model_bym2_empty$dic$dic
model_bym2_empty$waic$waic




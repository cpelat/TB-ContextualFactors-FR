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
library(flextable)


set_flextable_defaults(font.size = 10)
setEPS() # Option to generate EPS file with the postscript function


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



#---------------- Univariable models ------------------
# NB: All continuous variables were log-transformed -except the FDEP-
# then grouped in 25 classes for INLA to be able to run RW2 models
list_var_cont_all <- c( 
  "crowded_house", 
  "unemploy",
  "house_income", 
  "manual_workers",
  "high_school_grad"
) %>% paste0("_log") %>%
  c("fdep") %>%
  paste0("_g")

# Create en empty list to stock the models
list_bym2_univ_rw2 <- list()

# Start with the continuous variables
for( v in list_var_cont_all){  
  print( v )
  list_bym2_univ_rw2[[v]] <- get_univ_model( 
    mydata=map_PMSI21_2,
    Ei="Ei_knn",
    form0 = 'n ~ 
          f(ID_PMSI21, 
           model="bym2",
           graph = W,
           scale.model = TRUE,
           constr = TRUE,
           hyper = prior.prec.bym2
           )',
    form_v = 'f({v}, 
              model = "rw2", 
              scale.model = TRUE, 
              constr = TRUE
              )')
}

# Add the population density level (categorial variable)
v <- "dens_2010"
list_bym2_univ_rw2[[v]] <- get_univ_model( 
  mydata=map_PMSI21_2,
  Ei="Ei_knn",
  form0 = 'n ~ 
                                  f( ID_PMSI21, 
                                     model="bym2",
                                     graph = W,
                                     scale.model = TRUE,
                                     constr = TRUE,
                                     hyper = prior.prec.bym2
                                  )', 
  form_v = '{v}')






#------------- Univariable models with period * variable interaction ----------
map_PMSI21_2 %<>%
  mutate(period_index = as.integer(period))

list_bym2_univ_rw2_inter <- list()

for(v in list_var_cont_all){
  print( v )
  list_bym2_univ_rw2_inter[[v]] <- get_univ_model( 
    mydata=map_PMSI21_2,
    Ei="Ei_knn",
    form0 = 'n ~ 
             f(ID_PMSI21, 
             model="bym2",
             graph = W,
             scale.model = TRUE,
             constr = TRUE,
             hyper = prior.prec.bym2
             )',
    form_v = 'f( {v}, 
              model = "rw2", 
              scale.model = TRUE, 
              constr = TRUE,
              replicate = period_index
              )'
  )
}

v <- "dens_2010"
list_bym2_univ_rw2_inter[[v]] <- get_univ_model( 
  mydata=map_PMSI21_2,
  Ei="Ei_knn",
  form0 = 'n ~ 
                                  f( ID_PMSI21, 
                                     model="bym2",
                                     graph = W,
                                     scale.model = TRUE,
                                     constr = TRUE,
                                     hyper = prior.prec.bym2
                                  )',
  form_v ="-1 + period:dens_2010" )  


#------------------------------------------------------------------------------#
#--- Show model results :
# Fixed effects
model_bym2_empty$summary.fixed

# Hyperparameters
model_bym2_empty$summary.hyperpar

# Variance of the spatial random effect
var_bym2_empty <- model_bym2_empty %>% get_marginal_var_sample() %>% get_marginal_mean()
var_bym2_empty

# Goddness-of-fit metrics: DIC and WAIC 
model_bym2_empty$dic$dic
model_bym2_empty$waic$waic






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
library(cowplot)
library(ggpubr)
library(flextable)

set_flextable_defaults(  font.size = 10 )
setEPS()  



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


#------ Read shapefile of PMSI21 codes (ZIP codes) for metropolitan France -----
map_PMSI21 <- readRDS( datapath( "France_shapefiles/map_PMSI21_simplified.rds" ) )  %>%
  arrange(PMSI21_CODE)

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

#-------------------- "Empty" model -----------------------------
# Contains solely a spatial BYM2 random effect at the PMSI21 (ZIP) code level 
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
#model_bym2_empty <- readRDS( "model_bym2_empty.rds" %>% respath)






#---------------- Sensitivity analysis: empty Leroux model -------
# see https://becarioprecario.bitbucket.io/inla-gitbook/ch-spatial.html
Q_diag <- Diagonal(x = sapply(tb.adj, length))
Q <- Q_diag - W 
C <- Diagonal(x = 1, n = nrow(map_PMSI21_2)/2) - Q

form_bym2_vide_leroux <- n ~ f(  ID_PMSI21, 
                                 model = "generic1",
                                 Cmatrix = C)  

model_leroux_empty <- inla(form_bym2_vide_leroux,
                           family="poisson",
                           data = map_PMSI21_2,
                           E = Ei_knn,
                           control.compute = list(dic = TRUE, 
                                                  waic = TRUE, 
                                                  cpo = TRUE,
                                                  config=TRUE, 
                                                  openmp.strategy="huge"),
                           control.predictor = list(compute = TRUE),
                           verbose=TRUE
)

#saveRDS(model_leroux_empty, file="model_leroux_empty.rds" %>% respath)
#model_leroux_empty <- readRDS( "model_leroux_empty.rds" %>% respath)






#------------------------------------------------------------------------------#
#--- Show model results :
# Fixed effects
model_bym2_empty$summary.fixed
model_leroux_empty$summary.fixed

# Hyperparameters
model_bym2_empty$summary.hyperpar
model_leroux_empty$summary.hyperpar

# Variance of the spatial random effect
var_bym2_empty <- model_bym2_empty %>% get_marginal_var_sample() %>% get_marginal_mean()
var_bym2_empty

var_leroux_empty <- model_leroux_empty %>% get_marginal_var_sample() %>% get_marginal_mean()
var_leroux_empty


# Goddness-of-fit metrics: DIC and WAIC 
model_bym2_empty$dic$dic
model_leroux_empty$dic$dic

model_bym2_empty$waic$waic
model_leroux_empty$waic$waic







#---------- Maps from the empty models --------------
# - map of the standardized notification rates (SNR: the model fitted values). 
# - map of the spatial random effects

# First create the dataset of the SNR and the spatial random effects
map_PMSI21_3 <- map_PMSI21 %>%
  bind_cols(  model_bym2_empty$summary.fitted.values[1:n_areas,]  %>% 
                dplyr::select( SNR_bym2=mean, SNR_bym2_low=`0.025quant`, SNR_bym2_up=`0.975quant` ),
              model_bym2_empty$summary.random[["ID_PMSI21"]][1:n_areas,] %>%
                dplyr::select( Ui_empty_bym2=mean, Ui_empty_bym2_low=`0.025quant`, Ui_empty_bym2_up=`0.975quant` ),
              
              model_leroux_empty$summary.fitted.values[ 1:n_areas, ]  %>% 
                dplyr::select( SNR_leroux=mean, SNR_leroux_low=`0.025quant`, SNR_leroux_up=`0.975quant` ),
              model_leroux_empty$summary.random[["ID_PMSI21"]][1:n_areas,] %>%
                dplyr::select( Ui_empty_leroux=mean, Ui_empty_leroux_low=`0.025quant`, Ui_empty_leroux_up=`0.975quant` )
  ) 


#----------- Create the SNR maps ----------
breakvalues <- c( 0.5, 1, 1.5 )
labels <- breakvalues
labels[ 1 ] <- paste0( "<", labels[ 1 ]  )
labels[ length(labels) ] <- paste0( ">", labels[ length(labels) ]  )
colorvalues <- rev(RColorBrewer::brewer.pal( n = 9, name = "RdYlGn"))  

p1_bym2 <- plotmap_inset( var= "SNR_bym2", 
                          mytitle="SNR, BYM2 model", 
                          mymap=map_PMSI21_3,
                          mymap_reg = map_reg,
                          mymap_dep= map_dep, 
                          mycolors=colorvalues, 
                          mybreaks=breakvalues, 
                          mylabels=labels, 
                          legend.title = "Standardized notification rate",
                          legend.position = "bottom",  
                          col.dep="black",
                          legend.text.size = 12
)


p1_leroux <- plotmap_inset( var= "SNR_leroux", 
                            mytitle="SNR, Leroux model", 
                            mymap=map_PMSI21_3,
                            mymap_reg = map_reg,
                            mymap_dep= map_dep, 
                            mycolors=colorvalues, 
                            mybreaks=breakvalues, 
                            mylabels=labels, 
                            legend.title = "Standardized notification rate",
                            legend.position = "bottom",  
                            col.dep="black",
                            legend.text.size = 12
)


#----------- Create the spatial random effects (Ui) maps ----------
breakvalues <- c( -1, -0.5, 0, 0.5, 1 )
labels <- breakvalues
labels[ 1 ] <- paste0( "<", labels[ 1 ]  )
labels[ length(labels) ] <- paste0( ">", labels[ length(labels) ]  )
colorvalues <- rev(RColorBrewer::brewer.pal( n = 9, name = "RdYlGn"))  

p2_bym2 <- plotmap_inset( var= "Ui_empty_bym2", 
                          mytitle="Spatial random effect, BYM2 model",  
                          mymap=map_PMSI21_3,
                          mymap_dep=map_dep, 
                          mymap_reg = map_reg,
                          mycolors=colorvalues, 
                          mybreaks=breakvalues, 
                          mylabels=labels, 
                          legend.title = "Spatial random effect",
                          legend.position = "bottom", 
                          col.dep="black",
                          legend.text.size = 12
)

p2_leroux <- plotmap_inset( var= "Ui_empty_leroux", 
                            mytitle="Spatial random effect, Leroux model", 
                            mymap=map_PMSI21_3,
                            mymap_dep=map_dep, 
                            mymap_reg = map_reg,
                            mycolors=colorvalues, 
                            mybreaks=breakvalues, 
                            mylabels=labels, 
                            legend.title = "Spatial random effect",
                            legend.position = "bottom", 
                            col.dep="black",
                            legend.text.size = 12
)





#---- Plot the SNR and Ui maps for the BYM2 model (central analysis) ------
ptot <- ggarrange(p1_bym2, p2_bym2, nrow=1, ncol=2)  

png( "SNR_Ui_empty_model_bym2.png" %>% respath, width=1000, height=600)
ptot
dev.off()

#----- Plot the SNR and Ui maps for the BYM2 and Leroux model ----
ptot_sens <- ggarrange(p1_bym2, p1_leroux, 
                       p2_bym2, p2_leroux, common.legend = TRUE, nrow=2, ncol=2)  
ptot_sens

png( "SNR_Ui_empty_model_bym2_leroux.png" %>% respath, width=1000, height=1000)
  ptot_sens
dev.off()



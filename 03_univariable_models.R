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

set_flextable_defaults(  font.size = 10 )

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


# Save the model list
saveRDS(list_bym2_univ_rw2, respath("list_bym2_univ_rw2.rds"))



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

# Save the model list
saveRDS(list_bym2_univ_rw2_inter, respath("list_bym2_univ_rw2_inter.rds"))





#-------------- Table of univariable models performance -----------------------
model_bym2_empty <- readRDS( "model_bym2_empty.rds" %>% respath)

# 1) DIC and WAIC
list_bym2_univ_rw2_0 <- c(
  list_bym2_univ_rw2[names(list_bym2_univ_rw2) %in% c(list_var_cont_all, "dens_2010")], 
  list_bym2_univ_rw2_inter %>% add_suffix_to_name)
list_bym2_univ_rw2_0[["Empty model"]] <- model_bym2_empty
res_rw <- get.DIC.table( model.list = list_bym2_univ_rw2_0, digits=2) 
res_rw

saveRDS(res_rw, "res_rw.rds" %>% respath)


# 2) Inter decile standardized rate ratio IdRR (for the models without period interaction)  
# 2.1) Find the decile values
q1 <- 0.1
q2 <- 0.9

list_q <- list()
for(v in list_var_cont_all){
  print(v)
  list_q[[v]] <- map_PMSI21_2[[v]] %>% quantile(probs = c(q1, q2)) 
}

# For the density, the deciles are the first and last levels :
map_PMSI21_2[["dens_2010"]] %>% table %>% prop.table
# sparse intermediate        dense 
#0.6349837    0.2641850    0.1008312 

list_q[["dens_2010"]] <- map_PMSI21_2[["dens_2010"]] %>% levels %>% extract(c(1,3))

list_q

# IdRR computation  
res_univ <- NULL
for(v in c(list_var_cont_all, "dens_2010")){
  print(v)
  # Generate a sample of size 100 from the posterior distribution
  samples <- inla.posterior.sample(100, list_bym2_univ_rw2[[v]], verbose=TRUE )
  
  # Calculate the IdRR
  res_univ <- bind_rows(res_univ,  get_IqRR_stats(v, samples, list_q) )
}

res_univ %<>%
  mutate(
    IdRR = str_glue("{formatC(IqRR_mean, format='f', digits=2, flag='0')} ({IqRR_0.025 %>% formatC(format='f', digits=2, flag='0')}, {IqRR_0.975 %>% formatC(format='f', digits=2, flag='0')})")) 

res_univ


#----- stopped here
# 3) Merger + ajouter le nom des variables
res_rw2 <- res_rw %>%
  filter(v0 %in% c(res_univ$variable, "Empty model")) %>% 
  left_join(res_univ %>% select(variable, IdRR, IqRR_mean), by=c(v0 = "variable")) %>% 
  mutate(Variable = v0 %>%
           str_remove_all('_log') %>% 
           str_remove_all('_g') %>% 
           get_var_name(nom_type="nom") %>%
           str_replace_all( fixed("\\n"), " " ),
         DIC = formatC(DIC, big.mark = ",", format="d")
  ) %>%
  arrange(DIC) %>%
  select(v0, 
         "Univariable model" = Variable, 
         DIC,
         "Inter-decile standardized rate ratio \n(95% Credible Interval)" = IqRR
  )  
res_rw2

saveRDS(res_rw2, "res_rw2.rds" %>% respath)
#res_rw2 <- readRDS("res_rw2.rds" %>% respath)

tt <- qflextable(res_rw2 %>% select(-v0)) 
tt
save_as_docx( tt, 
              values = NULL, 
              path="IqRR_DIC_model_univ_90.docx"  %>% respath, 
              pr_section = prop_section(page_size = page_size(orient = "portrait")), 
              align = "center")



#------------ Graphique des contributions des variables dans les modeles univariés sans slope aléatoire -----
listvar <- res_rw2$v0[!(res_rw2$v0 %in% c("dens_2010", "Empty model"))] 

titre <- get_var_name( listvar  %>% 
                         str_remove_all( "_g" )  %>% 
                         str_remove_all( "_log" )  ,
                       nom_type="nom" ) %>% 
  str_replace_all( fixed("\\n"), " " )
names( titre ) <- listvar

y_title <- ""

plotlist <- list()
v <- "dens_2010"
nomvar = "Population density level"
plotlist[[v]] <- trace_forestplot(mod=list_bym2_univ_rw2[[ v ]],
                                  v,
                                  nomvar=nomvar, 
                                  y_title=y_title,
                                  ylim=c(-0.5, 0.6),
                                  axis.title.size = 14)

for(i in seq_along(listvar)){
  v <- listvar[i]
  plotlist[[v]] <- trace_ribbon(tmp=list_bym2_univ_rw2[[ v ]]$summary.random[[v]],
                                v, 
                                titre[v], 
                                y_title=y_title, ylim=c(-0.5, 0.6),
                                axis.title.size = 14
  )
}

# Sorting by increasing DIC
plotlist <- plotlist[res_rw2$v0[!(res_rw2$v0 %in% c("Empty model"))] ]

ggfinal <- ggarrange(  plotlist = plotlist,
                       ncol = 2,
                       nrow = ceiling(length(plotlist)/2), 
                       labels = "AUTO",
                       hjust = -1.2
) %>%
  annotate_figure(  
    left = text_grob( "Contribution of each variable to the linear predictor (logarithm of TB standardized notification rates)", rot = 90, size=14), #, just ='top'
  )
ggfinal

CairoPNG("var_contrib_bym2_univ.png" %>% respath, width=600, height=800)
print(  ggfinal )
dev.off()

postscript("var_contrib_bym2_univ.eps" %>% respath, width=8, height=10)
print(  ggfinal )
dev.off()




#------------ Graphique des contributions des variables dans les modeles univariés avec interaction periode * variable -----
titre <- get_var_name( listvar  %>% 
                         str_remove_all( "_g" )  %>% 
                         str_remove_all( "_log" )  ,
                       nom_type="nom" ) %>% 
  str_replace_all( fixed("\\n"), " " )
names( titre ) <- listvar

y_title <- ""

plotlist <- list()
v <- "dens_2010"
nomvar = "Population density level"
plotlist[[v]] <- trace_forestplot_inter(
  mod=list_bym2_univ_rw2_inter[[v]], 
  v=v,
  nomvar=nomvar,  
  y_title=y_title,
  axis.title.size = 14)

for(i in seq_along(listvar)){
  v <- listvar[i]
  plotlist[[v]] <- trace_ribbon_inter(tmp=list_bym2_univ_rw2_inter[[v]]$summary.random[[v]],
                                      v=v, 
                                      nomvar=titre[v],
                                      y_title=y_title, 
                                      ylim=c(-0.5, 0.6),
                                      axis.title.size = 14)
}

# Sorting by increasing DIC
plotlist <- plotlist[res_rw2$v0[!(res_rw2$v0 %in% c("Empty model"))] ]

ggfinal <- ggarrange(  plotlist = plotlist,
                       ncol = 2,
                       nrow = ceiling(length(plotlist)/2),
                       labels = "AUTO",
                       hjust = -1.2
) %>%
  annotate_figure(  
    left = text_grob( "Contribution of each variable to the linear predictor  (logarithm of TB standardized notification rates)", rot = 90, size=14))
ggfinal

CairoPNG("var_contrib_bym2_univ_inter.png" %>% respath, width=600, height=800)
print(  ggfinal )
dev.off()

ggsave("var_contrib_bym2_univ_inter.pdf" %>% respath, plot=ggfinal, width=8, height=10)


# Ecological analysis of the association between tuberculosis incidence and
# ZIP code-level socioeconomic variables, 2008--2019
# code last revised on 2025-07-18

library(tidyverse)
library(magrittr)
library(sf)
library(Cairo)
library(INLA)
library(ggpubr)
library(corrgram)
library(corrplot)
library(flextable)
library(Hmisc)
library(plotly)
library(RColorBrewer)
library(PerformanceAnalytics)


setEPS() # Option to generate EPS file with the postscript function
set_flextable_defaults(  font.size = 10 )

source( "include/include_functions.R", encoding = "UTF-8" )



#----------- Read TB cases dataset -------------
d_PMSI21 <- arrow::read_parquet("TB_cases_PMSI21.parquet" %>% datapath) 

#---------- Read TB expected rates -------------
Ei_PMSI21 <- arrow::read_parquet("expected_cases_PMSI21.parquet" %>% datapath())

#---------------- Read explanatory variables ---------------
expl_var_PMSI21 <- arrow::read_parquet("explanatory_variables_PMSI21.parquet" %>% datapath())


#-------------------------- Read correspondence files --------------------------
# Variable names
corresp_variable_nom <- read_csv2("corresp_variable_name.csv" %>% datapath) %>%
  column_to_rownames("variable")

# Correspondence French departments (districts) - regions
corresp_dep_reg <- read_delim("depts2016.txt" %>% datapath) %>%
  dplyr::select( REGION, DEP )


#--------- Read shapefiles régions et départements------------
# Import shapefile of PMSI21 codes for metropolitan France
map_PMSI21 <- readRDS( datapath( "France_shapefiles/map_PMSI21_simplifiee_fusionBihorel.rds" ) ) 

map_dep <- readRDS( datapath( "France_shapefiles/map_PMSI21_dep_fusion_Bihorel.rds" ) )  

map_reg <- readRDS( datapath( "France_shapefiles/map_PMSI21_reg_fusion_Bihorel.rds" ) )



#---------- Read the neighboring matrix of PMSI21 codes ---------
# Neighbour list object
tb.adj <- readRDS("neighboring_structure_PMSI21.rds" %>% datapath)

# Create the adjacency matrix
W <- spdep::nb2mat( neighbours = tb.adj, style = "B", zero.policy = TRUE )



#------------------------- Prepare the analysis dataset ------------------------
n_areas <- nrow(map_PMSI21)

# Merge PMSI21 map shapefile with TB data counts + expected counts + explanatory variables
map_PMSI21_2 <- 
  map_PMSI21 %>% 
  cross_join(tibble(period=unique(d_PMSI21$period))) %>%
  left_join(d_PMSI21) %>%
  left_join(Ei_PMSI21) %>%
  left_join(expl_var_PMSI21)  %>%
  arrange(PMSI21_CODE, period)

# Replace missing TB counts with 0 
map_PMSI21_2 %<>%  
  mutate( 
    n = replace_na( n, 0 )
  )

# Adding region codes
map_PMSI21_2 %<>% 
  left_join( corresp_dep_reg %>% rename(reg = "REGION", dep_pmsi="DEP") )




#------------------ Pre-treatment of the analysis dataset ---------------------
# List of all continuous variables
list_var_cont_all <- c( 
  "crowded_house", 
  "unemploy",
  "house_income", 
  "fdep",
  "manual_workers",
  "high_school_grad"
)

# Categorical variable (density level)
list_var_cat0 <- c( 
  "dens_2010" 
)


# Verify there is no NA in the explanatory variables
map_PMSI21_2 %>% 
  st_drop_geometry() %>%
  summarise(across( all_of( list_var_cont_all ), function( x ){ return( sum( is.na(x)) )} ))


# Log transformation of all the continuous variables except the FDEP
# NB : function mylog defined in include/ to deal with 0 in the original variables
listvar_log0 <- list_var_cont_all[list_var_cont_all != "fdep"] 
listvar_log <- paste(listvar_log0, "log", sep="_")
  
map_PMSI21_2 %<>% 
  mutate( across( .cols=all_of(listvar_log0), .fns=mylog, .names="{.col}_log" )) 

# Standardize (mean 0, sd 1) the continuous variables by period
map_PMSI21_2 %<>% 
  group_by(period) %>%
  mutate(across(.cols=all_of(c(list_var_cont_all, listvar_log)), .fns=cr, .names="{.col}_cr" )) %>%
  ungroup()  



#------------------ Histograms of the variable distributions -------------------
# Log and not log
tmp <- map_PMSI21_2[, c(str_glue("{list_var_cont_all}_cr"), str_glue("{listvar_log}_cr"), list_var_cat0 ) ] %>%
  st_drop_geometry()

names(tmp) <- corresp_variable_nom[colnames(tmp) %>% str_remove("_cr"), ]$name

plotlist <- list()
for( v in names(tmp)[names(tmp) != "dens_2010"] ){
  cat( "\n\n", v, "\n" )
  
  plotlist[[v]] <-  
    ggplot( tmp, aes(x =!!sym(v) ) ) + 
    geom_histogram()  +
    scale_fill_grey() + 
    theme_minimal() +
    labs( y="" )
}

v = "Population density level"
plotlist[[v]] <- ggplot( tmp, aes(x =!!sym(v) ) ) +
  geom_bar() + 
  theme_minimal() +
  labs( y="" )
  
# Histogramme
ncol = 2
nrow = (length(plotlist) / ncol) %>% ceiling
ggfinal <- ggarrange(  plotlist = plotlist[sort(names(plotlist))],
                       ncol = ncol,
                       nrow = nrow,
                       common.legend = TRUE
) %>%
  annotate_figure(  
    left = text_grob( "Count", rot = 90, size=10 ))

ggfinal

png(file="distribution_variables.png" %>% respath(), width=600, height=900)
ggfinal
dev.off()


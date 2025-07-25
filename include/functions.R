# Fonctions utiles
library( tidyverse )
library( magrittr )
library( sf )
library( readr )
library( data.table )
library(yaml)

conf <- read_yaml('conf/global_conf.yml')

#---------- Function that adds the data path to a filename------
datapath <- function(file){
  return( paste0(conf$data_path, file))
}

#---------- Function that adds the result path to a filename------
respath <- function(file){
  return( paste0(conf$res_path, file))
}


#---- Fonction pour renommer les éléments d'une liste -----
add_suffix_to_name <- function(mylist,  suffix='_inter') {
  names(mylist) <- paste0(names(mylist), suffix)
  return(mylist)
}

#---- Fonction pour centrer réduire -----
cr <- function( x ){
  return( ( x - mean( x, na.rm=TRUE ) ) / sd( x, na.rm=TRUE ) )
}


#---- Fonction pour mutliplier par un sd et rajouter un mean -----
cr_inv <- function( x, mymean, mysd ){
  return( x * mysd + mymean )
}


myround <- function( d, digits=2 ){
  sapply( d, FUN=function(x){ formatC( x, digits=digits, format="f", flag = "0" ) } ) %>% unname()
}


#----- Transforme un vecteur en log en ajoutant un petit increment xc si il y a des 0
mylog <- function(x){
  # https://aosmith.rbind.io/2018/09/19/the-log-0-problem/
  #square of the first quartile divided by the third quartile
  if(any(x==0)){
    xc <- sqrt(quantile(x, probs=0.25) / quantile(x, probs=0.75))
    x <- x + xc
  }  
  return(log(x))
}



#----- Correspondance function variable name - variable label-----
get_var_name <- function( v_vect , tab_corresp = corresp_variable_nom, nom_type="nom_court" ){
  tab_corresp <- tab_corresp %>% 
    full_join( tibble( v=v_vect ) ) %>%
    mutate( nom = ifelse( is.na( nom ), v, nom ),
            nom_court = ifelse( is.na( nom_court ), v, nom_court )
    )

  tt <- tab_corresp %>% dplyr::slice( match( v_vect, tab_corresp$v  ) )
  
  return( tt[ , nom_type ] %>% unlist() %>% unname )
}


#------------------- Inter-quantile range ---------------------------
get_IqValues <- function( x, pr=probs, weights=NULL, type="quantile" ){
  if( is.null( weights ) ){
    res <- x %>% quantile( probs=pr )
  } else {
    res <- x %>% wtd.quantile( weights=tmp$pop, probs=pr, type=type )
  }
  return( res )
}


get_IqRR <- function( x, pr=probs, weights=NULL, type="quantile" ){
  return( get_IqValues(  x, pr, weights, type ) %>% diff %>% exp %>% unname )
}




#----------------- Function to deal with INLA outputs --------------------------



# create a margin marge by sampling from the marginal distribution of variances of a parameter v : 
# https://gkonstantinoudis.github.io/INLA/StrokeSheffield.html
# https://www.paulamoraga.com/book-geospatial/sec-inla.html
get_marginal_var_sample <- function( modele, v="ID_PMSI21" ){
  marg.hyper <- modele$marginals.hyperpar[[ str_glue( "Precision for {v}" ) ]]
  marg <-  inla.tmarginal( 
    function(x) 1/x,
    marg.hyper
  )
  return( marg )
}

# Echantilloner dans la distribution marginale des sd d'un paramètre v : créer une marge 
get_marginal_sd_sample <- function( modele, v="ID_PMSI21" ){
  marg_sd <-  inla.tmarginal( 
    function(x) sqrt( 1/x ),
    modele$marginals.hyperpar[[ str_glue( "Precision for {v}" ) ]]
  )
  return( marg_sd )
}

# Calculer la moyenne sur une marge (échantillon issu de la distriubtion marginale d'un paramètre )
get_marginal_mean <- function( marge ){
  inla.emarginal( function(x) x, marge  )
}







#-------- Table de résultats pour les modeles : WAIC, variance, proportional change in variance (PCV) --
# @ model.list : liste de modeles
# @spatial_var_ref : la varaince spataile de référence (modele vide)
get.res.table <- function( model.list = list.prepend( model_bym2[ c( "dep", list_var ) ], model_bym2_vide ), 
                           model.list.names = NULL,
                           digits=2,
                           spatial_var_ref = m_vide
){
  if( is.null( model.list.names ) ){
    model.list.names <- names( model.list )
  }
  
  res <- tibble(
    v0 =  model.list.names,
    WAIC = sapply( model.list, FUN=function(x){ return( round( x$waic$waic ) ) } ),
    DIC = sapply( model.list, FUN=function(x){ return( round( x$dic$dic ) ) } ),
    spatial_var =  sapply( model.list, FUN=get.spatial.var )  ,
    phi = str_glue( "{sapply( model.list, FUN=get.phi )} [{sapply( model.list, FUN=get.IC.phi )}]" )
    #variable_prec = mapply( FUN=get.prec, model.list, names( model.list ), SIMPLIFY = FALSE  ) %>% unlist,
  ) %>%
    mutate( 
      percent_change_var = round( ( spatial_var - spatial_var_ref ) / spatial_var_ref  * 100, digits ) ,
      spatial_var = round( spatial_var, digits=digits ),
      spatial_var = str_glue( "{ spatial_var } [{sapply( model.list, FUN=get.IC.spatial.var )}]" )
       
      
      ) %>%
    arrange( DIC ) %>%
    dplyr::select(  v0, DIC, WAIC, spatial_var, percent_change_var, phi, everything() )
  
  res
} 



#-------- Table de résultats pour les modeles : DIC --
# @ model.list : liste de modeles
get.DIC.table <- function( model.list, 
                           model.list.names = NULL,
                           digits=2
){
  if( is.null( model.list.names ) ){
    model.list.names <- names( model.list )
  }
  
  res <- tibble(
    v0 =  model.list.names,
    DIC = sapply( model.list, FUN=function(x){ return( round( x$dic$dic ) ) } ),
    WAIC = sapply( model.list, FUN=function(x){ return( round( x$waic$waic ) ) } )
  ) %>%
    arrange( DIC ) 
  
  res
} 






#------------- Function to create a linear univariate BYM2 model ------------
get_univ_model <- function(     mydata = map_PMSI21_2,
                                Ei = "Ei_knn",
                                form0 = 'n ~ 
                                  f( ID_area, 
                                     model="bym2",
                                     graph = W,
                                     scale.model = TRUE,
                                     constr = TRUE,
                                     hyper = prior.prec.bym2
                                  )',
                                form_v = 'f( {v}, 
                                 model = "rw2", 
                                 scale.model = TRUE, 
                                 constr = TRUE
                              )'
){

  
  form <- str_glue( '{form0} + {form_v}' ) %>% 
    str_glue %>% 
    str_replace_all("[\r\n]", "") %>% 
    as.formula()
  
  print( form )

  model_bym2_rw2 <- inla(
    form,
    family = "poisson",
    data = mydata,
    E = mydata[[Ei]],
    control.predictor = list( compute=TRUE ),
    control.compute = list( config=TRUE, waic=TRUE, dic=TRUE, cpo=FALSE, openmp.strategy="huge"),
    verbose=TRUE
  )

  return( model_bym2_rw2 )
}






#---- carto avec ggplot ------

#------- Fonction qui crée une carte de valeurs catégorielles avec échelle discrétisée avec ggplot ------
plotmap_cat <- function( mymap = map_PMSI21_2  %>% filter( annee == "2010" ), 
                          mymap_dep = map_dep, 
                          mymap_reg = map_reg,
                          var="dens_2016",  
                          mylabels = NULL, 
                          mycolors = RColorBrewer::brewer.pal( n = 3, name = "Blues" ),
                          col.dep="darkgrey",  #"#00554b", 
                          #legend.position = c(0.9, 0.5), 
                          legend.position = "bottom", 
                          legend.direction = "horizontal",
                          legend.title="Population density in 2016",
                          mytitle = "",
                          size = 0,
                          size.dep = 0.5,
                          show.legend = TRUE,
                          legend.text.size = 10
                          #color.scale = "gradient"
){
  
  
  # Labels
  if( is.null( mylabels ) & is.factor( mymap[[ var ]]  ) ) {
    mylabels <- levels( mymap[[ var ]] )
  }
  
  # Faire une échelle de couleurs discrète manuelle dans ggplot : https://stackoverflow.com/questions/65947347/r-how-to-manually-set-binned-colour-scale-in-ggplot
  
  
  gg <- ggplot() + 
    geom_sf( data = mymap, aes(fill =  .data[[ var ]] ), color=NA, size=size, show.legend = show.legend 
    )  +
  scale_fill_manual( values = mycolors, drop = FALSE, na.value = "transparent" )
  
  gg <- gg  +
    geom_sf( data = mymap_dep, color=col.dep, fill=NA, size = size.dep ) +
    theme_void() + 
    theme( legend.position = legend.position,
           legend.direction = legend.direction,
           plot.title = element_text(hjust = 0.5),
           legend.text = element_text( size = legend.text.size ),
           legend.title = element_text( size = legend.text.size )
    ) +
    labs( fill=legend.title ) +
    ggtitle( mytitle )
  
  if( !is.null( mymap_reg ) ){
    gg <- gg + geom_sf( data = mymap_reg, color="black", fill=NA, size = size.dep*3 ) 
  }
  
  return( gg )
}



#------- Fonction qui crée une carte de valeurs continues avec échelle discrétisée avec ggplot ------
plotmap_cont <- function( mymap = map_PMSI21_2, 
                     mymap_dep = map_dep, 
                     mymap_reg = map_reg,
                     var="rapport", 
                     mybreaks = breakvalues, 
                     mylabels = mybreaks, 
                     mylimits = range( mybreaks ),
                     mycolors = colorvalues,
                     show.limits = FALSE, 
                     col.dep="darkgrey",  
                     legend.position = "bottom", 
                     legend.direction = "horizontal",
                     legend.title="SIR",
                     mytitle = "",
                     size = 0,
                     size.dep = 0.5,
                     color.scale = "gradient",
                     show.legend = TRUE,
                     legend.text.size = 10
                     ){
  # Faire une échelle de couleurs discrète manuelle dans ggplot : https://stackoverflow.com/questions/65947347/r-how-to-manually-set-binned-colour-scale-in-ggplot
  
  if( class( mybreaks ) == "character" ){
    mybreaks <- do.call( mybreaks, list( x = mymap[[ var  ]], na.rm=TRUE )  ) %>% as.numeric() %>% round( 2 )
    mylabels <- mybreaks
  }
  
  gg <- ggplot() + 
    geom_sf( data = mymap, aes(fill =  .data[[ var ]] ), color=NA, size=size, show.legend = show.legend 
    ) 
  
  
  if( color.scale == "binned" ){
    gg <- gg  +
      binned_scale(aesthetics = "fill",
                   scale_name = "stepsn",
                   palette = function(x) mycolors,
                   breaks = mybreaks,
                   limits = mylimits,
                   show.limits = show.limits,
                   guide = "colorsteps"
      ) 
  } else {
     gg <- gg +
        scale_fill_gradientn( colours=mycolors,
                         na.value = "transparent",
                         limits = range( mybreaks ),
                         breaks = mybreaks,
                         labels = mylabels,
                         oob = scales::squish )
     
  }
  
   gg <- gg  +
    geom_sf( data = mymap_dep, color=col.dep, fill=NA, size = size.dep ) +
    theme_void() + 
    theme( legend.position = legend.position,
           legend.direction = legend.direction,
           plot.title = element_text(hjust = 0.5),
           legend.text = element_text( size = legend.text.size ),
           legend.title = element_text( size = legend.text.size )
           ) +
    labs( fill=legend.title ) +
    ggtitle( mytitle )
  
  if( !is.null( mymap_reg ) ){
    gg <- gg + geom_sf( data = mymap_reg, color="black", fill=NA, size = size.dep*3 ) 
  }
  
  return( gg )
}


#---- Fonction de choroplethe pour variable classe ou continue -------
plotmap <- function( 
                          var=variables[1], #"rapport", 
                          mytitle = "",
                          legend.title =  "SIR",
                          mymap = map_PMSI21_2, 
                          mymap_dep = map_dep, 
                          mymap_reg = map_reg,
                          mybreaks = NULL, # = breakvalues, 
                          mylabels = mybreaks, 
                          mylimits = range( mybreaks ),
                          show.limits = FALSE, 
                          mycolors = colorvalues,
                          col.dep = "darkgrey",  
                          legend.position = "bottom", 
                          legend.direction = "horizontal",
                          
                          size = 0,
                          size.dep = 0.5,
                          color.scale = "gradient",
                          show.legend = TRUE,
                          legend.text.size = 10
){
  
  if( is.numeric( mymap[[ var ]] )  ){
    gg <- plotmap_cont(
      mymap = mymap, 
      mymap_dep = mymap_dep, 
      mymap_reg = mymap_reg,
      var = var,
      mybreaks = mybreaks, 
      mylabels = mylabels,  
      mylimits = mylimits, 
      show.limits = show.limits, 
      mycolors = mycolors,
      col.dep = col.dep, 
      legend.position = legend.position, 
      legend.direction = legend.direction,
      legend.title = legend.title,
      mytitle = mytitle,
      size = size,
      size.dep = size.dep,
      color.scale = color.scale,
      show.legend =  show.legend,
      legend.text.size = legend.text.size
    )
    
    } else {
      gg <- plotmap_cat(
        mymap = mymap, 
        mymap_dep = mymap_dep, 
        mymap_reg = mymap_reg,
        var = var,
        mylabels = mylabels, 
        mycolors = mycolors,
        col.dep = col.dep, 
        legend.position = legend.position, 
        legend.direction = legend.direction,
        legend.title = legend.title,
        mytitle = mytitle,
        size = size,
        size.dep = size.dep,
        show.legend =  show.legend,
        legend.text.size = legend.text.size
      )
    }
  
  return( gg )
  
}



plotmap_inset <- function(  var= "rapport", 
                            mytitle="test", 
                            mymap=map_PMSI21_2, 
                            mymap_dep = map_dep,
                            mymap_reg = map_reg,
                            mycolors=colorvalues, 
                            mybreaks=breakvalues,  
                            mylabels = mybreaks,
                            size = 0, 
                            size.dep = 0.5, 
                            col.dep="darkgrey",
                            legend.title = "SIR",
                            color.scale = "gradient",
                            legend.position = "bottom", 
                            legend.direction = "horizontal",
                            legend.text.size = 10
                            ){
  
  
  #---- Si variable continue et fonction de création de breaks ------
  if( is.numeric( mymap[[ var ]] ) & class( mybreaks ) == "character" ){
    
    mybreaks <- mapsf::mf_get_breaks(x=mymap[[ var  ]], nbreaks=length(mycolors), breaks=mybreaks ) %>% round( 2 )
    
    mylabels <- mybreaks
  }
  
  ggm1 <- plotmap( var = var, 
                   mymap = mymap,
                   mymap_dep = mymap_dep,
                   mymap_reg = mymap_reg,
                   mycolors = mycolors, 
                   mybreaks = mybreaks,
                   mytitle = mytitle,
                   size=size,
                   size.dep=size.dep,
                   col.dep=NA, #col.dep,
                   legend.title = legend.title,
                   mylabels = mylabels,
                   color.scale = color.scale,
                   legend.position = legend.position,
                   legend.direction = legend.direction,
                   legend.text.size = legend.text.size
                   ) 
  ggm2 <- plotmap(  var = var, 
                    mymap = mymap %>% filter( dep_pmsi %in% c( 75, 92, 93, 94 )  ),
                    mymap_dep = mymap_dep %>% filter( dep_pmsi %in% c( 75, 92, 93, 94 ) ),
                    mymap_reg = NULL,
                    mycolors = mycolors,
                    mybreaks = mybreaks,
                    legend.position = "none",
                    size=size,
                    size.dep=size.dep,
                    col.dep=col.dep,
                    mylabels = mylabels,
                    color.scale = color.scale,
                    legend.text.size = legend.text.size
                    )  + theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1))
  
  gg_inset_map1 <- ggdraw(ggm1) +
    #draw_plot( ggm1 ) +
    draw_plot( ggm2, x = 0.7, y = 0.8, width = 0.2, height = 0.2) 
  
  gg_inset_map1
}




library(ggside)

trace_ribbon <- function(tmp = mod$summary.random[[v]], v, nomvar, y_title, ylim, 
                         axis.title.size=12,
                         axis.text.size=12){
  
  ggplot(data = tmp, aes(x=`ID`)) +
    geom_ribbon(aes(ymin=`0.025quant`, ymax=`0.975quant`), fill = "grey70") +
    geom_line( aes(y=mean), linewidth=0.8) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    theme_classic() +
    xlab(str_glue("{nomvar}*")) +
    coord_cartesian(ylim=ylim) +
    theme(axis.title = element_text(size=axis.title.size), 
          axis.text = element_text(size=axis.text.size)) +
    geom_xsideboxplot(
      data = map_PMSI21_2 %>% 
        dplyr::select(all_of(c("myvar" = v))) %>% st_drop_geometry(), 
      aes(x=myvar),
      orientation = "y",
      outlier.size = 1) +
    ylab(y_title) + 
    theme_ggside_void(base_line_size=0) +
    scale_y_continuous() + scale_xsidey_discrete()
}


trace_ribbon_inter <- function(tmp = mod$summary.random[[v]], 
                               v,
                               nomvar, 
                               y_title, ylim,
                               axis.title.size =12,
                               axis.text.size=12){
  
  n <- nrow(tmp) / 2
  tmp1 <- tmp[1:n,]
  tmp2 <- tmp[n+(1:n),]
  
  ggplot(data = tmp1, aes(x=`ID`)) +
    geom_ribbon(aes(ymin=`0.025quant`, ymax=`0.975quant`), fill = "grey70", alpha=0.5) +
    geom_line( aes(y=mean), linewidth=0.8) +
    geom_ribbon(aes(ymin=`0.025quant`, ymax=`0.975quant`), fill = "pink", data = tmp2, alpha=0.5) +
    geom_line( aes(y=mean), linewidth=0.8, data = tmp2, col="red") +
    
    geom_hline(yintercept = 0, linetype = "dashed") +
    theme_classic() +
    xlab(str_glue("{nomvar}*")) +
    coord_cartesian(ylim=ylim) +
    theme(axis.title = element_text(size=axis.title.size), 
          axis.text = element_text(size=axis.text.size)) +
    geom_xsideboxplot( 
      data = map_PMSI21_2 %>% 
        dplyr::select(all_of(c("myvar" = v))) %>% st_drop_geometry(), 
      aes(x=myvar),
      orientation = "y",
      outlier.size = 1) + 
    ylab(y_title) + 
    theme_ggside_void(base_line_size=0) +
    scale_y_continuous() + scale_xsidey_discrete()
}




# Ribbon plot for the central analysis and the two sensitivity analysis models 
trace_ribbon_sensitivity <- function(tmp, 
                                     col_models = c("black", "red", "#14a45c"),
                                     fill_models = c("grey70", "pink", "lightgreen"),
                                     v, nomvar, 
                                     y_title, ylim,
                                     axis.title.size = 12,
                                     axis.text.size = 12
                                     ){

  ggplot(data = tmp, aes(x=`ID`)) +
    geom_ribbon(aes(ymin=`0.025quant`, ymax=`0.975quant`, fill = model), alpha=0.5) +
    geom_line( aes(y=mean, col = model), linewidth=0.8) +
    scale_fill_manual(values = fill_models) +
    scale_color_manual(values = col_models) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    theme_classic() +
    xlab("{nomvar}*" %>% str_glue) + 
    coord_cartesian(ylim=ylim) +
    theme(axis.title = element_text(size=axis.title.size), 
          axis.text = element_text(size=axis.text.size)) +
    geom_xsideboxplot( 
      data = map_PMSI21_2 %>% 
        dplyr::select(all_of(c("myvar" = v))) %>% st_drop_geometry(), 
      aes(x=myvar),
      orientation = "y",
      outlier.size = 1) + 
    ylab(y_title) + 
    theme_ggside_void(base_line_size=0) +
    scale_y_continuous() + scale_xsidey_discrete() +
    theme(legend.position = "none")
}



trace_forestplot <-  function(mod, v, nomvar, y_title, ylim,
                              axis.title.size = 12,
                              axis.text.size = 12){
  tmp <- mod$summary.fixed %>%
    rownames_to_column("var") %>%
    dplyr::filter( str_detect(var, v) ) %>%
    mutate(var = str_remove(var, v)) %>%
    tibble()
  
  tmp <- tibble(var = map_PMSI21_2[[v]] %>% levels ) %>%
    left_join( tmp ) %>%
    replace_na(list(mean=0, `0.025quant`=0, `0.975quant`=0)) %>%
    mutate( var = var %>% factor(levels=levels(map_PMSI21_2[[v]])))  
  
  # Densité : forest plot 
  tmp |>
    ggplot(aes(x = var)) + 
    theme_classic() +
    geom_point(aes(y=mean), shape=15, size=2) +
    geom_linerange(aes(ymin=`0.025quant`, ymax=`0.975quant`), na.rm = FALSE) +
    geom_hline(yintercept = 0, linetype="dashed") +
    xlab(nomvar) +
    ylab(y_title) + 
    coord_cartesian(ylim=ylim) +
    theme(axis.title = element_text(size=axis.title.size), 
          axis.text = element_text(size=axis.text.size)) +
    geom_xsidecol(data = map_PMSI21_2 %>% 
                    dplyr::select(all_of(c("myvar" = v))) %>%
                    st_drop_geometry() %>%
                    count(myvar),
                  aes(x=myvar, y=n)) + 
    theme_ggside_void(base_line_size=0) 
}


# Tracer le forestplot pour une variable catégorielle en interaction avec la période
trace_forestplot_inter <-  function(mod, v, nomvar, y_title,
                                    axis.title.size = 12,
                                    axis.text.size = 12){
  tmp <- mod$summary.fixed %>%
    rownames_to_column("var") %>%
    dplyr::filter( str_detect(var, v) ) %>%
    mutate(
      var = str_remove(var, v),
      var = str_remove(var, "annee"),
      annee = str_extract(var, "\\d+"),
      var = var %>%  str_extract(regex("[a-zA-z]+"))
    ) %>%
    tibble()

  # Densité : forest plot 
  pos <- position_dodge(width = 0.5)
  
  tmp |>
    ggplot(aes(x = var, col=annee)) + 
    theme_classic() +
    geom_point(aes(y=mean), shape=15, size=2, position=pos) +
    geom_linerange(aes(ymin=`0.025quant`, ymax=`0.975quant`), na.rm = FALSE, position=pos) +
    xlab(nomvar) + 
    ylab(y_title) + 
    theme(axis.title = element_text(size=axis.title.size), 
          axis.text = element_text(size=axis.text.size),
          legend.position = "none"
    ) +
    geom_xsidecol(data = map_PMSI21_2 %>% 
                    dplyr::select(all_of(c("myvar" = v))) %>%
                    st_drop_geometry() %>%
                    count(myvar),
                  aes(x=myvar, y=n), col="grey") + 
    theme_ggside_void(base_line_size=0) +
    scale_color_manual(values=c("2010"="black", "2016"="red"))

}


# Tracer le forestplot for the central analysis + an alternate model (sensitivty analysis)
# @models : model list (length 2 ?)
trace_forestplot_sensitivity <-  function(models=list(model_bym2_chom_suroc_rev_dens, 
                                                      model_bym2_chom_suroc_rev_dens_ventil, 
                                                      model_bym2_chom_suroc_rev_dens_RF),
                                          name_models=c("kNN", "Pro rata", "RF"), 
                                          col_models = c("black", "red", "darkgreen"),
                                          v, 
                                          nomvar, 
                                          y_title,
                                          ylim = c(-0.5, 0.6),
                                          axis.title.size = 12,
                                          axis.text.size = 12
                                          ){
  
  tmp_all <- NULL
  
  for(i in seq_along(models)){
    print(i)
    tmp <- models[[i]]$summary.fixed %>%
      rownames_to_column("var") %>%
      dplyr::filter( str_detect(var, v) ) %>%
      mutate(var = str_remove(var, v)) %>%
      tibble()
    
    tmp <- tibble(var = map_PMSI21_2[[v]] %>% levels ) %>%
      left_join( tmp ) %>%
      replace_na(list(mean=0, `0.025quant`=0, `0.975quant`=0)) %>%
      mutate( var = var %>% factor(levels=levels(map_PMSI21_2[[v]])))
    
    tmp$model <- name_models[i]
    
    tmp_all <- bind_rows(tmp_all, tmp) 
  }
  
  # Densité : forest plot 
  pos <- position_dodge(width = 0.5)

  tmp_all |>
    ggplot(aes(x = var, col=model)) + 
    theme_classic() +
    geom_point(aes(y=mean), shape=15, size=2, position=pos) +
    geom_linerange(aes(ymin=`0.025quant`, ymax=`0.975quant`), na.rm = FALSE, position=pos) +
    geom_hline(yintercept = 0, linetype="dashed") +
    xlab(nomvar) + 
    ylab(y_title) + 
    coord_cartesian(ylim=ylim) +
    theme(axis.title = element_text(size=axis.title.size), 
          axis.text = element_text(size=axis.text.size),
          legend.position = "none"
    ) +
    geom_xsidecol(data = map_PMSI21_2 %>% 
                    dplyr::select(all_of(c("myvar" = v))) %>%
                    st_drop_geometry() %>%
                    count(myvar),
                  aes(x=myvar, y=n), col="grey") + 
    theme_ggside_void(base_line_size=0) +
    scale_color_manual(values=col_models, labels=name_models)
}





#-------- Graphique des contributions des variables dans ce modele BYM2 multivarié ------------
graphic_contrib_mult <- function(modele_name = "model_bym2_chom_suroc_rev_dens"){
  modele <- get(modele_name)
  
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
  
  i <- 1
  plotlist[[v]] <- trace_forestplot(modele,
                                    v,
                                    nomvar,  
                                    y_title, 
                                    ylim=c(-0.4, 0.6),
                                    axis.title.size = 14)
  
  for(i in seq_along(listvar)){
    myletter <- letters[i+1] 
    v <- listvar[i]
    plotlist[[v]] <- trace_ribbon(tmp=modele$summary.random[[v]],
                                  v, 
                                  titre[v], 
                                  y_title=y_title, 
                                  ylim=c(-0.4, 0.6),
                                  axis.title.size = 14)
  }
  
  ggfinal <- ggarrange(  plotlist = plotlist,
                         ncol = 2,
                         nrow = ceiling(length(plotlist)/2),
                         labels = "AUTO",
                         hjust = -1.2
  ) %>%
    annotate_figure(  
      left = text_grob( "Contribution of each variable to the linear predictor (log of TB standardized notification rates)", rot = 90, size=14))
  ggfinal
  
  CairoPNG("resultat/article/var_contrib_{modele_name}.png" %>% str_glue, width=600, height=600)
  print(  ggfinal )
  dev.off()
  
  postscript("resultat/article/var_contrib_{modele_name}.eps" %>% str_glue(), width=8, height=8)
  print(  ggfinal )
  dev.off()
}
  



#--------------------- Create Inter-quanile Rate Ratios --------------------
get_IqRR <- function(mysample, v, list_q){
  
  tt <- mysample$latent %>% 
    as.data.frame() %>%
    rownames_to_column(var="name") %>%
    filter(str_detect(name, v)) 
  
  if(is.numeric(map_PMSI21_2[[v]])){
    tt %<>%
      mutate(ID = map_PMSI21_2[[v]] %>% unique %>% sort) 
  } else if(is.factor(map_PMSI21_2[[v]])) {
    tt %<>%
      mutate(ID = name %>% str_remove_all(v) %>% str_remove_all(":1") ) 
    
  } else {
    stop("The column v must be either numeric or factor")
  }
  
  tt %<>%
    filter(ID %in% list_q[[v]]) 
  
  if(nrow(tt) > 1){
    IqRR <- tt %>%
      summarise(IqRR = diff(V1) %>% exp)  %>%
      extract(1,1)
  } else {
    IqRR <- tt$V1 %>% exp
  }
  
  return(IqRR)
}

get_IqRR_stats <- function(v, samples, list_q){
  res <- map_dbl(.x = samples, .f=get_IqRR, v=v, list_q=list_q)
  
  tibble(
    variable = v,
    IqRR_mean = mean(res), 
    IqRR_sd = sd(res), 
    IqRR_0.025= quantile(res, 0.025),
    IqRR_0.975= quantile(res, 0.975))
}

#--------------- IQRR modèle multivarié ----------------------
IQRR_mult <- function(modele_name = "model_bym2_chom_suroc_rev_dens"){

  modele <- get(modele_name)
  samples <- inla.posterior.sample(100, modele, verbose=TRUE )
  
  # IqRR variables continues
  listvar <- c("perc_chom_log_g", "perc_suroc_log_g", "revenu_median_log_g") 
  res_multi <- map(listvar, get_IqRR_stats, samples=samples, list_q=list_q) %>% 
    bind_rows()
  
  # IqRR variables catégorielles
  res_multi <- bind_rows(res_multi,  get_IqRR_stats("dens_2010", samples, list_q) )
  
  res_multi %<>%
    arrange(-abs(1-IqRR_mean)) %>%
    mutate(IC = str_glue("({IqRR_0.025 %>% formatC(format='f', digits=2, flag='0')}; {IqRR_0.975 %>% formatC(format='f', digits=2, flag='0')})"),
           IqRR_mean =  formatC(IqRR_mean, format='f', digits=2, flag='0')
    ) %>%
    mutate( v = str_remove_all( variable, "_g|_cr|_log" )) %>% 
    left_join( corresp_variable_nom ) %>%
    mutate(nom = str_replace_all(nom, fixed("\\n"), " "))
  
  res_multi2 <- res_multi %>% 
    mutate("IdRR (95 CrI)" = paste(IqRR_mean, IC)) %>%
    select(Variable=nom, `IdRR (95 CrI)`)

  names(res_multi2)[2] <-  "IdRR (95% CrI)"
  
  tt <- qflextable(res_multi2) 
  tt
  
  fich <- "resultat/article/IqRR_{modele_name}_90.docx" %>% str_glue
  cat("\n---------- Writing IQRR in file {fich}----------------------\n" %>% str_glue())
  save_as_docx( tt, 
                values = NULL, 
                path=fich, 
                pr_section = prop_section(page_size = page_size(orient = "portrait")), 
                align = "center")
  
  return(res_multi)
}


# Capitalize the first letter
cap <- function( s ) {
  paste( toupper( substring(s, 1, 1) ),  substring(s, 2), sep="" )
}

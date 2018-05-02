#' Plot a set of PARAFAC models to compare the single components
#'
#' @description Three plots are returned:
#' \enumerate{
#'     \item plot of numer of components vs. model fit
#'     \item plot of different components as colour maps
#'     \item plot of different components as peak lines
#' }
#' The plots are intended to help with a suitable number of components.
#'
#' @param pfres list of several objects of class parafac
#'
#' @return 3 objects of class ggplot
#' @export
#'
#' @seealso \code{\link[staRdom]{eempf_fits}}, \code{\link[staRdom]{eempf_plot_comps}}
#'
#' @examples
#' \donttest{
#' data(pfres_comps1)
#'
#' eempf_compare(pfres_comps)
#' }
eempf_compare <- function(pfres){
  p1 <- eempf_fits(pfres)
  p2 <- eempf_plot_comps(pfres,type=1)
  p3 <- eempf_plot_comps(pfres,type=2)
  p1 %>% print()
  p2 %>% print()
  p3 %>% print()
  return(list(p1,p2,p3)) %>% invisible()
}

#' Fits vs. components of PARAFAC models are plotted
#'
#' @param pfres list of objects of class parafac
#'
#' @return object of class ggplot
#' @export
#'
#' @import ggplot2
#'
#' @examples
#' data(pfres_comps1)
#'
#' eempf_fits(pfres_comps)
eempf_fits <- function(pfres){
  pl <- data.frame(comps=lapply(pfres,"[[","A") %>% lapply(ncol) %>% unlist(),
                   fit=lapply(pfres,"[[","Rsq") %>% unlist()) %>%
    ggplot(aes(x=comps,y=fit))+
    geom_point(size=5,shape=4)
  #pl %>% print()
  pl
}

#' Plot all components of PARAFAC models
#'
#' @description The components can be plottet in two ways: either as a colour map or as two lines (emission, excitation wavelengths) intersecting at the component maximum.
#'
#' @param pfres list of PARAFAC models
#' @param type 1 for a colour map and 2 for peak lines
#'
#' @return object of class ggplot
#' @export
#'
#' @import dplyr
#' @import ggplot2
#' @importFrom grDevices rainbow
#'
#' @examples
#' data(pfres_comps1)
#'
#' eempf_plot_comps(pfres_comps)
#' eempf_plot_comps(pfres_comps,type=2)
eempf_plot_comps <- function(pfres,type=1){
  #pf_fits <- cp_out
  c <- pfres %>% lapply(eempf_comp_mat)
  tab <- lapply(c,function(c1){
    nc1 <- length(c1)
    nc2 <- 0
    lapply(c1,function(c2){
      nc2 <<- nc2 + 1
      c2 <- c2 %>%
        mutate(comps = nc1, comp = paste0("Comp.",nc2))
    }) %>%
      bind_rows()
  }) %>%
    bind_rows()
  fill_max <- tab$value %>% max(na.rm=TRUE)
  vals <- seq(from=0,to=fill_max,length.out = 51)
  vals <- (vals - min(vals))/diff(range(vals))
  breaksx <- tab %>% select(ex) %>% unique() %>% filter(as.numeric(ex)%%50 == 0) %>% unlist()
  breaksy <- tab %>% select(em) %>% unique() %>% filter(as.numeric(em)%%50 == 0) %>% unlist()
  if(type==2){
    breaks <- c(breaksx,breaksy) %>% unique() %>% sort()
    #options(warn=-1)
    pl <- tab %>%
      group_by(comps,comp) %>%
      mutate(max_pos = which.max(value), max_em = em[max_pos],max_ex = ex[max_pos]) %>%
      mutate(exn = ifelse(em == max_em,ex,NA),emn = ifelse(ex == max_ex,em,NA)) %>%
      filter(!is.na(emn) | !is.na(exn)) %>%
      ungroup() %>%
      ggplot()+
      geom_line(aes(x=exn,y=value),colour="lightblue",group="excitation", na.rm=TRUE)+
      geom_line(aes(x=emn,y=value),colour="darkblue",group="emission", na.rm=TRUE)+
      facet_grid(comp ~ comps)+
      scale_x_discrete(breaks = breaks) +
      labs(x="wavelength")
  } else {
    pl <- tab %>%
      ggplot()+
      geom_raster(aes(x=ex,y=em,fill=value))+ #,interpolate=TRUE
      scale_fill_gradientn(colours=rainbow(75)[51:1],values=vals,limits = c(tab$value %>% min(na.rm=TRUE),fill_max))+
      scale_x_discrete(breaks = breaksx) +
      scale_y_discrete(breaks = breaksy) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1))+
      facet_grid(comp ~ comps)
  }
  #print(pl)
  pl
}


#' Plot leverage of emission wavelengths, excitation wavelengths and samples.
#'
#' @description The labels to be shown can be selected by adding the quatile of samples with highest leverages to be labeled.
#'
#' @param cpl leverage, outpout from \code{\link[staRdom]{eempf_leverage}}
#' @param qlabel optional, quantile of which labels are shown (1 = all, 0 = no labels)
#'
#' @return ggplot
#' @export
#'
#' @import ggplot2
#' @import dplyr
#'
#' @seealso \code{\link[staRdom]{eempf_leverage_ident}}
#'
#' @examples
#' data(pfres_comps1)
#'
#' leverage <- eempf_leverage(pfres_comps[[3]])
#' eempf_leverage_plot(leverage)
eempf_leverage_plot <- function(cpl,qlabel=0.1){
  breaks <- cpl$x[cpl$mode != "sample"] %>% as.numeric() %>% na.omit() %>% unique() %>% .[.%%50 == 0]
  pl <- eempf_leverage_data(cpl,qlabel=qlabel) %>%
    ggplot(aes(x=x,y=leverage))+
    geom_point(alpha = 0.4)+
    geom_text(aes(label=label),vjust="inward",hjust="inward", na.rm=TRUE)+
    scale_x_discrete(breaks = breaks) +
    labs(x="mode name") +
    facet_wrap( ~ mode,scales = "free")
  pl
}

#' Plot leverage of emission wavelengths, excitation wavelengths and samples.
#'
#' @description Plot is interactive where you can select values with your mouse. A list of vectors is returned to remove this outliers in a further step from your samples. The labels to be shown can be selected by adding the quatile of samples with highest leverages to be labeled.
#'
#' @param cpl leverage, outpout from \code{\link[staRdom]{eempf_leverage}}
#' @param qlabel optional, quantile of which labels are shown (1 = all, 0 = no labels)
#'
#' @return list of three vectors containing the names of selected samples
#'
#' @seealso \code{\link[staRdom]{eempf_leverage_plot}}
#'
#' @export
#'
#' @importFrom graphics plot
#' @importFrom graphics text
#' @importFrom graphics identify
#' @import dplyr
#' @importFrom stats setNames
#'
#' @examples
#' data(pfres_comps1)
#'
#' leverage <- eempf_leverage(pfres_comps[[2]])
#' outliers <- eempf_leverage_ident(leverage)
eempf_leverage_ident <- function(cpl,qlabel=0.1){
  pl <- eempf_leverage_data(cpl,qlabel=qlabel) %>%
    mutate(label = ifelse(is.na(label),"",label))
  exclude <- lapply(pl$mode %>% unique(),function(mod){
    #mod="sample"
    data <- pl %>% filter(mode == mod)
    plot(data$x %>% factor(),data$leverage,xlab=mod,ylab="leverage")
    text(data$x %>% factor(),data$leverage,label=data$label)
    ide <- identify(data$x %>% factor(),data$leverage,label=data$label)
    ide <- data$x %>% factor() %>% .[ide] %>% as.character()
  }) %>%
    setNames(pl$mode %>% unique())
}

#' Plot components from a PARAFAC model
#'
#' @description Additionally a bar plot with the amounts of each component in each sample is produced.
#'
#' @param cp_out object of class parafac
#'
#' @return ggplot
#' @export
#'
#' @seealso \code{\link[staRdom]{ggeem}}, \code{\link[staRdom]{eempf_load_plot}}
#'
#' @examples
#' data(pfres_comps1)
#'
#' eempf_comp_load_plot(pfres_comps[[2]])
eempf_comp_load_plot <- function(cp_out){
  pl1 <- ggeem(cp_out)
  pl2 <- eempf_load_plot(cp_out)
  #pl1 %>% print()
  #pl2 %>% print()
  list(pl1,pl2)
}


#' Plot amount of each component in each sample as bar plot
#'
#' @param cp_out parafac model
#'
#' @return ggplot
#' @export
#'
#' @import dplyr
#' @importFrom tibble rownames_to_column
#'
#' @examples
#' data(pfres_comps1)
#'
#' eempf_load_plot(pfres_comps[[2]])
eempf_load_plot <- function(cp_out){
  cp_out <- norm2A(cp_out)
  (cp_out$A) %>%
    data.frame() %>%
    rownames_to_column("sample") %>%
    #mutate(sample = names[[3]]) %>%
    gather(comp,amount,-sample) %>%
    ggplot()+
    geom_bar(aes(x=sample,y=amount,fill=comp),stat="identity",width=0.8)+
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
}

#' 3D plots of PARAFAC components
#'
#' @description Interactive 3D plots are created using plotly.
#'
#' @param cp_out object of class parafac
#' @param which optional, if numeric selects certain component
#'
#' @return plotly plot
#' @export
#'
#' @importFrom plotly plot_ly
#' @importFrom plotly layout
#' @importFrom plotly add_surface
#' @importFrom tibble column_to_rownames
#' @importFrom tibble remove_rownames
#' @importFrom grDevices rainbow
#'
#' @examples
#' \dontrun{
#' data(pfres_comps1)
#'
#' eempf_comps3D(pfres_comps[[3]])
#' }
eempf_comps3D <- function(cp_out,which=NULL){
  data <- cp_out %>% eempf_comp_mat()
  z <- lapply(data,function(mat){
    #mat <- data[[1]]
    mat %>%
      data.frame() %>%
      #mutate_all(as.numeric()) %>%
      spread(em,value) %>%
      remove_rownames() %>%
      column_to_rownames("ex") %>%
      #select(-ex) %>%
      as.matrix()
  })
  ex <- lapply(data,function(mat){
    mat$ex %>% unique() %>% as.numeric() %>% as.vector()
  })
  em <- lapply(data,function(mat){
    mat$em %>% unique() %>% as.numeric() %>% as.vector()
  })
  scene <- list(xaxis=list(title="em"),
                yaxis=list(title="ex"))
  lapply(1:length(ex),function(comp){
    if(is.null(which) | comp %in% which){
      plot_ly(x=em[[comp]],y=ex[[comp]], z=z[[comp]],colors = rainbow(12)[9:1]) %>%
        layout(scene=scene) %>%
        add_surface()
    }
  })
}


#' Plot correlations of components in samples
#'
#' @description A pair plot showing correlations between samples is created.
#'
#' @param cp_out object of class parafac
#' @param lower style of lower plots, see \code{\link[GGally]{ggpairs}}
#' @param mapping aesthetic mapping, see \code{\link[GGally]{ggpairs}}
#' @param ... passed on to \code{\link[GGally]{ggpairs}}
#'
#' @return object of class ggplot
#' @export
#'
#' @seealso \code{\link[GGally]{ggpairs}}
#' @importFrom GGally ggpairs
#'
#' @examples
#' \donttest{
#' data(pfres_comps1)
#' eempf_corplot(pfres_comps[[1]])
#' }
#'
eempf_corplot <- function(cp_out,lower=list(continuous="smooth"),mapping=aes(alpha=0.2),...){
  cp_out %>%
    norm2A() %>%
    .$A %>%
    data.frame() %>%
    ggpairs(lower=lower,mapping=mapping,...)
}

#' Plot samples by means of whole sample, each single component and residuum
#'
#' @description  A raster of plots is created. Each column shows one sample. The top n rows show the n components from the model according their occurance in the certain samples. The second last row shows the residual, not covered by any component in the model and the last row shows the whole sample.
#'
#'
#' @param cp_out object of class parafac containing the generated model
#' @param eem_list object of class eemlist with all the samples that should be plotted
#' @param res_data optional, data of sample residuals related to the model, output from \code{\link[staRdom]{eempf_residuals}}
#' @param spp optional, samples per plot
#' @param select optional, character vector of samples you want to plot
#' @param residuals_only plot only residuals
#' @param cores number of cores to use for parallel processing
#'
#' @details eem_list may contain samples not used for modelling. Calculation is done by \code{\link[staRdom]{A_missing}}. This especially interesting if outliers are excluded prior modelling and should be evaluated again afterwards.
#'
#' @return several ggplot objects
#' @export
#'
#' @import dplyr
#' @import tidyr
#' @import ggplot2
#' @import parallel
#' @importFrom grDevices rainbow
#'
#' @examples
#' \donttest{
#' data(eem_list)
#' data(pfres_comps1)
#'
#' eempf_residuals_plot(pfres_comps[[3]],eem_list)
#' }
#'
eempf_residuals_plot <- function(cp_out,eem_list,res_data = NULL, spp = 5, select=NULL, residuals_only = FALSE , cores = parallel::detectCores(logical = FALSE)/2){
  #cp_out,eem_list,select=eem_names(eem_list)[10:19]
  if(is.null(res_data)){
    res_data <- eempf_residuals(cp_out,eem_list,select=select,cores = cores)
  }
  if (!is.null(select)){
    res_data <- res_data %>% filter(Sample %in% select)
  }
  breaksx <- res_data %>% select(ex) %>% unique() %>% filter(as.numeric(ex)%%50 == 0) %>% unlist()
  breaksy <- res_data %>% select(em) %>% unique() %>% filter(as.numeric(em)%%50 == 0) %>% unlist()
  #if(!is.numeric(fill_max)){
  if(residuals_only){
    res_data <- res_data %>%
      filter(type == "residual")
  }
  fill_max <- res_data$value %>% max(na.rm=TRUE)
  #}
  vals <- c(res_data$value %>% min(na.rm=TRUE),seq(from=0,to=fill_max,length.out = 50))
  vals <- (vals - min(vals))/diff(range(vals))
  ppp <- res_data$Sample %>% unique() %>% length() /spp
  ov_plot <- lapply(1:ceiling(ppp),function(pos){
    #pos <- 1
    pl <- res_data %>%
      filter(Sample %in% (res_data$Sample %>% unique() %>% .[(spp*(pos-1)+1):(spp*pos)])) %>%
      ggplot()+
      geom_raster(aes(x=ex,y=em,fill=value))+ #,interpolate=TRUE
      scale_fill_gradientn(colours=c(rainbow(70)[62],rainbow(70)[50:1]),values=vals,limits = c(res_data$value %>% min(na.rm=TRUE),fill_max))+
      scale_x_discrete(breaks = breaksx) +
      scale_y_discrete(breaks = breaksy) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1))
    if(residuals_only){
      pl <- pl +
        facet_wrap(~ Sample)
    } else {
      pl <- pl +
        facet_grid(type ~ Sample)
    }
    pl
  })
  ov_plot
}

#' Plot results from a splithalf analysis
#'
#' @description Graphs of all components of all models are plotted to be compared.
#'
#' @param fits list of components data
#'
#' @return ggplot
#' @export
#'
#' @seealso \code{\link[staRdom]{splithalf}}
#'
#' @examples
#' data(sh)
#'
#' splithalf_plot(sh)
splithalf_plot <- function(fits){
  sel <- 0
  table <- lapply(fits,function(fit){
    sel <<- sel + 1
    c <- fit %>% lapply(eempf_comp_mat)
    tab <- lapply(c,function(c1){
      nc1 <- length(c1)
      nc2 <- 0
      lapply(c1,function(c2){
        nc2 <<- nc2 + 1
        c2 <- c2 %>%
          mutate(comps = nc1, comp = paste0("Comp.",nc2))
      }) %>%
        bind_rows()
    }) %>%
      bind_rows() %>%
      mutate(selection = sel) %>%
      group_by(comps,comp) %>%
      mutate(max_pos = which.max(value), max_em = em[max_pos],max_ex = ex[max_pos]) %>%
      mutate(exn = ifelse(em == max_em,ex,NA),emn = ifelse(ex == max_ex,em,NA)) %>%
      filter(!is.na(emn) | !is.na(exn)) %>%
      mutate(ex=exn,em=emn) %>%
      select(-exn,-emn,-max_pos,-max_em,-max_ex) %>%
      ungroup()
  }) %>%
    bind_rows()
  breaks <- table %>% select(ex,em) %>% unlist() %>% na.omit() %>% unique() %>% as.numeric() %>% .[.%%50 == 0]
  pl1 <- table %>%
    mutate(selection = factor(selection,ordered=FALSE)) %>%
    ggplot()+
    geom_line(data = . %>% filter(!is.na(ex)),aes(x=ex,y=value,colour=selection,group=selection),linetype=2)+
    geom_line(data = . %>% filter(!is.na(em)),aes(x=em,y=value,colour=selection,group=selection),linetype=1)+
    facet_grid(. ~ comp)+
    scale_x_discrete(breaks = breaks) +
    labs(x="wavelength") +
    theme(legend.position="none")
  pl1 %>% print()
}

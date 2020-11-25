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
#' @param ... arguments passed on to \code{\link[staRdom]{eempf_fits}} and \code{\link[staRdom]{eempf_plot_comps}}
#'
#' @return 3 objects of class ggplot
#' @export
#'
#' @seealso \code{\link[staRdom]{eempf_fits}}, \code{\link[staRdom]{eempf_plot_comps}}
#'
#' @examples
#' \donttest{
#' data(pf_models)
#'
#' eempf_compare(pf4)
#' }
eempf_compare <- function(pfres,...){
  #pfres <- pf4
  p1 <- eempf_fits(pfres,...)
  p2 <- eempf_plot_comps(pfres,type=1,...)
  p3 <- eempf_plot_comps(pfres,type=2,...)
  p1 %>% print()
  p2 %>% print()
  p3 %>% print()
  return(list(p1,p2,p3)) %>% invisible()
}

#' Fits vs. components of PARAFAC models are plotted
#'
#' @param pfres list of objects of class parafac
#' @param ... arguments passed on to ggplot
#'
#' @return object of class ggplot
#' @export
#'
#' @import ggplot2
#'
#' @examples
#' data(pf_models)
#'
#' eempf_fits(pf4)
eempf_fits <- function(pfres,...){
  #pfres <- pf1
  if(is.null(names(pfres))) names <- paste0("model",1:length(pfres)) else names <- names(pfres) #paste0("model",1:length(pfres))
  pl <- data.frame(comps=lapply(pfres,"[[","A") %>% lapply(ncol) %>% unlist(),
                   fit=lapply(pfres,"[[","Rsq") %>% unlist(), mod_name = names) %>%
    rowwise() %>%
    mutate(comps = ifelse(is.na(mod_name),paste0(comps, " comps"),paste0(mod_name," (",comps, " comps)"))) %>%
    ggplot(aes(x=comps,y=fit),...)+
    labs(x="model", y="model fit (Rsq)")+
    geom_point(size=5,shape=4)
  pl
}

#' Plot all components of PARAFAC models
#'
#' @description The components can be plottet in two ways: either as a colour map or as two lines (emission, excitation wavelengths) intersecting at the component maximum. If the list of provided models is named, these names are shown in the plot. Otherwise, the models are automatically named by "model#".
#'
#' @param pfres list of PARAFAC models
#' @param type 1 for a colour map and 2 for em and ex wavelength loadings
#' @param names logical, whether names of components should be written into the plot
#' @param contour in case of 3 dimensional component plots, contours are added
#' @param colpal "default" to use a subset of the rainbow palette or any custom vector of colors. A gradient will be produced from this vector. Larger vectors (e.g. 50 elements) can produce smoother gradients.
#' @param ... arguments passed on to other functions, e.g.
#'
#' @return object of class ggplot
#' @export
#'
#' @import dplyr
#' @import ggplot2
#' @importFrom grDevices rainbow
#'
#' @examples
#' data(pf_models)
#'
#' eempf_plot_comps(pf4, type = 1)
#'
#' # use a different colour scheme:
#' # eempf_plot_comps(pf4, type = 1, colpal = heat.colors(50))
#' eempf_plot_comps(pf4, type = 2)
#' eempf_plot_comps(list(pf4[[1]],pf4[[1]]), type=1)
#'
eempf_plot_comps <- function(pfres, type = 1, names = TRUE, contour = FALSE, colpal = "default", ...){
  #pfres <- pf3
  if(colpal[1] == "default"){
    colpal <- rainbow(75)[53:1]
    # warning("using rainbow colour palette")
  } else if(!is.vector(colpal)) {
    stop("Please provide a vector of colours!")
  }
  c <- pfres %>% lapply(eempf_comp_mat)
  if(is.null(names(c))) names(c) <- paste0("model",seq(1:length(c)))
  tab <- lapply(1:length(c),function(n){
    c1 <- c[[n]]
    mod_name <- names(c)[n]
    nc1 <- length(c1)
    nc2 <- 0
    lapply(c1,function(c2){
      nc2 <<- nc2 + 1
      c2 <- c2 %>%
        mutate(comps = nc1, comp = paste0("Comp.",nc2), modname = mod_name)
    }) %>%
      bind_rows()
  }) %>%
    bind_rows() %>%
    mutate(modname = factor(modname,levels = names(c))) %>%
    mutate_at(vars(ex,em,value),as.numeric)
  fill_max <- tab$value %>% max(na.rm=TRUE)
  vals <- seq(from = 0, to = fill_max,length.out = length(colpal))
  vals <- (vals - min(vals))/diff(range(vals))
  if(type==2){
    plot <- tab %>%
      group_by(modname,comp) %>%
      mutate(max_pos = which.max(value), max_em = em[max_pos],max_ex = ex[max_pos]) %>%
      mutate(exn = ifelse(em == max_em,ex,NA),emn = ifelse(ex == max_ex,em,NA)) %>%
      filter(!is.na(emn) | !is.na(exn)) %>%
      ungroup() %>%
      ggplot()+
      geom_line(aes(x=exn,y=value),colour="lightblue",group="excitation", na.rm=TRUE)+
      geom_line(aes(x=emn,y=value),colour="darkblue",group="emission", na.rm=TRUE)+
      theme(axis.text.x = element_text(angle = 90, hjust = 1))+
      facet_grid(comp ~ modname)+
      labs(x="wavelength (nm)")
  } else {
    diffs <- tab %>%
      select(-value,-comps) %>%
      gather("spec","wl", -comp, -modname) %>%
      group_by(comp,modname,spec) %>%
      unique() %>%
      summarise(slits = diff(wl) %>% n_distinct()) %>% #View()
      .$slits != 1

    plot <- tab %>%
      ggplot(aes(x = ex, y = em, z = value))

    if(any(diffs)){
      plot <- plot +
        #geom_raster(aes(fill = value), interpolate = interpolate)
        layer(mapping = aes(colour = value, fill = value),
              geom = "tile", stat = "identity", position = "identity")
    } else {
      plot <- plot +
        layer(mapping = aes(fill = value),
              geom = "raster", stat = "identity", position = "identity")
    }

    plot <- plot +
      scale_fill_gradientn(colours = colpal, values = vals, limits = c(tab$value %>% min(na.rm=TRUE),fill_max), aesthetics = c("fill", "colour"))+
      theme(axis.text.x = element_text(angle = 90, hjust = 1))+
      labs(x = "Excitation (nm)", y = "Emission (nm)") +
      facet_grid(comp ~ modname)
    if(contour){
      plot <- plot +
        geom_contour(colour = "black", size = 0.3)
    }
  }
  #print(pl)
  plot
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
#' data(pf_models)
#'
#' leverage <- eempf_leverage(pf4[[1]])
#' eempf_leverage_plot(leverage)
eempf_leverage_plot <- function(cpl, qlabel = 0.1){
  cpl <- eempf_leverage_data(cpl,qlabel=qlabel) %>%
    ungroup() %>%
    mutate(mode = str_replace(mode,"em","Emission (nm)") %>% str_replace("ex","Excitation (nm)") %>% str_replace("sample","Sample"))
  breaks <- seq(0,1000,50)
  breaks2 <- cpl$x[cpl$mode == "Emission (nm)"] %>%
    as.numeric() %>%
    range(na.rm = TRUE)
  breaks2 <- breaks[breaks >= breaks2[1] & breaks <= breaks2[2]]
  vals2 <- cpl$x[cpl$mode == "Emission (nm)"][findInterval(breaks2, cpl$x[cpl$mode =="Emission (nm)"])]
  breaks1 <- cpl$x[cpl$mode == "Excitation (nm)"] %>%
    as.numeric() %>%
    range(na.rm = TRUE)
  breaks1 <- breaks[breaks >= breaks1[1] & breaks <= breaks1[2]]
  vals <- cpl$x[cpl$mode == "Excitation (nm)"][findInterval(breaks1, cpl$x[cpl$mode =="Excitation (nm)"])]
  pl <- cpl %>%
    ggplot(aes(x=x,y=leverage))+
    geom_point(alpha = 0.4)+
    geom_text(aes(label=label),vjust="inward",hjust="inward", na.rm=TRUE, check_overlap = TRUE)+
    scale_x_discrete(labels = c(breaks1,breaks2), breaks = c(vals, vals2)) +
    labs(x="Variables (wavelengths or samples)", y = "Leverage") +
    facet_wrap( ~ mode, scales = "free")
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
#' data(pf_models)
#'
#' leverage <- eempf_leverage(pf4[[1]])
#' outliers <- eempf_leverage_ident(leverage)
eempf_leverage_ident <- function(cpl,qlabel=0.1){
  pl <- eempf_leverage_data(cpl,qlabel=qlabel) %>%
    mutate(label = ifelse(is.na(label),"",label))
  exclude <- lapply(pl$mode %>% unique(),function(mod){
    data <- pl %>% filter(mode == mod)
    plot(data$x %>% factor(), data$leverage, xlab = mod, ylab = "leverage")
    text(data$x %>% factor(), data$leverage, label = data$label)
    ide <- identify(data$x %>% factor(), data$leverage) #, labels = data$label
    ide <- data$x %>% factor() %>% .[ide] %>% as.character()
  }) %>%
    setNames(pl$mode %>% unique())
}

#' Plot components from a PARAFAC model
#'
#' @description Additionally a bar plot with the amounts of each component in each sample is produced.
#'
#' @param pfmodel object of class parafac
#' @param ... attributes passe don to \code{\link[staRdom]{ggeem}}
#'
#' @return ggplot
#' @export
#'
#' @seealso \code{\link[staRdom]{ggeem}}, \code{\link[staRdom]{eempf_load_plot}}
#'
#' @examples
#' data(pf_models)
#'
#' eempf_comp_load_plot(pf4[[1]])
eempf_comp_load_plot <- function(pfmodel,...){
  pl1 <- ggeem(pfmodel,...)
  pl2 <- eempf_load_plot(pfmodel)
  #pl1 %>% print()
  #pl2 %>% print()
  list(pl1,pl2)
}


#' Plot amount of each component in each sample as bar plot
#'
#' @param pfmodel parafac model
#'
#' @return ggplot
#' @export
#'
#' @import dplyr
#' @importFrom tibble rownames_to_column
#'
#' @examples
#' data(pf_models)
#'
#' eempf_load_plot(pf4[[1]])
eempf_load_plot <- function(pfmodel){
  pfmodel <- norm2A(pfmodel)
  (pfmodel$A) %>%
    data.frame() %>%
    rownames_to_column("sample") %>%
    #mutate(sample = names[[3]]) %>%
    gather(comp, amount, -sample) %>%
    ggplot()+
    geom_bar(aes(x = sample, y = amount, fill = comp), stat = "identity", width = 0.8)+
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
}

#' 3D plots of PARAFAC components
#'
#' @description Interactive 3D plots are created using plotly.
#'
#' @param pfmodel object of class parafac
#' @param which optional, if numeric selects certain component
#'
#' @return plotly plot
#' @export
#'
#' @importFrom tibble column_to_rownames
#' @importFrom tibble remove_rownames
#' @importFrom grDevices rainbow
#'
#' @examples
#' \dontrun{
#' data(pf_models)
#'
#' eempf_comps3D(pf4[[1]])
#' }
eempf_comps3D <- function(pfmodel, which = NULL){
  data <- pfmodel %>% eempf_comp_mat()
  z <- lapply(data,function(mat){
    mat %>%
      data.frame() %>%
      spread(em,value) %>%
      remove_rownames() %>%
      column_to_rownames("ex") %>%
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
      plotly::plot_ly(x=em[[comp]],y=ex[[comp]], z=z[[comp]],colors = rainbow(12)[9:1]) %>%
        plotly::layout(scene=scene) %>%
        plotly::add_surface()
    }
  })
}

#' Plot correlations of components in samples
#'
#' @description A pair plot showing correlations between samples is created.
#'
#' @param pfmodel object of class parafac
#' @param normalisation logical, whether normalisation is undone or not
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
#' data(pf_models)
#' eempf_corplot(pf4[[1]])
#' }
#'
eempf_corplot <- function(pfmodel,normalisation=FALSE,lower=list(continuous="smooth"),mapping=aes(alpha=0.2),...){
  if(normalisation) pfmodel <- norm2A(pfmodel)
  pfmodel %>%
    .$A %>%
    data.frame() %>%
    ggpairs(lower=lower,mapping=mapping,...)
}

#' Plot samples by means of whole sample, each single component and residuum
#'
#' @description  A raster of plots is created. Each column shows one sample. The top n rows show the n components from the model according their occurance in the certain samples. The second last row shows the residual, not covered by any component in the model and the last row shows the whole sample.
#'
#'
#' @param pfmodel object of class parafac containing the generated model
#' @param eem_list object of class eemlist with all the samples that should be plotted
#' @param res_data optional, data of sample residuals related to the model, output from \code{\link[staRdom]{eempf_residuals}}
#' @param spp optional, samples per plot
#' @param select optional, character vector of samples you want to plot
#' @param residuals_only plot only residuals
#' @param cores number of cores to use for parallel processing
#' @param contour logical, states whether contours should be plotted
#' @param colpal "default" to use a subset of the rainbow palette or any custom vector of colors. A gradient will be produced from this vector. Larger vectors (e.g. 50 elements) can produce smoother gradients.
#'
#' @details eem_list may contain samples not used for modelling. Calculation is done by \code{\link[staRdom]{A_missing}}. This especially interesting if outliers are excluded prior modelling and should be evaluated again afterwards.
#' Usually, residuals contain negative values, while these is the exception in samples and PARAFAC components. Therefore, we decided to use a similar colour palette as in the other plot functions but adding a purple tone for negative values.
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
#' data(pf_models)
#'
#' eem_list <- eem_extract(eem_list, 1:10)
#'
#' eem_list <- eem_rem_scat(eem_list, rep(TRUE, 4), c(15,10,16,12))
#'
#' eempf_residuals_plot(pf4[[1]], eem_list, cores = 2)
#'
#' # use other colour schemes:
#' # eempf_residuals_plot(pf4[[1]], eem_list, colpal = c("blue",heat.colors(50)))
#' # plots <- eempf_residuals_plot(pf4[[1]], eem_list)
#' # lapply(plots, function(pl){
#' #   pl +
#' #     scale_fill_viridis_c() +
#' #     scale_colour_viridis_c()
#' # })
#'
#' }
#'
eempf_residuals_plot <- function(pfmodel, eem_list, res_data = NULL, spp = 5, select = NULL, residuals_only = FALSE , cores = parallel::detectCores(logical = FALSE), contour = FALSE, colpal = "default"){
  if(is.null(res_data)){
    res_data <- eempf_residuals(pfmodel,eem_list,select=select,cores = cores)
  }
  if (!is.null(select)){
    res_data <- res_data %>% filter(Sample %in% select)
  }
  res_data <- res_data %>%
    mutate_at(vars(ex,em,value),as.numeric)
  if(residuals_only){
    res_data <- res_data %>%
      filter(type == "residual")
  }
  if(colpal[1] == "default"){
    colpal <- c(rainbow(70)[62], rainbow(70)[50:1])
    # warning("using rainbow colour palette")
  } else if(!is.vector(colpal)) {
    stop("Please provide a vector of colours!")
  }
  fill_max <- res_data$value %>% max(na.rm=TRUE)
  vals <- c(res_data$value %>% min(na.rm = TRUE), seq(from = 0, to = fill_max, length.out = length(colpal) - 1))
  vals <- (vals - min(vals))/diff(range(vals))
  ppp <- res_data$Sample %>% unique() %>% length() /spp
  ov_plot <- lapply(1:ceiling(ppp),function(pos){
    pl <- res_data %>%
      filter(Sample %in% (res_data$Sample %>% unique() %>% .[(spp*(pos-1)+1):(spp*pos)])) %>%
      ggplot(aes(x=ex,y=em,z=value))

    diffs <- res_data %>%
      select(-value, -type) %>%
      gather("spec","wl", -Sample) %>%
      group_by(Sample,spec) %>%
      unique() %>%
      #arrange(sample,spec,wl) %>%
      #mutate(diffs = wl - lag(wl))
      summarise(slits = diff(wl) %>% n_distinct()) %>%
      .$slits != 1

    if(any(diffs)){
      pl <- pl +
        #geom_raster(aes(fill = value), interpolate = interpolate)
        layer(mapping = aes(colour = value, fill = value),
              geom = "tile", stat = "identity", position = "identity") +
        scale_colour_gradientn(colours = colpal, values = vals, limits = c(res_data$value %>% min(na.rm = TRUE), fill_max))
    } else {
      pl <- pl +
        layer(mapping = aes(fill = value),
              geom = "raster", stat = "identity", position = "identity")
    }
    pl <- pl +
      #geom_raster(aes(fill=value))+ #,interpolate=TRUE
      scale_fill_gradientn(colours = colpal, values = vals, limits = c(res_data$value %>% min(na.rm = TRUE), fill_max))+
      theme(axis.text.x = element_text(angle = 90, hjust = 1))+
      labs(x = "Excitation (nm)", y = "Emission (nm)", fill = "fluorescence")
    if(contour){
      pl <- pl +
        geom_contour(colour = "black", size = 0.3)
    }
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
#' str(sh)
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
      ungroup() %>%
      mutate_at(vars(em,ex,value),as.numeric)
  }) %>%
    bind_rows()
  pl1 <- table %>%
    mutate(selection = factor(selection,ordered=FALSE)) %>%
    ggplot()+
    geom_line(data = . %>% filter(!is.na(ex)),aes(x = ex, y = value, colour = selection, group = selection), linetype = 2)+
    geom_line(data = . %>% filter(!is.na(em)),aes(x = em, y = value, colour = selection, group = selection), linetype = 1)+
    labs(x="Wavelength (nm)",y="Loading") +
    theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1))+
    facet_grid(. ~ comp)
  pl1 %>% print()
}


#' Create a html report of a PARAFAC analysis
#'
#' @param pfmodel PARAFAC model
#' @param export path to exported html file
#' @param eem_list optional EEM data
#' @param absorbance optional absorbance data
#' @param meta optional meta data table
#' @param metacolumns optional column names of metadata table
#' @param splithalf optional logical, states whether split-half analysis should be included
#' @param shmodel optional results from split-half analysis. If this data is not supplied but EEM data is available the split-half analysis is calculated on the creation of the report. Calculating the split-half analysis takes some time!
#' @param performance calculating model performance: \code{\link[staRdom]{eempf_eemqual}}
#' @param residuals logical, whether residuals are plotted in the report
#' @param spp plots per page for loadgins and residuals plot
#' @param cores cores to be used for the calculation
#' @param ... arguments to or from other functions
#'
#' @return TRUE if report was created
#' @export
#'
#' @examples
#' \donttest{
#' folder <- system.file("extdata/EEMs", package = "staRdom") # load example data
#' eem_list <- eem_read(folder, recursive = TRUE, import_function = eem_csv)
#'
#' abs_folder <- system.file("extdata/absorbance", package = "staRdom") # load example data
#' absorbance <- absorbance_read(abs_folder, cores = 2)
#'
#' metatable <- system.file("extdata/metatable_dreem.csv",package = "staRdom")
#' meta <- read.table(metatable, header = TRUE, sep = ",", dec = ".", row.names = 1)
#'
#' checked <- eem_checkdata(eem_list, absorbance, metadata = meta,
#' metacolumns = "dilution", error = FALSE)
#'
#' eem_names(eem_list)
#' pfm <- A_missing(eem_list,pf4[[1]], cores = 2)
#' eempf_report(pfm, export = "~/pf_report.html", eem_list = eem_list,
#'              absorbance = absorbance, meta = metatable, metacolumns = "dilution", cores = 2)
#'
#' }
eempf_report <- function(pfmodel, export, eem_list = NULL, absorbance = NULL, meta = NULL, metacolumns = NULL, splithalf = FALSE, shmodel = NULL, performance = FALSE, residuals = FALSE, spp = 5, cores = parallel::detectCores(logical=FALSE), ...){
  #shmodel <- sh
  # rm(shmodel)
  #splithalf = TRUE
  rmdfile <- system.file("PARAFAC_report.Rmd",package="staRdom")
  dir <- dirname(export)
  file <- basename(export)
  if(splithalf | performance){
    if(!is.null(eem_list) & is.null(shmodel)){
      splithalf <- splithalf(eem_list, comps = ncol(pfmodel$A), normalise = !is.null(attr(pfmodel,"norm_factors")),...)
      tcc <- splithalf_tcc(shmodel)
    } else if (!is.null(shmodel)){
      splithalf <- shmodel
      tcc <- splithalf_tcc(shmodel)
    } else {
      tcc <- NULL
      warning("Split-half analysis and/or model performance could not be incorporated due to missing EEM data or an already calculated split-half analysis.",fill=TRUE)
    }
  }
  if(performance){
    if(!is.null(eem_list)){
      performance <- eempf_eemqual(pfmodel,eem_list,splithalf)
    } else{
      warning("For a performance calculation, EEM data is needed!")
    }
  }
  imgwidth <- nrow(pfmodel$A)/8
  rmarkdown::render(rmdfile, output_file = file, output_dir = dir, params = list(pfmodel = pfmodel, eem_list = eem_list, absorbance = absorbance, meta = meta, tcc = tcc, metacolumns = metacolumns, splithalf = splithalf, performance = performance, residuals = residuals, spp = spp, imgwidth = imgwidth, cores = cores))
  TRUE
}

#' Plot results from an SSC check
#'
#' @param ssccheck outpout from \code{\link{eempf_ssccheck}}
#'
#' @return ggplot element
#' @export
#'
#' @import dplyr
#'
#' @examples
#' \donttest{
#' data(pf_models)
#'
#' ssccheck <- eempf_ssccheck(pfmodels = pf3[1:3], cores = 2)
#' eempf_plot_ssccheck(ssccheck)
#' }
eempf_plot_ssccheck <- function(ssccheck){
  ssccheck %>%
    mutate(excitation = B,  emission = C) %>%
    select(-B,-C) %>%
    gather("spectrum", "TCC", excitation, emission) %>%
    group_by(comp, spectrum) %>%
    mutate(mean = mean(TCC), min = min(TCC), max = max(TCC)) %>%
    ggplot()+
    geom_point(aes(x=comp + 0.1 * ifelse(spectrum == "excitation", 1, -1), y = mean, colour = comp, group = comp), shape = 21, size = 3)+
    geom_point(aes(x=comp + 0.1 * ifelse(spectrum == "excitation", 1, -1), y = TCC, colour = comp, group = comp), alpha = 0.4)+
    geom_errorbar(aes(x=comp + 0.1 * ifelse(spectrum == "excitation", 1, -1), ymin = min, ymax = max, colour = comp, group = comp, linetype = spectrum))+
    scale_color_viridis_c(guide = FALSE)+
    labs(x = "Component", y = attr(ssccheck,"method"), linetype = "", shape = "")+
    scale_x_continuous(breaks = c(1:max(ssccheck$comp)))
}

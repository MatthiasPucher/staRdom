##
#' EEM spectra plotted with ggplot2
#'
#' @description Plots from EEM spectra of class \code{ggplot}. In case you work with a larger number of EEMs and want to show then in several plots, you can use \code{\link{eem_overview_plot}}.
#'
#' @param data eem, eemlist, parafac or data.frame. The details are given under 'Details'.
#' @param fill_max sets the maximum fluorescence value for the colour scale. This is mainly used by other functions, and makes different plots visually comparable.
#' @param colpal  "default" to use the viridis colour palette, "rainbow" to use a subset of the rainbow palette, any custom vector of colors or a colour palette. A gradient will be produced from this vector. Larger vectors (e.g. 50 elements) can produce smoother gradients.
#' @param contour logical, whether contours should be plotted (default FALSE), see \code{\link[ggplot2]{geom_contour}}
#' @param interpolate logical, whether fluorescence should be interpolated, see \code{\link[ggplot2]{geom_raster}}
#' @param redneg deprecated! logical, whether negative values should be coloured discreet.
#' @param eemlist_order logical, in case of an eemlist, the order of samples in the plot is the same as in the eemlist, alphabetically otherwise
#' @param ... parameters passed on to \code{\link[ggplot2]{ggplot}}.
#'
#' @details The data can be of different sources:
#'     eem: a single EEM spectrum is plotted
#'     eemlist: all spectra of the samples are plotted, arranged in a grid
#'     data.frame: a data.frame containing EEM data. Can be created by e.g. \code{as.data.frame.eem}
#'     parafac: a PARAFAC model, the components are plotted then.
#'
#'     Using redneg you can give negative values a reddish colour. This can help identifying these parts in samples or components. Negative values are physically not possible and can only be the result of measuring errors, model deviations and problems with interpolated values.
#'
#'      Interpolation (interpolate = TRUE) leeds to smoother plots. The default is FALSE because it might cover small scale inconsistencies.
#'
#'      Contours (contour = TRUE)can be added to the EEM plots.
#'
#'     A colour palette can be specified using the argument colpal.
#'
#'     Plotting distinct samples can be done using \code{\link[eemR]{eem_extract}}. Please see example.
#'
#' @return a ggplot object
#'
#' @export
#' @import ggplot2 dplyr tidyr eemR
#' @importFrom grDevices rainbow col2rgb
#'
#' @examples
#' ## plotting two distinct samples
#' data(eem_list)
#' eem_names(eem_list)
#' eem <- eem_extract(eem_list,c("^d667sf$", "^d661sf$"),keep=TRUE)
#' ggeem(eem)
#'
#' # the former redneg argument is deprecated, please see a similar looking example below!
#' #ggeem(eem, redneg = TRUE)
#' ggeem(eem, colpal = c(rainbow(75)[58],rainbow(75)[53:1]))
#'
#' # use any custom colour palette
#' ggeem(eem, colpal = heat.colors(50))
#' # needs package matlab to be installed:
#' # ggeem(eem, colpal = matlab::jet.colors(50))
#' # or by adding ggplot2 colour and fill functions:
#' # ggeem(eem)+
#' #   scale_fill_viridis_c()+
#' #   scale_color_viridis_c()
#'
#' ggeem(eem, interpolate = TRUE)
#' ggeem(eem, contour = TRUE)
ggeem <- function(data, fill_max=FALSE, ...) UseMethod("ggeem")

#' @rdname ggeem
#' @export
ggeem.default <- function(data, fill_max=FALSE, ...){
  stop("data is not of a suitable format!")
}

#' @rdname ggeem
#' @export
ggeem.eemlist <- function(data, fill_max = FALSE, eemlist_order = TRUE, ...)
{
  table <- data %>% lapply(as.data.frame) %>% bind_rows()
  if(eemlist_order) table$sample <- table$sample %>% factor(levels = table$sample %>% unique())
  #filename <- paste0('EEM_spectra_',suffix,format(Sys.time(), "%Y%m%d_%H%M%S"))
  ggeem(table,fill_max=fill_max,...)
}

#' @rdname ggeem
#' @export
ggeem.eem <- function(data, fill_max = FALSE, ...)
{
  table <- data %>% as.data.frame()
  #filename <- paste0('EEM_spectrum_',table$sample[1],"_",format(Sys.time(), "%Y%m%d_%H%M%S"))
  ggeem(table,fill_max=fill_max,...)
}

#' @rdname ggeem
#' @export
ggeem.parafac <- function(data, fill_max = FALSE, ...)
{
  table <- data %>% eempf_comp_mat() #eem_list
  table <- lapply(table %>% names(),function(name){
    table[[name]] %>% mutate(sample = name)
  }) %>% bind_rows() %>%
    mutate(sample = factor(sample, levels = colnames(data$A)))
  #filename <- paste0('EEM_PARAFAC_components_',suffix,format(Sys.time(), "%Y%m%d_%H%M%S"))
  ggeem(table,fill_max=fill_max,...)
}

#' @rdname ggeem
#' @export
ggeem.data.frame <- function(data, fill_max=FALSE, colpal = "default", contour = FALSE, interpolate = FALSE, redneg = NULL, ...)
{
  #data <- table
  if(!is.null(redneg)){
    warning("redneg is deprecated and will be ignored! Please use the argument 'colpal = c(rainbow(75)[58],rainbow(75)[51:1])' to produce similar behaviour.")
  }
  if(is.vector(colpal)){
    if(colpal[1] == "rainbow"){
      colpal <- rainbow(75)[53:1]
      # warning("using rainbow colour palette")
    } else if (colpal[1] == "default"){
      colpal <- viridisLite::viridis(50)
    }
  } else if(is.function(colpal) & class(try(col2rgb(colpal(1)),silent=TRUE))[1] != "try-error"){
    colpal <- colpal(50)
    # warning("using rainbow colour palette")
  } else {
    stop("Please provide a palette or a vector of colours as argument colpal!")
  }
  table <- data %>%
    mutate_at(vars(ex,em,value),as.numeric)
  if(!is.numeric(fill_max)){
    fill_max <- table$value %>% max(na.rm=TRUE)
  }
  #values_fill <- seq(0,fill_max,length.out = 55)
  #diff()
  diffs <- table %>%
    select(-value) %>%
    gather("spec","wl", -sample) %>%
    group_by(sample,spec) %>%
    unique() %>%
    #arrange(sample,spec,wl) %>%
    #mutate(diffs = wl - lag(wl))
    summarise(slits = diff(wl) %>% n_distinct()) %>%
    .$slits != 1

  plot <- table %>%
    ggplot(aes(x = ex, y = em, z = value))+
    labs(x = "Excitation (nm)", y = "Emission (nm)")

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
    #geom_tile(aes(fill = value, colour = value))+
    #scale_fill_gradient2(low="blue",mid="yellow",high="red")
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))+
    facet_wrap(~ sample)
  if(contour){
    plot <- plot +
      geom_contour(colour = "black", size = 0.3, ...)
  }
  if(table$value %>% min(na.rm=TRUE) < 0){
    vals <- c(table$value %>% min(na.rm = TRUE), seq(from = 0, to = fill_max, length.out = length(colpal) - 1))
    vals <- (vals - min(vals))/diff(range(vals))
    plot <- plot +
      scale_fill_gradientn(colours = colpal, values = vals, limits = c(table$value %>% min(na.rm=TRUE),fill_max))+
      scale_colour_gradientn(colours = colpal, values = vals, limits = c(table$value %>% min(na.rm=TRUE),fill_max))
  } else {
    plot <- plot +
      scale_fill_gradientn(colours = colpal, limits = c(0,fill_max))+
      scale_colour_gradientn(colours = colpal, limits = c(0,fill_max))
  }
  plot
}


#' Plot fluorescence data from several samples split into several plots.
#'
#' @param data fluorescence data of class eemlist
#' @param spp number of samples per plot or a vector with the numbers of rows and columns in the plot.
#' @param ... arguments passed on to \code{\link[staRdom]{ggeem}}
#'
#' @return list of ggplots
#' @export
#'
#' @examples
#' \donttest{
#' data(eem_list)
#' eem_overview_plot(eem_list, spp = 9)
#'
#' # define number of rows and columns in plot
#' eem_overview_plot(eem_list, spp = c(3, 5))
#' }
eem_overview_plot <- function(data, spp = 8,...){
  #data <- eem_list
  ppp <- data %>% length()/prod(spp)
  fill_max <- data %>% eem_scale_ext() %>% .[2]
  #print(fill_max)
  ov_plot <- lapply(1:ceiling(ppp),function(pos){
    #pos <- 1
    pl <- data[(prod(spp)*(pos-1)+1):(prod(spp)*pos)] %>%
      `attr<-`("class", "eemlist") %>%
      ggeem(fill_max=fill_max,...)#,...
    if(length(spp) == 2){
      pl <- pl +
        facet_wrap(sample ~ ., nrow = spp[1], ncol = spp[2])
    }
    pl
  })
  ov_plot
}

#' Mark EEM scatter bands
#'
#' @description geom_eemscatter draws dashed lines at the locations where scatter is expected in EEMs from water samples.
#'
#' @param scatter logical vector of size 4 stating the scatter bands to be marked. The order is Raman 1, Raman 2, Rayleigh 1, Rayleigh 2.
#' @param ... additional arguments to be passed on to \code{\link[ggplot2]{geom_function}}
#'
#' @returns a layer to a ggplot2
#' @export
#'
#' @import ggplot2
#'
#' @examples
#' require(tidyr)
#'
#' eem_list %>%
#'   eem_extract(eem_names((eem_list))[1], keep = TRUE) %>%
#'   ggeem() +
#'   geom_eemscatter() +
#'   coord_cartesian(xlim = c(250,455), ylim = c(290,578))
geom_eemscatter <- function(scatter = rep(TRUE, 4), ...){
  raman <- function(ex, order = 1){
    order * 1/(1/ex - 3.4e-4)
  }

  rayleigh <- function(ex, order = 1){
    order*ex
  }

  layers <- list()
  if (scatter[1]){
    layers <- c(layers,
                geom_function(fun = raman, mapping = aes(z = NULL), ...))
  }
  if (scatter[2]){
    layers <- c(layers,
                geom_function(fun = raman, mapping = aes(z = NULL), args = list(order = 2), ...))
  }
  if (scatter[3]){
    layers <- c(layers,
                geom_function(fun = rayleigh, mapping = aes(z = NULL), ...))
  }
  if (scatter[3]){
    layers <- c(layers,
                geom_function(fun = rayleigh, mapping = aes(z = NULL), args = list(order = 2), ...))
  }

  return(layers)
}


#' Mark common EEM peaks
#'
#' @description geom_eempeakloc marks the locations of commonly used peaks to EEMs of water samples
#'
#' @param data data.frame containing information about the peaks, such as names, locations and ranges
#' @param ... additional arguments to be passed on to \code{\link[ggplot2]{geom_function}}
#'
#' @details data, the data.frame with the peaks, is included in the package and accessible using data(peaks).
#' It is possible to alter that data.frame and e.g. provide altered peaks or a selection of peaks using common functions on data.frames.
#' geom_eempeakloc plots several labels which can be close together in smaller plots. By default, labels are omitted to keep it readable, but this means important information is lost.
#' If the ggrepel package is installed, it is used for plotting the labels and the result is much more convenient. You have to install ggrepel manually to use it with this function.
#'
#' @returns a layer to a ggplot2
#' @export
#'
#' @import ggplot2
#'
#' @examples
#' require(tidyr)
#'
#' eem_list %>%
#'   eem_extract(eem_names((eem_list))[1], keep = TRUE) %>%
#'   ggeem() +
#'   geom_eempeakloc()
#'
#' # We plot only component 3 of the PARAFAC model to avoid letter overload
#' pf4[[1]] %>%
#'   (\(x, comp = 3) { #please change comp to the index of the comp. to plot.
#'     x$A <- x$A[, comp, drop = FALSE]
#'     x$B <- x$B[, comp, drop = FALSE]
#'     x$C <- x$C[, comp, drop = FALSE]
#'     x
#'   })() %>%
#'   ggeem() +
#'   geom_eempeakloc() +
#'   coord_cartesian(xlim = c(250,455), ylim = c(290,578))
geom_eempeakloc <- function(data = peaks, ...){
  list(
    geom_point(mapping = aes(x = ex, y = em), data = peaks, ...),
    geom_errorbar(mapping = aes(x = ex, ymin = em_min, ymax = em_max), data = peaks, ...),
    geom_errorbar(mapping = aes(y = em, xmin = ex_min, xmax = ex_max), data = peaks, ...),
    if (requireNamespace("ggrepel", quietly = TRUE)) {
      ggrepel::geom_text_repel(mapping = aes(x = ex, y = em, label = name), data = peaks, ...)
    } else {
      geom_text(mapping = aes(x = ex, y = em, label = name), data = peaks, check_overlap = TRUE, ...)
    }
  )
}

#' Mark common reagions to EEMs to show molecular groups
#'
#' @description geom_eemregions marks the locations of commonly destinguished molecular groups in EEMs
#'
#' @param lim limit the extend of the lines separating the groups
#' @param detail logical, whether the reagions are given with numbers only or with more detailed information
#' @param ... additional arguments to be passed on to \code{\link[ggplot2]{geom_function}}, \code{\link[ggplot2]{geom_segment}} and \code{\link[ggplot2]{geom_text}}
#'
#' @details data, the data.frame with the peaks, is included in the package and accessible using data(peaks).
#' It is possible to alter that data.frame and e.g. provide altered peaks or a selection of peaks using common functions on data.frames.
#'
#' @returns a layer to a ggplot2
#' @export
#'
#' @import ggplot2
#'
#' @examples
#' require(tidyr)
#'
#' eem_list %>%
#'   eem_extract(eem_names((eem_list))[1], keep = TRUE) %>%
#'   ggeem() +
#'   geom_eemregions()
geom_eemregions <- function(lim = c(ex_min = 200, ex_max = 500, em_min = 250, em_max = 700), detail = TRUE, ...){

  label <- c("I", "II", "III", "IV", "V")

  if(detail){
    label2 <- c("Aromatic protein I", "Aromatic protein II", "Fulvic acid-like", "Soluble microbial\nby-product-like", "Humic acid-like")
    label <- paste(label, label2, sep = "\n")
  }



  text <- data.frame(ex = c(rep(mean(c(250, lim["ex_min"])),3), mean(c(rep(250, 2), min(lim["ex_max"],342))), min(320, lim["ex_max"])),
                     em = c(mean(c(lim["em_min"],330)), 355, mean(c(lim["em_max"],380)), mean(c(lim["em_min"],380)), mean(c(lim["em_max"],380))),
                     label = label,
                     value = rep(0,5))

  list(
    geom_segment(aes(y = max(lim["em_min"],260), yend = lim["em_max"], x = 250, xend = 250), linetype = 3, ...),
    geom_segment(aes(y = 330, yend = 330, x = lim["ex_min"], xend = 250), linetype = 3, ...),
    geom_segment(aes(y = 380, yend = 380, x = lim["ex_min"], xend = min(lim["ex_max"],342)), linetype = 3, ...),
    geom_function(fun = function(x) {y <- -67.016 + 1.308 * x}, mapping = aes(z = NULL), ...),
    geom_text(mapping = aes(x = ex, y = em, label = label), data = text, ...)
  )
}

#' Add layers of scatter bands, molecular regions and common peaks to EEM plots
#'
#' @description This function is a wrapper for \code{\link[staRdom]{geom_eemregions}}, \code{\link[staRdom]{geom_eempeakloc}}, \code{\link[staRdom]{geom_eemscatter}} and can limit the plot extend to the area of the original EEM.
#' Therefore, it is not added using '+' but the plot has to be supplied as an argument.
#'
#' @param plot the ggplot where the layers are added
#' @param scatter logical vector of size 4 stating the scatter bands to be marked. The order is Raman 1, Raman 2, Rayleigh 1, Rayleigh 2.
#' @param regions logical, whether molecular regions are marked
#' @param peakloc logical, whether common EEM peaks are marked
#' @param limit logical, whether the plot is limited to the original EEM
#' @param ... additional arguments passed to geom_eemregions(), geom_eemscatter() and geom_eempeakloc().
#'
#' @returns an altered ggplot2
#' @export
#'
#' @examples
#' require(tidyr)
#'
#' eem_list %>%
#'  eem_extract(eem_names((eem_list))[1], keep = TRUE) %>%
#'  ggeem() %>%
#'  ggeem_overlay()
ggeem_overlay <- function(plot, scatter = rep(TRUE, 4), regions = TRUE, peakloc = TRUE, limit = TRUE, ...){

  if(limit){
    built <- ggplot_build(plot)
    panel_ranges <- built$layout$panel_params[[1]]
    xlim <- panel_ranges$x.range
    ylim <- panel_ranges$y.range

    plot <- plot +
      coord_cartesian(xlim = xlim, ylim = ylim)
  }

  if(any(scatter)){
    plot <- plot +
      geom_eemscatter(scatter = scatter, linetype = 2, ...)
  }

  if(regions){
    plot <- plot +
      geom_eemregions(...)
  }

  if(peakloc){
    plot <- plot +
      geom_eempeakloc(...)
  }

  plot

}


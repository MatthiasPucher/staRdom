##
#' EEM spectra plotted with ggplot2
#'
#' @description Plots from EEM spectra of class \code{ggplot}. In case you work with a larger number of EEMs and want to show then in several plots, you can use \code{\link{eem_overview_plot}}.
#'
#' @param data eem, eemlist, parafac or data.frame. The details are given under 'Details'.
#' @param fill_max sets the maximum fluorescence value for the colour scale. This is mainly used by other functions, and makes different plots visually comparable.
#' @param colpal "default" to use a subset of the rainbow palette or any custom vector of colors. A gradient will be produced from this vector. Larger vectors (e.g. 50 elements) can produce smoother gradients.
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
#'     Plotting distinct samples can be done using \code{\link{eem_extract}}. Please see example.
#'
#' @return a ggplot object
#'
#' @export
#' @import ggplot2 dplyr tidyr eemR
#' @importFrom grDevices rainbow
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
  if(colpal[1] == "default"){
    colpal <- rainbow(75)[53:1]
    # warning("using rainbow colour palette")
  } else if(!is.vector(colpal)) {
    stop("Please provide a vector of colours!")
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

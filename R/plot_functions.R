##
#' EEM spectra plotted with ggplot2
#'
#' @description Plots from EEM spectra of class \code{ggplot}. In case you work with a larger number of EEMs and want to show then in several plots, you can use \code{\link{eem_overview_plot}}.
#'
#' @param data eem, eemlist, parafac or data.frame. The details are given under 'Details'.
#' @param fill_max sets the maximum fluorescence value for the colour scale. This is mainly used by other functions, and makes different plots visually comparable.
#' @param redneg logical, whether negative values should be coloured discreet.
#' @param ... parameters passed on to \code{ggplot}.
#'
#' @details The data can be of different sources:
#'     eem: a single EEM spectrum is plotted
#'     eemlist: all spectra of the samples are plotted, arranged in a grid
#'     data.frame: a data.frame containing EEM data. Can be created by e.g. \code{as.data.frame.eem}
#'     parafac: a PARAFAC model, the components are plotted then.
#'
#'     Using redneg you can give negative values a reddish colour. This can help identifying these parts in samples or components. Negative values are physically not possible and can only be the result of measuring errors, model deviations and problems with interpolated values.
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
#' eem <- eem_extract(eem_list,c("^dreem_667sf$", "^dreem_661sf$"),keep=TRUE)
#' ggeem(eem)
ggeem <- function(data, fill_max=FALSE, ...) UseMethod("ggeem")

#' @rdname ggeem
#' @export
ggeem.default <- function(data, fill_max=FALSE, ...){
  stop("data is not of a suitable format!")
}

#' @rdname ggeem
#' @export
ggeem.eemlist <- function(data,fill_max=FALSE,...)
{
  table <- data %>% lapply(as.data.frame) %>% bind_rows()
  #filename <- paste0('EEM_spectra_',suffix,format(Sys.time(), "%Y%m%d_%H%M%S"))
  ggeem(table,fill_max=fill_max,...)
}

#' @rdname ggeem
#' @export
ggeem.eem <- function(data,fill_max=FALSE,...)
{
  table <- data %>% as.data.frame()
  #filename <- paste0('EEM_spectrum_',table$sample[1],"_",format(Sys.time(), "%Y%m%d_%H%M%S"))
  ggeem(table,fill_max=fill_max,...)
}

#' @rdname ggeem
#' @export
ggeem.parafac <- function(data,fill_max=FALSE,...)
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
ggeem.data.frame <- function(data,fill_max=FALSE,redneg = FALSE, ...)
{
  if(!exists("colpal")){
    colpal <- rainbow
    # warning("using rainbow colour palette")
  }
  table <- data
  breaksx <- table %>% select(ex) %>% unique() %>% filter(as.numeric(ex)%%50 == 0) %>% unlist()
  breaksy <- table %>% select(em) %>% unique() %>% filter(as.numeric(em)%%50 == 0) %>% unlist()
  if(!is.numeric(fill_max)){
    fill_max <- table$value %>% max(na.rm=TRUE)
  }
  #values_fill <- seq(0,fill_max,length.out = 55)
  plot <- table %>% ggplot()+
    geom_raster(aes(x=ex,y=em,fill=value))+ #,interpolate=TRUE
    #scale_fill_gradient2(low="blue",mid="yellow",high="red")
    scale_x_discrete(breaks = breaksx) +
    scale_y_discrete(breaks = breaksy) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))+
    facet_wrap(~sample)
  if(table$value %>% min(na.rm=TRUE) < 0){
    vals <- c(table$value %>% min(na.rm=TRUE),seq(from=0,to=fill_max,length.out = 51))
    vals <- (vals - min(vals))/diff(range(vals))
    if(redneg){
      plot <- plot +
        scale_fill_gradientn(colours=c(colpal(75)[58],colpal(75)[51:1]),values=vals,limits = c(table$value %>% min(na.rm=TRUE),fill_max))
    } else {
      plot <- plot +
        scale_fill_gradientn(colours=colpal(75)[52:1],values=vals,limits = c(table$value %>% min(na.rm=TRUE),fill_max))
    }
  } else {
    plot <- plot +
      scale_fill_gradientn(colours=colpal(75)[51:1],limits = c(0,fill_max))
  }
  plot
}


#' Plot fluorescence data from several samples split into several plots.
#'
#' @param data fluorescence data of class eemlist
#' @param spp number of samples per plot
#' @param ... arguments passed on to \code{\link[staRdom]{ggeem}}
#'
#' @return list of ggplots
#' @export
#'
#' @examples
#' \donttest{
#' data(eem_list)
#' eem_overview_plot(eem_list,spp=9)
#' }
eem_overview_plot <- function(data,spp = 8,...){
  #data <- eem_list
  ppp <- data %>% length()/spp
  fill_max <- data %>% eem_scale_ext() %>% .[2]
  #print(fill_max)
  ov_plot <- lapply(1:ceiling(ppp),function(pos){
    #pos <- 1
    data[(spp*(pos-1)+1):(spp*pos)] %>%
      `attr<-`("class", "eemlist") %>%
      ggeem(fill_max=fill_max,...)#,...
  })
  ov_plot
}

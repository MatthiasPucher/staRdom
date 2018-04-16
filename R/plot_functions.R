##
#' EEM spectra plotted with ggplot2
#'
#' @description \code{ggeem} creates nice plots from EEM spectra of class ggplot. Plots can be modified as any ggplot by adding layers and/or elements with "+".
#'
#' @param data eem, eemlist, parafac or data.frame. The details are given under 'Details'
#' @param fill_max FALSE
#' @param ... parameters passed on to \code{ggplot}
#'
#' @details The data can be of different sources:
#'     eem: a single EEM pectrum is plotted
#'     eemlist: all spectra of the samples are plotted in one facet plot
#'     data.frame: a data.frame containing EEM data. Can be created by e.g. \code{as.data.frame.eem}
#'
#' @return a ggplot object
#'
#' @export
#' @import ggplot2 dplyr tidyr eemR
#' @importFrom grDevices rainbow
#'
#' @examples
#' ## plotting one distinct sample
#' data(eem_list)
#' eem <- eem_extract(eem_list,c("sample6","sample7"),keep=TRUE)
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
  ggeem(table,fill_max=fill_max)
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
  table <- data  %>% eempf_comp_mat() #eem_list
  table <- lapply(table %>% names(),function(name){
    table[[name]] %>% mutate(sample = name)
  }) %>% bind_rows()
  #filename <- paste0('EEM_PARAFAC_components_',suffix,format(Sys.time(), "%Y%m%d_%H%M%S"))
  ggeem(table,fill_max=fill_max,...)
}

#' @rdname ggeem
#' @export
ggeem.data.frame <- function(data,fill_max=FALSE,...)
{
  table <- data
  breaksx <- table %>% select(ex) %>% unique() %>% filter(as.numeric(ex)%%50 == 0) %>% unlist()
  breaksy <- table %>% select(em) %>% unique() %>% filter(as.numeric(em)%%50 == 0) %>% unlist()
  if(!is.numeric(fill_max)){
    fill_max <- table$value %>% max(na.rm=TRUE)
  }
  #values_fill <- seq(0,fill_max,length.out = 55)
  plot <- table %>% ggplot(...)+
    geom_raster(aes(x=ex,y=em,fill=value))+ #,interpolate=TRUE
    scale_fill_gradientn(colours=rainbow(75)[51:1],limits = c(0,fill_max))+
    #scale_fill_gradient2(low="blue",mid="yellow",high="red")
    scale_x_discrete(breaks = breaksx) +
    scale_y_discrete(breaks = breaksy) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))+
    facet_wrap(~sample)
  if(table$value %>% min(na.rm=TRUE) < 0){
    vals <- c(table$value %>% min(na.rm=TRUE),seq(from=0,to=fill_max,length.out = 51))
    vals <- (vals - min(vals))/diff(range(vals))
    plot <- plot +
      scale_fill_gradientn(colours=c(rainbow(75)[58],rainbow(75)[51:1]),values=vals,limits = c(table$value %>% min(na.rm=TRUE),fill_max))
  }
  plot
}


#' Plot fluorescence data from several samples split into several plots.
#'
#' @param data fluorescence data of class eemlist
#' @param spp number of samples per plot
#'
#' @return list of ggplots
#' @export
#'
#' @examples
#' \donttest{
#' data(eem_list)
#' eem_overview_plot(eem_list,spp=3)
#' }
eem_overview_plot <- function(data,spp = 8){
  #data <- eem_list
  ppp <- data %>% length()/spp
  fill_max <- data %>% eem_scale_ext() %>% .[2]
  #print(fill_max)
  ov_plot <- lapply(1:ceiling(ppp),function(pos){
    #pos <- 2
    data[(spp*(pos-1)+1):(spp*pos)] %>%
      `attr<-`("class", "eemlist") %>%
      ggeem(fill_max=fill_max)
  })
  ov_plot
}

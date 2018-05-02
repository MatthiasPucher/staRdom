## function to do a full_join on a list of tables
#' Full join of a list of data frames.
#'
#' @param df_list list of data frames to by joined
#' @param by character vector containing information how to join data frames. Format to be according to by in \code{\link[dplyr]{full_join}}. Each data frame has to contain the column(s) used for joining.
#'
#' @return The joint data frame.
#'
#' @seealso \code{\link[dplyr]{full_join}}
#'
#' @import dplyr
#' @export
#'
#' @examples
#' a <- data.frame(what=letters[1:5],a=c(1:5))
#' b <- data.frame(what=letters[1:5],b=c(7:11))
#' c <- data.frame(what=letters[1:5],c=c(20:24))
#'
#' df_list <- list(a,b,c)
#'
#' list_join(df_list,by="what")
list_join <- function(df_list,by){
  df <- df_list[[1]]
  for(n in 2:length(df_list)){
    df <- full_join(df,df_list[[n]],by=by)
  }
  df
}

## make data.frame from eem-metrix and optionally gather values
#' Converting EEM data from class eem to data.frame.
#'
#' @param x blabla
#' @param row.names asfas
#' @param optional ignored
#' @param ... ignored
#' @param gather logical, says whether data.frame is returned with excitation wavelength as column names or as values of a column. If the data is gathered, the sample name is added as value in a calumn
#'
#' @return A data frame containing the EEM data.
#'
#' @import dplyr tidyr
#' @export
#'
#' @examples
#' data(eem_list)
#' as.data.frame(eem_list[[1]])
#' as.data.frame(eem_list[[1]],gather=FALSE)
as.data.frame.eem <- function(x, row.names = NULL, optional = FALSE, ..., gather = TRUE){
  eem_d <- x$x %>% data.frame()
  colnames(eem_d) <- x$ex
  row.names(eem_d) <- x$em
  if (gather == TRUE){
    eem_d <- eem_d  %>%
      mutate(em = rownames(.)) %>%
      gather(ex,value, -em) %>%
      mutate(sample = paste0(x$sample))
  }
  eem_d
}


## remove ".FD3" or something else from sample name
#' Replace matched patterns in sample names
#'
#' @description Sample names in eemlist can be altered.
#'
#' @param eem_list data of class eemlist
#' @param pattern character vector containing pattern to look for.
#' @param replacement character vector of replacements. Has to have the same length as pattern
#'
#' @details \code{\link[stringr]{str_replace_all}} from package stringr is used for the replacement. Please read the corresponding help for further options.
#'
#' @return An eemlist.
#'
#' @seealso \code{\link[stringr]{str_replace_all}}
#'
#' @import dplyr tidyr eemR
#' @importFrom stringr str_replace_all
#' @export
#'
#' @examples
#' eem_names(eem_list)
#'
#' eem_list <- eem_name_replace(eem_list,"sample","Sample")
#' eem_names(eem_list)
eem_name_replace <- function(eem_list, pattern, replacement){
  eem_list <- lapply(eem_list,function(eem){
    eem$sample <- eem$sample %>%
      stringr::str_replace_all(pattern,replacement) #%>%
    #stringr::str_replace_all(c("\\.\\." = ".", "\\.$" = ""))
    eem
  })
  class(eem_list) <- "eemlist"
  eem_list
}


## extract part of eem matrix
#' Cut EEM data matching a given wavelength range
#'
#' @param data EEM data as eemlist
#' @param ex optional desired range of excitation wavelength
#' @param em optional desired range of emission wavelength
#'
#' @return An eemlist of reduced spectra size.
#' @import eemR
#' @export
#'
#' @examples
#' data(eem_list)
#' eem_range(eem_list,ex = c(250,Inf),em = c(280,500))
eem_range <- function(data,ex = c(0,Inf),em = c(0,Inf)){
  data %>%
    eem_cut(ex = c(0,max(0,ex[1]-1)), em = c(0,max(0,em[1]-1)), exact=FALSE) %>%
    eem_cut(ex = c(ex[2]+1,Inf), em = c(em[2]+1,Inf), exact=FALSE)
}

## get biggest eem spectra covered by each of the samples
#' Determines the the biggest range of EEM spectrum where data is available from each sample.
#'
#' @param data eemlist
#'
#' @return list of numeric vector containing the biggest available range
#' @import dplyr
#' @export
#'
#' @examples
#' data(eem_list)
#' eem_getextreme(eem_list)
#'
#' eem_list <- eem_range(eem_list,ex = c(250,Inf),em = c(280,500))
#' eem_getextreme(eem_list)
eem_getextreme <- function(data){
  res <- lapply(c("em","ex"),function(e){
    range <- lapply(data,function(eem){
      eem[[e]]
    }) %>%
      lapply(range)
    min <- range %>% lapply(function(r) r[1]) %>% unlist() %>% max()
    max <- range %>% lapply(function(r) r[2]) %>% unlist() %>% min()
    c(min,max)
  })
  names(res) <- c("em","ex")
  res
}


## get range of excitation values of sample set
#' Determine the range of fluorescence values in a set of samples
#'
#' @param data eemlist containing the EEM data
#'
#' @return numeric vector
#' @import dplyr
#' @export
#'
#' @examples
#' data(eem_list)
#' eem_scale_ext(eem_list)
eem_scale_ext <- function(data){
  range <- lapply(data,function(eem){
    eem$x
  }) %>%
    lapply(range,na.rm=TRUE,finite=TRUE) %>%
    unlist() %>%
    range(na.rm=TRUE)
  range
}

## reduce all samples to biggest wavelength range covered by all samples
#' Reduce wavelength range of EEM spectra to widest available in the whole sample set.
#'
#' @param data data of EEM samples as eemlist
#' @param verbose states whether additional information is given in the command line
#'
#' @details This step is neccessary to perform a PARAFAC analysis which can only be calculated with spectra of similar range.
#'
#' @return eemlist with reduced spectral width
#'
#' @import dplyr
#' @export
#'
#' @examples
#' data(eem_list)
#' eem_red2smallest(eem_list)
eem_red2smallest <- function(data,verbose=FALSE){
  extr <- data %>% eem_getextreme()
  if(verbose) cat(paste0("Samples are cut to ex from ",extr[['ex']][1]," to ",extr[['ex']][2]," and em from ",extr[['em']][1]," to ",extr[['em']][2]," \n"))
  data %>% eem_range(ex=extr[["ex"]], em=extr[["em"]])
}

#' Load all eemlist obects saved in different Rdata files in folder
#'
#' @description Reads Rdata files with one eemlist each. eemlists are combined into one and returned.
#'
#' @param dir folder where RData files are saved
#'
#' @return eemlist
#' @export
#'
#' @import dplyr
#' @importFrom eemR eem_bind
#'
#' @examples
#' \dontrun{
#' # due to package size issues no example data is provided for this function
#' eem_import_dir("C:/some_folder/with_EEMS/only_Rdata_files")
#' }
eem_import_dir <- function(dir){
  eem_files <- dir(dir, pattern=".RData") %>%
    paste0(dir,"/",.)

  for(file in eem_files){
    file <- load(file)
    if(get(file) %>% class() == "eemlist"){
    if(exists("eem_list")) eem_list <- eem_bind(eem_list,get(file)) else eem_list <- get(file)
    } else {
      warning(paste0(file," is no object of class eemlist!"))
    }
    NULL
  }

  eem_list
}

#' Opens an R markdown template for an easy and userfriendly analysis of EEM data.
#'
#' @description In your default editor (e.g. RStudio), a Rmd file is opened. It consists of bloacks gathering the parameters and information needed and continues with a series of data corrections, peak picking and plots. Finally you get a report of your analysis, a table with the peaks and optional pngs of your fluorescence data. To continue working and keeping your settings, the file can be sa ved anywhere and reused anytime.
#'
#' @details Function does not work well in Windows. You might try file.edit(system.file("EEM_simple_analysis.Rmd", package = "staRdom"))
#'
#' @return A pdf report, a peak picking table and optional plots.
#'
#' @importFrom utils file.edit
#' @export
#'
#' @examples
#' \dontrun{
#' #
#' eem_easy()
#'
#' # this function fails very often, so you might use that:
#' file.edit(system.file("EEM_simple_analysis.Rmd", package = "staRdom"))
#' }
eem_easy <- function(){
  file.edit(system.file("EEM_simple_analysis.Rmd", package = "staRdom"))
}

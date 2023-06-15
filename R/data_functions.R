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
    df <- try(full_join(df,df_list[[n]],by=by), silent = TRUE)
    if(inherits(df,"try-error")) stop("Error in ",names(df_list)[n],": ",df)
  }
  df
}

## make data.frame from eem-matrix and optionally gather values
#' Converting EEM data from class eem to data.frame.
#'
#' @param x abc
#' @param row.names abc
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
#'
#' as.data.frame(eem_list[[1]])
#' as.data.frame(eem_list[[1]],gather=FALSE)
as.data.frame.eem <- function(x, row.names = NULL, optional = FALSE, gather = TRUE, ...){
  #x <- eem_list[[1]]
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

## Reduce all samples to wavelengths present in all samples.
#' Remove wavelengths, that are missing in at least one sample form the whole set.
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
#' require(dplyr)
#'
#' data(eem_list)
#'
#' eem_list_red <- eem_red2smallest(eem_list)
#'
#' # create an eemlist where data is missing
#' eem_list2 <- eem_exclude(eem_list,
#'     list("ex" = c(280,290,350),
#'          "em" = c(402,510),
#'          "sample" = c()))
#'
#' # modify names of samples with missing data
#' eem_names(eem_list2) <- paste0("x",eem_names(eem_list2))
#'
#' # combined the lists with and without missing data
#' eem_list3 <- eem_bind(eem_list,eem_list2)
#' #ggeem(eem_list3)
#'
#' # reduce the data in the whole sampleset to the smallest wavelengths that are present in all samples
#' eem_list4 <- eem_red2smallest(eem_list3)
#' #ggeem(eem_list4)
#'
eem_red2smallest <- function(data,verbose=FALSE){
  ems <- lapply(data,`[[`,"em")
  ems
  emspres <- ems %>%
    unlist() %>%
    unique()
  emsTF <- lapply(ems,function(em){
    emspres %in% em
  })
  emspresTF <- emsTF %>%
    unlist() %>%
    matrix(nrow = length(emsTF), byrow = TRUE) %>%
    apply(2,all)

  exs <- lapply(data,`[[`,"ex")
  exspres <- exs %>%
    unlist() %>%
    unique()
  exsTF <- lapply(exs,function(ex){
    exspres %in% ex
  })
  exspresTF <- exsTF %>%
    unlist() %>%
    matrix(nrow = length(exsTF), byrow = TRUE) %>%
    apply(2,all)

  eem_exclude(data,list(ex = exspres[!exspresTF], em = emspres[!emspresTF]))
}

#' Load all eemlist obects saved in different Rdata or RDa files in a folder.
#'
#' @description Reads Rdata and RDa files with one eemlist each. The eemlists are combined into one and returned.
#'
#' @param dir folder where RData files are saved
#' @param verbose logical, set TRUE to show more information during import
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
#' # eem_import_dir("C:/some_folder/with_EEMS/only_Rdata_files")
#' }
eem_import_dir <- function(dir, verbose = FALSE){
  eem_files <- dir(dir, pattern=".RData$|.RDa$", ignore.case = TRUE) %>%
    paste0(dir,"/",.)

  for(file in eem_files){
    if(verbose) cat(paste0(file, " is read..."),fill=TRUE)
    file <- load(file)
    if(get(file) %>% inherits("eemlist")){
      if(exists("eems")) eems <- eem_bind(eems,get(file)) else eems <- get(file)
      if(verbose) cat("completed.",fill=TRUE)
    } else {
      warning(paste0(file," is no object of class eemlist!"))
    }
    NULL
  }

  eems
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

#' Load original data from the drEEM tutorial and return it as eemlist
#'
#' @return eemlist
#' @export
#'
#' @import dplyr tidyr
#' @importFrom R.matlab readMat
#'
#' @examples
#' \donttest{
#' # Reading MATLAB files from recent versions like the demo dataset from drEEM
#' # can cause problems if the R installation lacks UTF8 support in iconv.
#' # Therefore, we use try() in the example. If you encounter related problems,
#' # please refer to the help for R.matlab::readMat() for details.
#'
#' eem_list <- try(eem_load_dreem(), silent = FALSE)
#' eem_list
#' }
eem_load_dreem <- function(){
  dreem_raw <- tempfile()
  download.file("https://gitlab.com/dreem/drEEM/-/raw/master/tutorials_demos/datasets/PortSurveyData_corrected.mat?inline=false",dreem_raw)
  dreem_data <- dreem_raw %>%
    readMat()
  unlink(dreem_raw)

  eem_list <- lapply(dreem_data$filelist.eem, function(file){
    #file <- dreem_data$filelist.eem[1]
    n <- which(dreem_data$filelist.eem == file)
    file <- file %>% gsub("^\\s+|\\s+$", "", .) %>% # trim white spaces in filenames
      sub(pattern = "(.*)\\..*$", replacement = "\\1", .) # remove file extension from sample name
    eem <- list(sample = file,x = dreem_data$XcRU[n,,] %>% as.matrix(),ex = dreem_data$Ex %>% as.vector(), em = dreem_data$Em.in %>% as.vector(),location = "drEEM/dataset/")
    class(eem) <- "eem"
    attr(eem, "is_blank_corrected") <- TRUE
    attr(eem, "is_scatter_corrected") <- FALSE
    attr(eem, "is_ife_corrected") <- TRUE
    attr(eem, "is_raman_normalized") <- TRUE
    attr(eem, "manufacturer") <- "unknown"
    eem
  }) %>%
    `class<-`("eemlist")
}

#' Import EEMs from generic csv tables (deprecated)
#'
#' @description This function is deprecate, please use \code{\link[eemR]{eem_read}}(..., import_function = \code{\link{eem_csv}}) or eem_read(..., import_function = \code{\link{eem_csv2}}) instead.
#' EEM data is loaded from generic files. First column and first row contains wavelength values. The other values are to be plain numbers. \code{\link[data.table]{fread}} is used to read the table. It offers a lot of helpful functions (e.g. skipping any number n of header lines by adding `skip = n`)
#'
#' @param path path to file(s), either a filename or a folder
#' @param col either "ex" or "em", what wavelengths are in the columns
#' @param recursive logical, whether directories are loaded recursively
#' @param is_blank_corrected logical, whether blank correction was done
#' @param is_scatter_corrected logical, wether scatters were corrected
#' @param is_ife_corrected logical, wether inner-filter effect correction was done
#' @param is_raman_normalized logical, wether raman normalisation applied
#' @param manufacturer string specifying manufacturer of instrument
#' @param ... parameters from other functions, currently not used
#'
#' @export
#'
#' @import dplyr tidyr
#'
#' @examples
#' eems <- system.file("extdata/EEMs",package="staRdom")
#' eem_list <- eem_read_csv(eems)
#'
#' eem_list
eem_read_csv <- function(path, col = "ex", recursive = TRUE, is_blank_corrected = FALSE, is_scatter_corrected = FALSE, is_ife_corrected = FALSE, is_raman_normalized = FALSE, manufacturer = "unknown", ...){
  #path <- "./inst/extdata/EEMs"
  warning("This function is deprecated, please use eem_read(..., import_function = eem_csv) or eem_read(..., import_function = eem_csv2) instead!")
  if(col == "em"){
    csv_func <- eem_csv2
  } else {
    csv_func <- eem_csv
  }
  eems <- eem_read(path, recursive = recursive, import_function = csv_func)

  lapply(eems, function(eem){
    class(eem) <- "eem"
    attr(eem, "is_blank_corrected") <- is_blank_corrected
    attr(eem, "is_scatter_corrected") <- is_scatter_corrected
    attr(eem, "is_ife_corrected") <- is_ife_corrected
    attr(eem, "is_raman_normalized") <- is_raman_normalized
    attr(eem, "manufacturer") <- manufacturer
    eem
  }) %>%
    `class<-`("eemlist")
}

#' Applying functions on EEMs
#'
#' @param data eemlist to be modified
#' @param func a function to be applied on the data.
#' @param return either "eemlist" or "value"
#' @param ... additional arguments passed on to func
#'
#' @details The EEMs are passed on as first argument to func. Additionally, the vector of excitation wavelengths is passed on as \code{ex} and the emission wavelengths as \code{em}. Therefore, the supplied function has to allow these arguments. The easiest way would be \code{...} (see example).
#'
#' @return eemlist or list
#' @export
#'
#' @examples
#' ## define a function, that would divide a matrix by its maximum
#' # more general, if you want to return a valid eemlist (see below),
#' # a matrix of the same size has to be returned
#' # ... is used as a placeholder for any argument, important: em and
#' # ex wavelengths are passed on, so the function needs to take them as arguments,
#' # even if they are not used
#' norm_max <- function(x, ...){
#'   x/max(x)
#' }
#'
#' # load example data
#' data("eem_list")
#'
#' # normalise eems by the function defined above
#' norm_eems <- eem_apply(eem_list,norm_max,"eemlist")
#'
#' # plot the results to see the difference
#' ggeem(norm_eems)
#'
#' # define another function. what values were used to
#' # multiply the eems with?
#' norm_fac <- function(x, ...){
#'   1/max(x)
#' }
#'
#' # return a list of factors used for normalisation
#' norm_factors <- eem_apply(eem_list,norm_fac,"value")
#'
#' unlist(norm_factors)
#'
#' # return list of em vectors.
#' # important: x needs to be in the first position, but
#' # is not used later!
#' extr_em <- function(x,em,...){
#'   em
#' }
#'
#' em_vectors <- eem_apply(eem_list,extr_em,"value")
#'
#' em_vectors
#'
eem_apply <- function(data, func, return = c("eemlist","value"),...){
  stopifnot(class(data) == "eemlist")
  #eem <- eem_list[[1]]
  ress <- lapply(data, function(eem){
    #extr_em(eem$x, ex = eem$ex, em = eem$em)
    eem$x <- func(eem$x, ex = eem$ex, em = eem$em,...)#,...
    if(return[1] == "eemlist"){
      res <- eem
    } else {
      res <- eem$x
    }
    res
  })
  if(return[1] == "eemlist"){
    class(ress) <- "eemlist"
  }
  ress
}

#' Export all samples of an eem_list
#'
#' @param file path to directory (csv format) or file (Matlab format)
#' @param format either "csv" or "mat" to specify export format
#' @param ... one or more eem_list objects
#'
#' @return 0 on successful export
#'
#' @import dplyr
#' @importFrom readr write_csv
#'
#' @export
#'
#' @examples
#' # create temporary directory to write out
#' file <- paste0(tempdir(),"/eem_export/")
#' dir.create(file)
#'
#' # run eem_export to write one csv file for each sample
#' eem_export(file, format = "csv", eem_list)
#'
#' # show content of output directory
#' dir(file)
#'
eem_export <- function(file, format = c("csv","mat"), ...){
  if(format[1] == "mat"){
    eem_export_matlab(file = file, ...)
  } else if(format[1] == "csv"){
    if(!file.exists(file) | !file.info(file)$isdir){
      stop(file, " is not a valid folder!")
    } else{
      eem <- eem_bind(...)
      if(!identical(eem_names(eem), eem_names(eem) %>% unique())){
        eem_names(eem) <- eem_names(eem) %>% make.unique()
        warning("Sample names were not unique! Unique names were produced automatically.")
      }
      lapply(eem, function(e){
        data.frame(e$x) %>%
          setNames(e$ex) %>%
          mutate(em = e$em) %>%
          select(em, as.character(e$ex)) %>%
          write_csv(file = paste0(file,"/",e$sample,".csv"))
      })
    }
  } else {
    stop("The supplied format is not (yet) available!")
  }
  invisible(0)
}

#' Export samples in an EEM list to a single csv files
#'
#' @param eem_list EEM data as eemlist
#' @param output path to folder where csv files are exported to
#' @param ... additional arguments
#'
#' @return returns the exported EEMs as a list of data.frames
#' @export
#'
#' @examples
#' data(eem_list)
#'
#' output <- tempdir()
#' output
#' a <- eem_write_csv(eem_list, output)
eem_write_csv <- function(eem_list, output, ...){
  stopifnot(dir.exists(output))
  stopifnot(inherits(eem_list,"eemlist"))
  lapply(eem_list, function(eem){
    cat("Exporting ",eem$sample,"...",fill=TRUE)
    tab <- data.frame(wavelength = eem$em) %>%
      bind_cols(eem$x %>%
      `colnames<-`(eem$ex))
    filename = paste0(eem$sample,".csv")
    fil <- paste0(output,"/",filename)
    if(file.exists(fil)){
      warning(fil," existed and was overwritten!")
    }
    write_csv(x = tab, file = fil)
  }) %>%
    invisible()
}


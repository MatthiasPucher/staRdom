## interpolation of missing values due to removing scatter, pchip algorithm done columnwise (for each excitation wavelength)
#' Missing values are interpolated within EEM data
#'
#' @description Missing EEM data can be interpolated. Usually it is the result of removing scatter or other parts where noise is presumed. Different interpolation algorithms can be used (see details).
#'
#' @param data object of class eemlist with spectra containing missing values
#' @param cores specify number of cores for parallel computation
#' @param type numeric 0 to 4 or TRUE which resembles type 1
#' @param nonneg logical, whether negative values should be replaced by 0
#' @param verbose logical, whether more information on calculation should be provided
#' @param extend logical, whether data is extrapolated using type 1
#' @param ... arguments passed on to other functions (pchip, na.approx, mba.points)
#'
#' @details The types of interpolation are (0) setting all NAs to 0, (1) spline interpolation with \code{\link[MBA]{mba.points}}, (2) excitation and emission wavelength-wise interpolation with \code{\link[pracma]{pchip}} and subsequent mean, (3) excitation wavelength-wise interpolation with \code{\link[pracma]{pchip}} and (4) linear interpolation in 2 dimensions with \code{\link[zoo]{na.approx}} and again subsequent mean calculation. Calculating the mean is a way of ensuring NAs are also interpolated where missing boundary values would make that impossible. Using type = 1, extrapolation can be suppressed by adding the argument extend = FALSE.
#'
#' @return object of class eemlist with interpoleted spectra.
#'
#' @seealso \code{\link[pracma]{pchip}}, \code{\link[MBA]{mba.points}}, \code{\link[zoo]{na.approx}}
#'
#' @references Elcoroaristizabal, S., Bro, R., García, J., Alonso, L. 2015. PARAFAC models of fluorescence data with scattering: A comparative study. Chemometrics and Intelligent Laboratory Systems, 142, 124-130
#' \url{https://doi.org/10.1016/j.chemolab.2015.01.017}
#'
#' @import dplyr tidyr doParallel
#' @importFrom stats na.omit
#' @importFrom pracma pchip
#' @importFrom zoo na.approx
#' @importFrom MBA mba.points
#' @export
#'
#' @examples
#' \donttest{
#' data(eem_list)
#' eem_list <- eem_list[1:6]
#' class(eem_list) <- "eemlist"
#'
#' remove_scatter <- c(TRUE, TRUE, TRUE, TRUE)
#'
#' remove_scatter_width = c(15,10,16,12)
#'
#' eem_list <- eem_rem_scat(eem_list,remove_scatter,remove_scatter_width)
#'
#' eem_list <- eem_interp(eem_list)
#'
#' ggeem(eem_list)
#'
#' eem_list2 <- eem_setNA(eem_list,ex=200:280,interpolate=FALSE)
#'
#' ggeem(eem_list2)
#'
#' eem_list3 <- eem_interp(eem_list2,type=1,extend = TRUE)
#'
#' ggeem(eem_list3)
#'
#' eem_list3 <- eem_interp(eem_list2,type=1,extend = FALSE)
#'
#' ggeem(eem_list3)
#'
#'
#' }
eem_interp <- function(data,cores = parallel::detectCores(logical = FALSE), type = TRUE, verbose = FALSE, nonneg=TRUE, extend = FALSE,...){
  cl <- makeCluster(spec = cores, type = "PSOCK")
  doParallel::registerDoParallel(cl)
  if(verbose){
    cat("interpolating missing data in",length(data),"EEMs", fill=TRUE)
  }
  eem_list <- foreach(i = 1:length(data)) %dopar% {
    #data <- eem_list
    #i <- 4
    eem <- data[[i]]
    if(type == 4){
      eem$x <- cbind(zoo::na.approx(eem$x,...),t(zoo::na.approx(t(eem$x),...))) %>% array(c(nrow(eem$x),ncol(eem$x),2)) %>%
        apply(1:2, mean, na.rm = TRUE)
            }
    if(type == 1 | type == TRUE){
      x <- eem$x %>%
        data.frame() %>%
        `colnames<-`(eem$ex) %>%
        `rownames<-`(eem$em) %>%
        mutate(em = eem$em) %>%
        gather("ex","z",-em) %>%
        mutate_all(as.numeric)
      x2 <- x %>%
        filter(!is.na(z))
      x3 <- MBA::mba.points(xyz = x2 %>% select(em,ex,z), xy.est = expand.grid(em = eem$em, ex = eem$ex), verbose = verbose, extend = extend, ...)
      eem$x[is.na(eem$x)] <- x3$xyz.est[,3] %>% matrix(nrow = nrow(eem$x), ncol = ncol(eem$x)) %>% .[is.na(eem$x)] #%>% pmin(max(eem$x,na.rm=TRUE)) %>% pmax(min(eem$x,na.rm=TRUE))#matrix(x3,nrow=length(eem$em), ncol = length(eem$ex))
      }
    if(type == 2){
    x1 <- try(eem$x %>% apply(1,function(row) pracma::pchip(xi=eem$ex[!is.na(row)],yi=row %>% na.omit(),x=eem$ex),...) %>% t(),silent=TRUE)

    x2 <- try(eem$x %>% apply(2,function(col) pracma::pchip(xi=eem$em[!is.na(col)],yi=col %>% na.omit(),x=eem$em),...),silent=TRUE)
    if(class(x1)=="try-error" & class(x2)=="try-error") warning(eem$sample," could not be interpolated!")
    if(class(x1)=="try-error") x1 <- matrix(NA,nrow(eem$x),ncol(eem$x))
    if(class(x2)=="try-error") x2 <- matrix(NA,nrow(eem$x),ncol(eem$x))
    eem$x <- cbind(x1,x2) %>% array(c(nrow(x1),ncol(x2),2)) %>%
      apply(1:2, mean, na.rm = TRUE)
    }
    if(type == 3){
      eem$x <- eem$x %>% apply(2,function(col) pracma::pchip(xi=eem$em[!is.na(col)],yi=col %>% na.omit(),x=eem$em,...))
    }
    if(type == 0){
      eem$x[is.na(eem$x)] <- 0
    }
    if(nonneg) eem$x[eem$x<0] <- 0
    eem
  }
  stopCluster(cl)

  class(eem_list) <- "eemlist"
  eem_list
}

#' set parts of specific samples to NA and optionally interpolate these parts
#'
#' @param eem_list EEMs as eemlist
#' @param sample optional, names or indices of samples to process
#' @param em optional, emission wavelengths to set NA
#' @param ex optional, excitation wavelengths to set NA
#' @param interpolate FALSE, 1 or 2, interpolate NAs or not, 2 different methods, see \code{\link[staRdom]{eem_interp}}
#' @param ... arguments passed on to \code{\link[staRdom]{eem_interp}}
#'
#' @details Samples and wavelengths are optional and if not set all of them are considered in setting data to NA. Wavelengths can be set as vectors containing more than the wavelengths present in the data. E.g. 230:250 removes all wavelengths between 230 and 250 if present. Data is best interpolated if it does not reach data boundaries. Please check the results otherwise as in some cases the interpolation might not produce meaningful data.
#'
#' @return eemlist
#' @export
#'
#' @import dplyr
#'
#' @examples
#' data(eem_list)
#' eem <- eem_list[1:9]
#' class(eem) <- "eemlist"
#'
#' ggeem(eem)
#'
#' eem_list2 <- eem_setNA(eem,ex=200:280,em=500:600, interpolate=FALSE)
#' ggeem(eem_list2)
eem_setNA <- function(eem_list,sample=NULL,em=NULL,ex=NULL,interpolate = TRUE,...){
  if(is.null(sample)){
    sample <- eem_names(eem_list)
  }
  if(is.numeric(sample)){
    sample <- eem_names(eem_list)[sample]
  }
  eem_list <- lapply(eem_list,function(eem){
    if(eem$sample %in% sample){
      if(is.null(ex)) ex2 <- 1:ncol(eem$x) else ex2 <- which(eem$ex %in% ex)
      if(is.null(em)) em2 <- 1:nrow(eem$x) else em2 <- which(eem$em %in% em)
      eem$x[em2,ex2] <- NA
    }
    eem
  }) %>%
    `class<-`("eemlist")
  if(interpolate != FALSE){
    eem_list[which(eem_names(eem_list) %in% sample)] <- eem_interp(eem_list[which(eem_names(eem_list) %in% sample)],type = interpolate,...) %>%
      `class<-`("eemlist")
  }
  eem_list
}

#' Exclude complete wavelengths or samples form data set
#'
#' @description Outliers in all modes should be avoided. With this functions excitation or emission wavelengths as well as samples can be removed completely from your sample set.
#'
#' @param eem_list object of class eemlist
#' @param exclude list of three vectors, see details
#' @param verbose states whether additional information is given in the command line
#'
#' @details The argument exclude is a named list of three vectors. The names must be "ex", "em" and "sample". Each element contains a vector of wavelengths or sample names that are to be excluded from the data set.
#'
#' @return object of class eemlist
#' @export
#'
#' @examples
#' data(eem_list)
#'
#' exclude <- list("ex" = c(280,285,290,295),
#' "em" = c(),
#' "sample" = c("667sf", "494sf")
#' )
#'
#' eem_list_ex <- eem_exclude(eem_list, exclude)
eem_exclude <- function(eem_list, exclude = list,verbose=FALSE){
  ex_exclude <- exclude[["ex"]]
  em_exclude <- exclude[["em"]]
  sample_exclude <- exclude[["sample"]]
  if(!is.null(sample_exclude)){
    sample_exclude <- paste0("^",sample_exclude,"$")
    eem_list <- eem_extract(eem_list,sample_exclude,verbose=verbose)
  }
  eem_list <- lapply(eem_list,function(eem){
    eem$x <- eem$x[!eem$em %in% em_exclude,!eem$ex %in% ex_exclude]
    eem$ex <- eem$ex[!eem$ex %in% ex_exclude] #%>% length()
    eem$em <- eem$em[!eem$em %in% em_exclude] #%>% length()
    eem
  })
  if(!is.null(ex_exclude) & verbose) cat(paste0("Removed excitation wavelength(s): ",paste0(ex_exclude %>% sort(),collapse=", ")),fill=TRUE)
  if(!is.null(em_exclude) & verbose) cat(paste0("Removed emission wavelength(s): ",paste0(em_exclude %>% sort(),collapse=", ")),fill=TRUE)
  class(eem_list) = "eemlist"
  eem_list
}

## raman normalisation either by blank, by one value or by a set of values, one for each sample
#' Wrapper function to eem_raman_normalisation (eemR).
#'
#' @description Usually Raman normalisation is done with fluorescence data from a blank sample. Sometimes you already know a value for the Raman area. This function can do both.
#'
#' @param data fluorescence data of class eemlist
#' @param blank defines how Raman normalisation is done (see 'Details')
#'
#' @details Possible values for blank:
#'
#'     "blank": normalisation is done with a blank sample. Please refer to \code{\link[eemR]{eem_raman_normalisation}}.
#'
#'     numeric: normalisation is done with one value for all samples.
#'
#'     data frame: normalisation is done with different values for different samples. Values are taken from a data.frame with sample names as rownames and one column containing the raman area values.
#'
#' @return fluorescence data of class eemlist
#'
#' @import eemR
#' @export
#'
#' @examples
#' data(eem_list)
#' # correction by blank
#' eems_bl <- eem_raman_normalisation2(eem_list,blank="blank")
#'
#' # correction by value
#' eems_num <- eem_raman_normalisation2(eem_list,blank=168)
eem_raman_normalisation2 <- function(data, blank="blank"){
  if(blank == "blank"){
    res_list <- try(data %>% eem_raman_normalisation(),silent=TRUE)
    if (class(res_list) == "try-error") {
      warning(res_list)
      stop("There was a problem with raman normalisation, please check the presence of a blank sample and the parameter blank!")
    }
  } else if(is.data.frame(blank) & class(data) == "eemlist"){

    res_list <- lapply(data,function(eem){
      rar <- blank[eem$sample %>% make.names(),]
      if(!is.na(rar)){
        eem$x <- eem$x/rar
        attr(eem, "is_raman_normalized") <- TRUE
      } else {
        warning(paste0("Sample ",eem$sample," was not normalised because an according raman area was not found in the table!"))
      }
      class(eem) <- "eem"
      return(eem)
    })
    class(res_list) <- class(data)
  } else if(is.numeric(blank) & length(blank)==1 & class(data) == "eemlist"){
    res_list <- lapply(1:length(data),function(i){
      data[[i]]$x <- data[[i]]$x/blank
      attr(data[[i]], "is_raman_normalized") <- TRUE
      return(data[[i]])
    })
    class(res_list) <- class(data)
  } else {
    stop("First argument must be of class eemlist, second argument must be 'blank' for a blank correction, a number for correction with this number or a vector with a number for each sample")
  }
  res_list
}

#' Wrapper function to allow eem_inner_filter_effect (eemR) handling different cuvette lengths.
#'
#' @description Calls \code{\link[eemR]{eem_inner_filter_effect}} for each sample to use different cuvette lengths.
#'
#' @param data fluorescence data of class eemlist
#' @param abs_data absorbance data
#' @param cuvl length of cuvette of absorption measurment in cm. Either a number or a data frame. Row names of data frame have to be similar to sample names in data. This is ignored, if unit is "absorption".
#' @param unit unit of absorbance data. Either "absorbance" or "absorption".
#'
#' @return fluorescence data of class eemlist
#' @import eemR
#' @export
#'
#' @examples
#' folder <- system.file("extdata/cary/scans_day_1", package = "eemR") # load example data
#' eem_list <- eem_read(folder, import_function = "cary")
#' data(absorbance)
#'
#' eem_list <- eem_ife_correction(data = eem_list, abs_data = absorbance,
#'     cuvl = 5, unit = "absorbance")
eem_ife_correction <- function(data, abs_data, cuvl = NULL, unit = c("absorbance","absorption")){
  if(!unit[1] %in% c("absorbance","absorption")) stop("Unit must be either 'absorbance' or 'absorption'!")
  if(unit[1] == "absorbance" & !is.numeric(cuvl)) stop("Please specify a valid cuvette length!")
  if(unit[1] == "absorption"){
    abs_data <- abs_data %>%
      mutate_at(vars(-wavelength),`*`,(1/log(10)))
    cuvl <- 100
  }
  eem_list <- data %>% lapply(function(eem1){
    if(is.data.frame(cuvl)) cl <- cuvl[eemnam,] else cl <- cuvl
    eem1 <- list(eem1)
    nam <- eem1[[1]]$sample
    class(eem1) <- "eemlist"
    if(nam %in% colnames(abs_data)){
      eem1 <- eem1 %>% eem_inner_filter_effect(absorbance = na.omit(abs_data[c("wavelength",nam)]),pathlength=cl)
      eem1[[1]]$sample <- nam
    } else {
      warning(paste0("No ",unit," data was found for sample ",nam,"!"))
    }
    eem1
  }) %>%
    unlist(recursive = FALSE)
  class(eem_list) <- "eemlist"
  eem_list
}

## multiplying the matrix with a dilution factor
#' Modifying fluorescence data according to dilution.
#'
#' @description If samples were diluted before measuring, a dilution factor has to be added to the measured data. This function can do that by either multilpying each sample with the same value or using a data frame with different values for each sample.
#'
#' @param data fluorescence data with class eemlist
#' @param dilution dilution factor(s), either numeric value or data frame. Row names of data frame have to be similar to sample names in eemlist.
#'
#' @return fluorescence data with class eemlist
#' @export
#'
#' @examples
#' data(eem_list)
#'
#' eem_list <- eem_dilution(eem_list,dilution=5)
eem_dilution <- function(data,dilution=1){
  if(((is.numeric(dilution) & length(dilution)==1) | is.data.frame(dilution)) & class(data) == "eemlist"){
    res_list <- lapply(1:length(data),function(i){
      if(is.data.frame(dilution)){
        data[[i]]$x <- data[[i]]$x*dilution[data[[i]]$sample,]
      } else {
        data[[i]]$x <- data[[i]]$x*dilution
      }
      data[[i]]
    })
    class(res_list) <- class(data)
  } else {
    stop("First argument must be of class eemlist, second argument must be a number or a data.frame to multiply with for dilution.")
  }
  res_list
}

## smooth eem matrix by rolling mean along excitation wavelength
#' Smooth fluorescence data by calculating rolling mean along excitation wavelengths.
#'
#' @param data fluorescence data of class eemlist
#' @param n width of rolling mean window in nm
#'
#' @return eemlist with smoothed data
#'
#' @importFrom zoo rollmean
#' @export
#'
#' @examples
#' \donttest{
#' data(eem_list)
#'
#' eem_list <- eem_smooth(eem_list,n=4)
#' }
eem_smooth <- function(data,n = 4){
  n <- n/2
  data <- lapply(data,function(eem){
    k <- which(eem$em[1] + n >= eem$em) %>% max()
    eem$x <- eem$x %>% apply(2,function(col) col %>% rollmean(k=k,fill=c(0,0,0)))
    eem
  })
  class(data) <- "eemlist"
  data
}


## function to remove scattering
#' Remove Raman and Rayleigh scattering in fluorescence data
#'
#' @description Wrapper function to remove several scatterings in one step using \code{\link[eemR]{eem_remove_scattering}}.
#'
#' @param data object of class eemlist
#' @param remove_scatter logical vector. The meanings of the vector are "raman1", "raman2", "rayleigh1" and "rayleigh2" scattering. Set \code{TRUE} if certain scattering should be removed.
#' @param remove_scatter_width numeric vector containing width of scattering to remove. If there is only one element in this vector, each this is the width of each removed scattering. If there are 4 values, differnt widths are used ordered by "raman1", "raman2", "rayleigh1" and "rayleigh2".
#' @param interpolation logical, optionally states whether interpolation is done right away
#' @param cores optional, CPU cores to use for interpolation
#' @param verbose logical, provide additional information
#'
#' @return eemlist
#'
#' @import eemR dplyr tidyr
#' @export
#'
#' @examples
#' data(eem_list)
#'
#' remove_scatter <- c(TRUE, TRUE, TRUE, TRUE)
#'
#' remove_scatter_width = c(15,10,16,12)
#'
#' eem_rem_scat(eem_list,remove_scatter,remove_scatter_width)
eem_rem_scat <- function(data,remove_scatter,remove_scatter_width = 10, interpolation = FALSE, cores = parallel::detectCores(logical=FALSE), verbose = FALSE)
{
  if(data %>% class != "eemlist") stop("first argument has to be a list of eem samples!")
  if(!is.numeric(remove_scatter_width)) stop("removed scatter slot width has to be numeric!")
  if(length(remove_scatter_width) == 1) remove_scatter_width <- rep(remove_scatter_width,4)
  if(remove_scatter[1]) data <- data %>% eem_remove_scattering(type="raman",order=1,width=remove_scatter_width[1])
  if(remove_scatter[2]) data <- data %>% eem_remove_scattering(type="raman",order=2,width=remove_scatter_width[2])
  if(remove_scatter[3]) data <- data %>% eem_remove_scattering(type="rayleigh",order=1,width=remove_scatter_width[3])
  if(remove_scatter[4]) data <- data %>% eem_remove_scattering(type="rayleigh",order=2,width=remove_scatter_width[4])

  if(interpolation != FALSE) data <- eem_interp(data, cores = cores, type = interpolation, verbose = verbose)

  data
}


#' Calculate raman area of EEM samples
#'
#' @param eem_list An object of class eemlist.
#' @param blanks_only logical. States whether all samples or just blanks will be used.
#' @param average logical. States whether samples will be averaged before calculating the raman area.
#'
#' @return data frame containing sample names, locations and raman areas
#' @export
#'
#' @details Code based on \code{\link[eemR]{eem_raman_normalisation}}.
#'
#' @import dplyr
#' @importFrom eemR eem_extract
#' @importFrom pracma interp2
#' @importFrom stats setNames
#'
#' @examples
#' folder <- system.file("extdata/EEMs",package="staRdom")
#' eem_list <- eem_read(folder, recursive = TRUE, import_function = eem_csv)
#' blank <- eem_extract(eem_list,sample ="blank", keep = TRUE)
#'
#' eem_raman_area(blank)
eem_raman_area <- function(eem_list, blanks_only = TRUE, average = FALSE){
  if(blanks_only){
    blank_names <- c("nano", "miliq", "milliq", "mq", "blank")
    eem_list <- eem_extract(eem_list, blank_names, keep = TRUE, ignore_case = TRUE,
                            verbose = FALSE)
    if (average) {
      n <- length(eem_list)
      message("A total of ", n, " sample EEMs will be averaged.")
      X <- Reduce("+", lapply(eem_list, function(x) x$x))
      X <- X/n
      eem_list <- eem_list[1]
      eem_list[[1]]$x <- X
      class(eem_list) <- "eemlist"
    }
  }

  location <- lapply(eem_list, function(eem){
    eem$location
  }) %>% unlist()

  ram_area <- lapply(eem_list,function(eem){
    #eem <- unlist(eem, recursive = FALSE)
    em <- seq(371, 428, by = 2)
    ex <- rep(350, length(em))
    fluo <- pracma::interp2(eem$ex, eem$em, eem$x, ex,
                            em)
    if (any(is.na(em)) | any(is.na(fluo))) {
      stop("NA values found in the samples. Maybe you removed scattering too soon?",
           call. = FALSE)
    }
    area <- sum(diff(em) * (fluo[-length(fluo)] + fluo[-1])/2)
  }) %>%
    unlist() %>%
    data.frame(sample = eem_names(eem_list), location = location, raman_area = .)

  ram_area
}


#' Multiply all EEMs with a matrix
#'
#' @param eem_list EEM data as eemlist
#' @param matrix either a vactor containing "l" and/or "u" or a matrix, see details.
#' @param value in case matrices "l" or "u" are used, this specifies the value to use in this areas. Usually this is 0 (default) or NA but any numeric value can be used.
#'
#' @details All EEMs must be of the same size. If matrix is of type matrix, it is used right away to multiply the EEMs. It has to be of the same size as the EEMs. If matrix is a vector containing "l", values below 1st order Rayleigh scattering are set to 0. If matrix contains "u", values above 2nd order Raman scattering are set to 0. If you want to remove wavelength ranges, take into consideration to use \code{\link[eemR]{eem_cut}} or \code{\link[staRdom]{eem_range}}.
#'
#'
#' @return eemlist
#' @export
#'
#' @import dplyr
#'
#' @examples
#' \donttest{
#' data(eem_list)
#' eem <- eem_list[1:9]
#' class(eem) <- "eemlist"
#'
#' ggeem(eem)
#'
#' eem_list_cut <- eem_matmult(eem,matrix=c("l"), value= NA)
#' ggeem(eem_list_cut)
#' }
eem_matmult <- function(eem_list,matrix = NULL,value = 0){
  if(is.matrix(matrix)){
    mat <- matrix
  } else {
    mat <- matrix(data = 1, nrow = nrow(eem_list[[1]]$x), ncol = ncol(eem_list[[1]]$x))
    if("l" %in% matrix | "u" %in% matrix){
      ex_mat <- matrix(eem_list[[1]]$ex, nrow = nrow(eem_list[[1]]$x), ncol = ncol(eem_list[[1]]$x), byrow=TRUE)
      em_mat <- matrix(eem_list[[1]]$em, nrow = nrow(eem_list[[1]]$x), ncol = ncol(eem_list[[1]]$x), byrow=FALSE)
    }
    if("l" %in% matrix){
      mat <- mat * (ex_mat < em_mat)
      mat[mat == 0] <- value
    }
    if("u" %in% matrix){
      raman_peaks <- 1/(1/ex_mat - 0.00036)
      mat <- mat * (2*raman_peaks > em_mat)
      mat[mat == 0] <- value
    }
  }

  lapply(eem_list, function(eem){
    mat[mat == 0] <- value
    eem$x <- eem$x * mat
    eem
  }) %>%
    `class<-`("eemlist")
}


#' EEM sample data is extended to include all wavelengths in all samples
#'
#' @description Compared to the whole sample set, wavelengths missing in some samples are added and set NA or interpolated. This can be especially helpful, if you want to combine data measured with different wavelength intervals in a given range.
#'
#' @param eem_list eemlist
#' @param interpolation logical, whether added NAs should be interpolated
#' @param ... arguments passed to eem_interp
#'
#' @import dplyr
#'
#' @return eemlist
#' @export
#'
#' @examples
#' library(dplyr)
#' data(eem_list)
#' eem_list <- eem_exclude(eem_list[1:5] %>%
#' `class<-`("eemlist"), exclude = list(em = c(318,322,326,550,438), ex = c(270,275))) %>%
#' eem_bind(eem_list[6:15] %>% `class<-`("eemlist"))
#' ggeem(eem_list)
#'
#' eem_extend2largest(eem_list) %>%
#'   ggeem()
eem_extend2largest <- function(eem_list, interpolation = FALSE,...){
  exs <- lapply(eem_list, function(eem) eem$ex) %>%
    unlist() %>%
    unique() %>%
    sort()
  ems <- lapply(eem_list, function(eem) eem$em) %>%
    unlist() %>%
    unique() %>%
    sort()
  res <- lapply(eem_list, function(eem){
    eemx <- matrix(nrow = length(ems), ncol = length(exs))
    eemx[which(ems %in% eem$em),which(exs %in% eem$ex)] <- eem$x
    eem$x <- eemx
    eem$em <- ems
    eem$ex <- exs
    eem
  }) %>%
    `class<-`("eemlist")
  if(interpolation != FALSE){
    res <- eem_interp(res,type = interpolation, ...)
  }
  res
}


#' Multiply EEMs with spectral correction vectors (Emission and Excitation)
#'
#' @param eem_list eemlist
#' @param Excor data frame, first column wavelengths, second column excitation correction
#' @param Emcor data frame, first column wavelengths, second column emission correction
#'
#' @return eemlist
#' @export
#'
#' @examples
#' eems <- system.file("extdata/EEMs",package="staRdom")
#' eem_list <- eem_read(eems, recursive = TRUE, import_function = eem_csv)
#'
#' excorfile <- system.file("extdata/CorrectionFiles/xc06se06n.csv",package="staRdom")
#' Excor <- data.table::fread(excorfile)
#' emcorfile <- system.file("extdata/CorrectionFiles/mcorrs_4nm.csv",package="staRdom")
#' Emcor <- data.table::fread(emcorfile)
#'
#' # adjust range of EEMs to cover correction vectors
#' eem_list <- eem_range(eem_list,ex = range(Excor[,1]), em = range(Emcor[,1]))
#'
#' eem_list_sc <- eem_spectral_cor(eem_list,Excor,Emcor)
eem_spectral_cor <- function(eem_list,Excor,Emcor){
  x <- eem_getextreme(eem_list)
  xEx <- range(Excor[,1])
  xEm <- range(Emcor[,1])
  if(min(xEx) > min(x$ex) & max(xEx) < max(x$ex)) stop("Excitation correction does not cover EEM excitation spectrum!")
  if(min(xEm) > min(x$em) & max(xEm) < max(x$em)) stop("Emission correction does not cover EEM emission spectrum!")
  eem_list <- lapply(eem_list,function(eem){
    #eem <- eem_list[[1]]
    Excor1 <- approx(x=Excor[[1]],y=Excor[[2]],xout=eem$ex)$y
    Emcor1 <- approx(x=Emcor[[1]],y=Emcor[[2]],xout=eem$em)$y
    mcor <- Emcor1 %*% t(Excor1)
    eem$x <- eem$x * mcor
    eem
  }) %>%
    `class<-`("eemlist")
}

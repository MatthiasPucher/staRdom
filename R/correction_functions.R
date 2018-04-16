## interpolation of missing values due to removing scatter, pchip algorithm done columnwise (for each excitation wavelength)
#' Missing values are interpolated within EEM data
#'
#' @description Missing EEM data can be interpolated. Usually it is the result of removing scatter. It is done along each excitation wavelength. This step is recommended if you aim for a PARAFAC analysis. Interpolation is done with hermitean interpolation polynomials using \code{\link[pracma]{pchip}}.
#'
#' @param data object of class eemlist with spectra containing missing values
#' @param cores specify number of cores for parallel computation
#'
#' @return object of class eemlist with interpoleted spectra.
#'
#' @seealso \code{\link[pracma]{pchip}}
#'
#' @references Elcoroaristizabal, S., Bro, R., García, J., Alonso, L. 2015. PARAFAC models of fluorescence data with scattering: A comparative study. Chemometrics and Intelligent Laboratory Systems, 142, 124-130
#' \url{https://doi.org/10.1016/j.chemolab.2015.01.017}
#'
#' @import dplyr tidyr
#' @importFrom stats na.omit
#' @importFrom pracma pchip
#' @export
#'
#' @examples
#' \donttest{
#' data(eem_list)
#'
#' remove_scatter <- c()
#' remove_scatter["raman1"] = TRUE
#' remove_scatter["raman2"] = TRUE
#' remove_scatter["rayleigh1"] = TRUE
#' remove_scatter["rayleigh2"] = TRUE
#' remove_scatter_width = c(15,10,16,12)
#'
#' eem_list <- eem_rem_scat(eem_list,remove_scatter,remove_scatter_width)
#'
#' eem_list <- eem_interp(eem_list)
#' }
eem_interp <- function(data,cores = detectCores(logical = FALSE)/2){
  cl <- makeCluster(spec = cores, type = "PSOCK")
  registerDoParallel(cl)

  eem_list <- foreach(eem = data) %dopar% {
    eem$x[1,which(is.na(eem$x[1,]))] <- 0
    eem$x[nrow(eem$x),which(is.na(eem$x[nrow(eem$x),]))] <- 0
    eem$x <- eem$x %>% apply(2,function(col) pracma::pchip(xi=eem$em[!is.na(col)],yi=col %>% na.omit(),x=eem$em))
    eem
  }

  stopCluster(cl)

  class(eem_list) <- "eemlist"
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
eem_raman_normalisation2 <- function(data,blank="blank"){
  #norm <- ram_data
  if(blank == "blank"){
    res_list <- try(data %>% eem_raman_normalisation(),silent=TRUE)
    if (class(res_list) == "try-error") {
      stop("There was a problem with raman normalisation, please check the presence of a blank sample and the parameter blank!")
    }
  } else if(is.data.frame(blank) & class(data) == "eemlist"){
    #data <- eem_list
    #rownames(blank) <- rownames(blank) %>% make.names()
    res_list <- lapply(data,function(eem){
      #eem <- data[[3]]
      #blank <- ram_data
      rar <- blank[eem$sample %>% make.names(),]
      #cat("Raman area:", rar, "\n")
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
      #cat("Raman area:", blank, "\n")
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
#' @param cuvl length of cuvette of absorption measurment in cm. Either a number or a data frame. Row names of data frame have to be similar to sample names in data
#'
#' @return fluorescence data of class eemlist
#' @import eemR
#' @export
#'
#' @examples
#' data(eem_list)
#' data(abs_data)
#'
#' eem_ife_correction(eem_list,abs_data,5)
eem_ife_correction <- function(data,abs_data,cuvl){
  eem_list <- data %>% lapply(function(eem1){
    if(is.data.frame(cuvl)) cl <- cuvl[eemnam,] else cl <- cuvl
    #eem1 <- eem_list %>% .[eem_list %>% eem_names() == eemnam]
    #eem1 <- eem_list[[48]]
    eem1 <- list(eem1)
    #class(eem_list) <- "eemlist"
    nam <- eem1[[1]]$sample
    class(eem1) <- "eemlist"
    eem1 <- eem1 %>% eem_inner_filter_effect(absorbance=abs_data,pathlength=cl)
    eem1[[1]]$sample <- nam
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
#' data(eem_list)
#'
#' eem_list <- eem_smooth(eem_list,n=4)
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
#' @param remove_scatter named logical vector. The names of the vector must be out of "raman1", "raman2", "rayleigh1" and "rayleigh2". Set \code{TRUE} if certain scattering should be removed.
#' @param remove_scatter_width numeric vector containing width of scattering to remove. If there is only one element in this vector, each this is the width of each removed scattering. If there are 4 values, differnt widths are used ordered by "raman1", "raman2", "rayleigh1" and "rayleigh2".
#' @param interpolation logical, optionally states whether interpolation is done right away
#' @param cores optional, CPU cores to use for interpolation
#'
#' @return eemlist
#'
#' @import eemR dplyr tidyr
#' @export
#'
#' @examples
#' data(eem_list)
#'
#' remove_scatter <- c()
#' remove_scatter["raman1"] = TRUE
#' remove_scatter["raman2"] = TRUE
#' remove_scatter["rayleigh1"] = TRUE
#' remove_scatter["rayleigh2"] = TRUE
#' remove_scatter_width = c(15,10,16,12)
#'
#' eem_rem_scat(eem_list,remove_scatter,remove_scatter_width)
eem_rem_scat <- function(data,remove_scatter,remove_scatter_width = 10, interpolation = FALSE, cores = parallel::detectCores()/2)
{
  if(data %>% class != "eemlist") stop("first argument has to be a list of eem samples!")
  if(!all(remove_scatter %>% names() %in% c("raman1","raman2","rayleigh1","rayleigh2"))) stop("scatter bands you would like to remove are unknown!")
  if(!is.numeric(remove_scatter_width)) stop("removed scatter slot width has to be numeric!")
  if(length(remove_scatter_width) == 1) remove_scatter_width <- rep(remove_scatter_width,4)
  if(remove_scatter[["raman1"]]) data <- data %>% eem_remove_scattering(type="raman",order=1,width=remove_scatter_width[1])
  if(remove_scatter[["raman2"]]) data <- data %>% eem_remove_scattering(type="raman",order=2,width=remove_scatter_width[2])
  if(remove_scatter[["rayleigh1"]]) data <- data %>% eem_remove_scattering(type="rayleigh",order=1,width=remove_scatter_width[3])
  if(remove_scatter[["rayleigh1"]]) data <- data %>% eem_remove_scattering(type="rayleigh",order=2,width=remove_scatter_width[4])

  if(interpolation) data <- eem_interp(data, cores = cores)

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
#' data(blank)
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

## reads list of txt and csv files in directory and produces table with absorbance, sample name as column name
#' Reading absorbance data from txt and csv files.
#'
#' @param absorbance_path directory containing absorbance data files or path to single file. See details for format of absorbance data.
#' @param order logical, data is ordered according to wavelength
#' @param recursive read files recursive, include subfolders
#' @param ... additional arguments that are passed on to \code{\link[utils]{read.table}}.
#'
#' @details If absorbance_path is a directory, contained files that end on "csv" or "txt" are passed on to \code{read.table}. If the path goes to a file, this file is passed on. Tables can either contain data from one sample or from several samples in columns. The column header containig the wavelength must be either "wavelength" or "Wavelength". A multi-sample file must have sample names as column names. A single-sample file can have sample name as column name or sample name as file name and "Abs." as column name. All tables are combined to one with one wavelength column and one column for each sample containing the absorbance data.
#'
#' @return A data frame containing absorbance data
#'
#' @seealso \code{\link[utils]{read.table}}
#'
#' @import dplyr tidyr
#' @importFrom utils read.table
#' @importFrom stringr str_replace_all str_replace
#'
#' @export
#'
#' @examples
#' absorbance_path <- system.file("extdata", "absorbance_eemR", package = "staRdom")
#' absorbance_read(absorbance_path,sep = " ", dec = ".")
absorbance_read <- function(absorbance_path,order=TRUE,recursive=TRUE,...){
  if(dir.exists(absorbance_path)){
    #absorbance_path <- absorbance_dir
    abs_data <- list.files(absorbance_path, full.names = TRUE, recursive = recursive, no.. = TRUE, include.dirs = FALSE, pattern = "*.txt|*.csv", ignore.case = TRUE)
    abs_data <- abs_data[!file.info(abs_data)$isdir]
    #abs_data <- dir(absorbance_path, pattern=".txt|.csv") %>%
    #  paste0(absorbance_path,"/",.)
  } else if(file.exists(absorbance_path)){
    abs_data <- absorbance_path
  } else stop("Absorbance data was not found!")
  abs_data <- abs_data %>%
    lapply(function(tab) {
      #print(tab)
      #tab <- abs_data[1]
      #table <- read.table(tab,header=TRUE,row.names=NULL,sep = ",", dec = ".",skip=1) %>%
      table <- read.table(tab,header=TRUE,row.names=NULL,...) %>%
        rename_(.dots=setNames(names(.), names(.) %>%
                                 stringr::str_replace_all("^.*Wavelength.*$|^.*wavelength.*$","wavelength") %>%
                                 stringr::str_replace_all("Abs.",tab %>% basename() %>%
                                                            stringr::str_replace_all(".txt$|.csv$","") #%>%
                                                          #                                                            stringr::str_replace("(^[:digit:]){1}","X\\1")
                                 ) #%>%
                               #                                 make.names()
        ))
    })
  if(abs_data %>% length() == 1) abs_data <- abs_data[[1]] %>% as.data.frame() else abs_data <- abs_data %>%     list_join(by="wavelength")
  if(order) abs_data <- abs_data %>% arrange(wavelength)
  abs_data
}

#' Calculating slopes and slope ratios of a data frame of absorbance data.
#'
#' @param abs_data data frame containing absorbance data.
#' @param cuvle length of used cuvette in cm
#' @param limits list with vectors containig upper and lower bounds of wavelengeth ranges to be fitted
#' @param l_ref list with reference wavelengths, same length as limits
#' @param S logical, include slope parameter in the table
#' @param lref logical, include reference wavelength in the table
#' @param p logical, include ps of the coefficients in the table
#' @param model logical, include complete model in data frame
#' @param cores number of cores to be used for parallel processing
#'
#' @details The absorbance data is a data frame with the first column called "wavelength" containg the wavelength. Each other column contains the data from one sample. You can use \link{absorbance_read} to read in appropriate data.
#'
#' The following spectral parameters are calculated:
#' \itemize{
#'    \item $S_{275-295}$ slope between 275 and 295 nm calculated with nonlinear regression
#'    \item $S_{350-400}$ slope between 350 and 400 nm calculated with nonlinear regression
#'    \item $S_{300-700}$ slope between 275 and 295 nm calculated with nonlinear regression
#'    \item SR slope ratio, calculated by $S_{275-295}$/$S_{350-400}$
#'    \item E2:E3 ratio $a_{250}$/$a_{365}$
#'    \item E4:E6 ratio $a_{465}$/$a_{665}$
#'    \item $a_{254}$ absorbance at 254 nm
#'    \item $a_{300}$ absorbance at 300 nm
#'    }
#' Depending on available wavelength range, values might be NA.
#' Additionally other wavelength limits can be defined. The slope ratio might fail in this case.
#' For further details please refer to Helm et al. (2008).
#'
#' @return A data frame containing the adsorption slopes and slope ratios in column, one line for each sample.
#'
#' @references Helms, J., Kieber, D., Mopper, K. 2008. Absorption spectral slopes and slope ratios as indicators of molecular weight, source, and photobleaching of chromophoric dissolved organic matter. Limnol. Oceanogr., 53(3), 955–969
#' \url{http://onlinelibrary.wiley.com/doi/10.4319/lo.2008.53.3.0955/pdf}
#'
#' @import dplyr
#' @import doParallel
#' @import foreach
#' @importFrom stats coef
#' @import drc
#' @importFrom tibble rownames_to_column
#' @export
#'
#' @examples
#' \donttest{
#' data(abs_data)
#' abs_parms(abs_data[,1:5],5)
#' abs_parms(abs_data[,1:5],5,l_ref=list(NA,NA,NA), lref=TRUE) # fit lref as well
#' }
abs_parms <- function(abs_data,cuvle,limits=list(c(275,295),c(350,400),c(300,700)),l_ref=list(350,350,350),S=TRUE,lref=FALSE,p=FALSE,model=FALSE,cores = detectCores()/2){
  #samp <- names(abs_data)[2]
  #library(doParallel)
  #cores=2
  #cuvle = cuvl
  cl <- makeCluster(spec = cores, type = "PSOCK")
  registerDoParallel(cl)

  res <- foreach(samp = abs_data %>% names() %>% .[. != "wavelength"]) %dopar% {
    r <- lapply(seq(1,length(limits)),function(n){
      wllim <- limits[[n]]
      ref <- l_ref[[n]]
      do(abs_data,m1 = abs_fit_slope(.$wavelength,.[[samp]]/cuvle*100,wllim,ref))
    }) %>%
      bind_cols(as_tibble(samp),.) %>%
      setNames(c("sample",lapply(limits,paste0,collapse="_") %>% unlist() %>% paste0("model",.)))
  } %>%
    bind_rows()

  stopCluster(cl)

  if(S){
    res <- res %>%
      select(contains("model")) %>%
      rowwise() %>%
      mutate_at(.vars = vars(contains("model")), .funs = funs(ifelse(class(.) == "drc",coef(.)["S:(Intercept)"],NA))) %>%
      setNames(names(.) %>% stringr::str_replace("model","S")) %>%
      mutate(SR = S275_295/S350_400) %>%
      bind_cols(res,.)
  }

  if(S & p){
    res <- res %>%
      select(contains("model")) %>%
      rowwise() %>%
      mutate_at(.vars = vars(contains("model")), .funs = funs(ifelse(class(.) == "drc",summary(.)$coefficients["S:(Intercept)","p-value"],NA))) %>%
      setNames(names(.) %>% stringr::str_replace("model","p_S")) %>%
      bind_cols(res,.)
  }

  if(lref){
    res <- res %>%
      select(contains("model")) %>%
      rowwise() %>%
      mutate_at(.vars = vars(contains("model")), .funs = funs(ifelse(class(.) == "drc",ifelse(is.null(attr(., "lref")),coef(.)["lref:(Intercept)"],attr(., "lref")),NA))) %>%
      setNames(names(.) %>% stringr::str_replace("model","lref")) %>%
      bind_cols(res,.)
  }

  if(lref & p){
    res <- res %>%
      select(contains("model")) %>%
      rowwise() %>%
      mutate_at(.vars = vars(contains("model")), .funs = funs(ifelse(class(.) == "drc",ifelse(is.null(attr(., "lref")),summary(.)$coefficients["lref:(Intercept)","p-value"],NA),NA))) %>%
      setNames(names(.) %>% stringr::str_replace("model","p_lref")) %>%
      bind_cols(res,.)
  }

  if(!model){
    res <- res %>%
      select(-contains("model"))
  }


  res <- abs_data %>%
    filter(wavelength == 250 | wavelength == 254 | wavelength == 300 | wavelength == 365 | wavelength == 465 | wavelength == 665) %>%
    t() %>%
    data.frame() %>%
    setNames(paste0("a",.[1,])) %>%
    tibble::rownames_to_column(var="sample") %>%
    filter(sample != "wavelength") %>%
    mutate(E2_E3 = a250/a365, E4_E6 = a465/a665) %>%
    select(-a250,-a365,-a465,-a665) %>%
    full_join(res,by="sample")

  res
}

#' Fit absorbance data to exponential curve. \code{\link[drc]{drm}} is used for the fitting process.
#'
#' @param wl vector containing wavelengths
#' @param abs vector containing absorption in m^-1
#' @param lim vector containing lower and upper limits for wavelengths to use
#' @param l_ref numerical. reference wavelength, default is 350, if set to NA l_ref is fitted
#' @param control control parameters for drm, see \code{\link[drc]{drmc}}
#' @param ... parameters that are passed on to drm
#'
#' @return numeric exponential slope coefficient
#'
#' @seealso \code{\link[drc]{drm}}
#'
#' @importFrom drc drm drmc
#' @importFrom stats approx
#' @export
#'
#' @examples
#' data(abs_data)
#' abs_fit_slope(abs_data$wavelength,abs_data$sample1,lim=c(350,400),l_ref=350)
abs_fit_slope <- function(wl,abs,lim,l_ref = 350,control = drmc(errorm = FALSE, noMessage = TRUE),...){
  #library(minpack.lm)
  #library(nls2)
  if(!all(is.na(wl)) & !all(is.na(abs))){
    #abs <- abs_data[[3]]
    #wl <- abs_data[[1]]
    #lim <- limits[[1]]
    abs <- 2.303*abs
    al <- abs[wl > lim[1] & wl < lim[2]]
    l <- wl[wl > lim[1] & wl < lim[2]]
    dal <- diff(abs)/diff(wl)
    dal <- c(dal[1],dal)
    #lref <- l
    #library(drc)
    if(is.na(l_ref)){
    abs_curve_lr <- function(x, parm,al1=abs,ll=wl){approx(y=al1,x=ll,xout=parm[,1])[[2]]*exp(-parm[,2]*(x-parm[,1]))} ## parameters: Q,b,Q0
    abs_derivx_lr <- function(x, parm,al1=abs,ll=wl){approx(y=al1,x=ll,xout=parm[,1])[[2]]*exp(-parm[,2]*(x-parm[,1]))*(-parm[,2])} ## parameters: Q,b,Q0
    abs_deriv1_lr <- function(x, parm,al1=abs,ll=wl,d=dal){c(approx(y=d,x=ll,xout=parm[,1])[[2]]*exp(-parm[,2]*(x-parm[,1]))*(-parm[,2]) + approx(y=al1,x=ll,xout=parm[,1])[[2]]*exp(-parm[,2]*(x-parm[,1]))*(parm[,2]),approx(y=al1,x=ll,xout=parm[,1])[[2]]*exp(-parm[,2]*(x-parm[,1]))*(parm[,1]-x))} ## parameters: Q,b,Q0
    abs_ssfct_lr <- function(dframe){c(350,0.01)}
      lowerl = c(min(wl),0)
      upperl = c(max(wl),Inf)
    res <- try(suppressWarnings(drm(formula = al ~ l, data = data.frame(al,l),
                   fct = list(fct = abs_curve_lr, ssfct = abs_ssfct_lr, names = c("lref","S"), deriv1 = abs_deriv1_lr, derivx = abs_derivx_lr),
                   robust = "lts", start=c(350, 0.01), lowerl = lowerl, upperl = upperl, control = control, ...)),
               silent=TRUE)
    } else {
      abs_curve <- function(x, parm,al1=abs,ll=wl,lref = l_ref){approx(y=al1,x=ll,xout=lref)[[2]]*exp(-parm[,1]*(x-lref))} ## parameters: Q,b,Q0
      abs_derivx <- function(x, parm,al1=abs,ll=wl,lref = l_ref){approx(y=al1,x=ll,xout=lref)[[2]]*exp(-parm[,1]*(x-lref))*(-parm[,2])} ## parameters: Q,b,Q0
      abs_deriv1 <- function(x, parm,al1=abs,ll=wl,d=dal,lref = l_ref){c(approx(y=al1,x=ll,xout=lref)[[2]]*exp(-parm[,1]*(x-lref))*(lref-x))} ## parameters: Q,b,Q0
      abs_ssfct <- function(dframe){c(0.01)}
      res <- try(suppressWarnings(drm(formula = al ~ l, data = data.frame(al,l),
                                      fct = list(fct = abs_curve, ssfct = abs_ssfct, names = c("S"), deriv1 = abs_deriv1, derivx = abs_derivx),
                                      robust = "lts", start=c(0.01), lowerl = c(0), upperl = c(Inf), control = control, ...)),
                 silent=TRUE)
      attr(res, "lref") <- l_ref
    }
    #if(class(res) == "drc") res[[1,2]]$coefficients <- setNames(res[[1,2]]$coefficients,res[[1,2]]$parNames[[2]])
  } else {
    warning("Data columns are empty.")
    res <- NA
  }
  res
}

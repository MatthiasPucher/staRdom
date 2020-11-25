## reads list of txt and csv files in directory and produces table with absorbance, sample name as column name
#' Reading absorbance data from txt and csv files.
#'
#' @param absorbance_path directory containing absorbance data files or path to single file. See details for format of absorbance data.
#' @param order logical, data is ordered according to wavelength
#' @param recursive read files recursive, include subfolders
#' @param dec optional, either you set a decimal separator or the table is tested for . and ,
#' @param sep optional, either you set a field separator or it is tried to be determined automatically
#' @param verbose logical, provide more information
#' @param cores number of CPU cores to be used simultanuously
#' @param ... additional arguments that are passed on to \code{\link[data.table]{fread}}.
#'
#' @details If absorbance_path is a directory, contained files that end on "csv" or "txt" are passed on to \code{read.table}. If the path is a file, this file is read. Tables can either contain data from one sample or from several samples in columns. The first column is considered the wavelength column. A multi-sample file must have sample names as column names. All tables are combined to one with one wavelength column and one column for each sample containing the absorbance data.
#'     Column and decimal separators are guessed from the supplied data. In some cases, this can lead to strange results. Plaese set `sep` and `dec` manually if you encounter any problems.
#'
#' @return A data frame containing absorbance data. An attribute "location" contains the filenames where each sample was taken from.
#'
#' @seealso \code{\link[data.table]{fread}}
#'
#' @import dplyr tidyr parallel
#' @importFrom utils read.table
#' @import stringr
#' @importFrom data.table fread
#'
#' @export
#'
#' @examples
#' absorbance_path <- system.file("extdata", "absorbance", package = "staRdom")
#' absorbance <- absorbance_read(absorbance_path, verbose = TRUE, cores = 2)
absorbance_read <- function(absorbance_path, order = TRUE, recursive = TRUE, dec = NULL, sep = NULL, verbose = FALSE, cores = parallel::detectCores(logical = FALSE), ...){
  #absorbance_path = absorbance_dir
  if(dir.exists(absorbance_path)){
    abs_data <- list.files(absorbance_path, full.names = TRUE, recursive = recursive, no.. = TRUE, include.dirs = FALSE, pattern = "*.txt|*.csv", ignore.case = TRUE)
    abs_data <- abs_data[!file.info(abs_data)$isdir]
  } else if(file.exists(absorbance_path)){
    abs_data <- absorbance_path
  } else stop("Absorbance data was not found!")

  cl <- makeCluster(min(cores, length(abs_data)), type="PSOCK")
  clusterExport(cl, c("dec","sep","verbose"), envir=environment())
  clusterEvalQ(cl, require(dplyr))
  clusterEvalQ(cl, require(stringr))

  #abs_data <- lapply(abs_data, function(tab) {
  abs_data <- parLapply(cl,abs_data, function(tab){
    #if(verbose) warning("processing",tab,fill=TRUE)
    tryCatch({
    rawdata <- readLines(tab)
    data <- rawdata %>%
      sapply(str_remove,pattern="([^0-9]*$)")
    first_number <- min(which((substr(data,1,1) %>% grepl("[0-9]",.))))
    last_number <- max(which((substr(data,1,1) %>% grepl("[0-9]",.))))
    if(is.null(sep) | is.null(dec)){
      nsepdec <- data[first_number] %>% str_extract_all("[^-0-9eE]") %>% unlist()
      example_number <- data[first_number] %>% str_extract("([-]?[0-9]+[.,]?[0-9]+[eE]?[-0-9]+)$")
      if(is.null(dec) & length(nsepdec) > 1) dec <- example_number %>% str_replace("([-0-9eE]+)([.,]?)([-0-9eE]*)","\\2")
      if(is.null(sep)) sep <- gsub(pattern = dec, replacement = "", x = data[first_number], fixed = TRUE) %>%
        str_extract(paste0("[^-0-9eE",dec,"]"))
      if(verbose) warning("processing",tab,": using",sep,"as field separator and",dec,"as decimal separator.", fill=TRUE)
    }
    data <- str_split(data,sep)
    table <- data[(first_number):last_number] %>%
      unlist() %>%
      matrix(ncol = length(data[[first_number]]), byrow = TRUE) %>%
      data.frame(stringsAsFactors = FALSE) %>%
      mutate_all(gsub,pattern=ifelse(dec != "",dec,"."),replacement=".",fixed=TRUE)
    # if(any(as.numeric(table[1,]) %>% is.na())){
    #   table <- table %>%
    #   setNames(.[1,])
    # }
    table <- table %>%
      mutate_all(as.numeric)
    attr(table,"location") <- rep(tab,ncol(table) - 1)
    if(ncol(table) == 2) {samples <- tab %>%
      basename() %>%
      str_replace_all(regex(".txt$|.csv$", ignore_case = TRUE),"")
    } else {
      samples <- rawdata[[1]] %>%
        str_split(sep) %>%
        unlist() %>%
        matrix(ncol = length(.), byrow = TRUE) %>%
        data.frame(stringsAsFactors = FALSE) %>%
        .[-1]
    }
    table <- table %>%
      setNames(c("wavelength",samples)
      )},
    error = function(err){
      stop("Error while reading ",tab,": ",err)
      #warning("Error while reading ",tab,": ",err,"\nFile ",tab," was not included in the absorbance data!")
    })
  })
  stopCluster(cl)

  locations <- lapply(abs_data,function(tab){
    attr(tab,"location")
  }) %>%
    unlist()
  if(length(abs_data) == 1) abs_data <- abs_data[[1]] %>% as.data.frame() else abs_data <- abs_data %>% list_join(by="wavelength")
  if(order) abs_data <- abs_data %>% arrange(wavelength)
  attr(abs_data,"location") <- locations
  abs_data
}

#' Calculating slopes and slope ratios of a data frame of absorbance data.
#'
#' @param abs_data data frame containing absorbance data.
#' @param cuvle cuvette (path) length in cm, ignored if unit is absorption
#' @param unit unit of absorbance data: if "absorbance", absorbance data is multiplied by log(10) = 2.303 for slope calculations
#' @param add_as additionally to a254 and a300, absorbance at certain wavelengths can be added to the table
#' @param limits list with vectors containig upper and lower bounds of wavelengeth ranges to be fitted
#' @param l_ref list with reference wavelengths, same length as limits
#' @param S logical, include slope indices in the table
#' @param lref logical, include reference wavelength in the table
#' @param p logical, include ps of the coefficients in the table
#' @param model logical, include complete model in data frame
#' @param Sint logical, wether the spectral curve is calculated interval-wise (\code{\link[cdom]{cdom_spectral_curve}})
#' @param interval passed on to \code{\link[cdom]{cdom_spectral_curve}}
#' @param r2threshold passed on to \code{\link[cdom]{cdom_spectral_curve}}
#' @param cores number of cores to be used for parallel processing
#' @param verbose logical, additional information is provided
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
#' \url{https://aslopubs.onlinelibrary.wiley.com/doi/pdf/10.4319/lo.2008.53.3.0955}
#'
#' @import dplyr
#' @import doParallel
#' @import foreach
#' @importFrom stats coef
#' @import drc
#' @importFrom tibble rownames_to_column
#' @importFrom cdom cdom_spectral_curve
#' @export
#'
#' @examples
#' \donttest{
#' data(absorbance)
#'
#' a1 <- abs_parms(absorbance, cuvle = 5, verbose = TRUE, cores = 2)
#' a2 <- abs_parms(absorbance, cuvle = 5,l_ref=list(NA,NA,NA), lref=TRUE, cores = 2) # fit lref as well
#' }
abs_parms <- function(abs_data, cuvle = NULL, unit = c("absorbance", "absorption"), add_as = NULL, limits = list(c(275, 295), c(350, 400), c(300, 700)), l_ref = list(275, 350, 300), S = TRUE, lref = FALSE, p = FALSE, model = FALSE, Sint = FALSE, interval = 21, r2threshold = 0.8, cores = parallel::detectCores(logical = FALSE), verbose = FALSE){
  if(!unit[1] %in% c("absorbance","absorption")) stop("Unit must be either 'absorbance' or 'absorption'!")
  if(unit[1] == "absorbance" & !is.numeric(cuvle)) stop("Please specify a valid cuvette length!")
  if(S | lref | p | model){
    cl <- makeCluster(spec = min(cores,ncol(abs_data)), type = "PSOCK")
    doParallel::registerDoParallel(cl)

    if(verbose){
      cat("calculating slopes of",ncol(abs_data)-1,unit[1],"spectra...", fill=TRUE)
    }

    res <- foreach(i = 1:(ncol(abs_data)-1)) %dopar% { #, .options.snow = opts
      samp <- abs_data %>% names() %>% .[. != "wavelength"] %>% .[i]
      r <- lapply(seq(1,length(limits)),function(n){
        wllim <- limits[[n]]
        ref <- l_ref[[n]]
        do(abs_data,m1 = abs_fit_slope(.$wavelength,.[[samp]]*ifelse(unit[1] == "absorbance",100/cuvle*log(10),1),wllim,ref))
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
  }

  if(Sint){
    cl <- parallel::makeCluster(spec = min(cores,ncol(abs_data)), type = "PSOCK")
    doParallel::registerDoParallel(cl)

    if(verbose){
      cat("calculating spectral curves interval-wise",ncol(abs_data)-1,"absorbance spectra...", fill=TRUE)
    }

    resSint <- foreach(i = 1:(ncol(abs_data)-1)) %dopar% { # , .options.snow = opts
      samp <- abs_data %>% names() %>% .[. != "wavelength"] %>% .[i]
      r <- do(abs_data,Sint = cdom::cdom_spectral_curve(.$wavelength,.[[samp]]*ifelse(unit[1] == "absorbance",100/cuvle*log(10),1),interval = interval, r2threshold = r2threshold)) %>%
        bind_cols(as_tibble(samp),.) %>%
        setNames(c("sample","Sint"))
    } %>%
      bind_rows()

    stopCluster(cl)

    if(S | lref | p | model){
      res <- full_join(res,resSint,by="sample")
    } else {
      res <- resSint
    }

    res
  }

  as_vec <- c(250,254,300,365,465,665,add_as) %>% unique()

  as <- abs_data[-1] %>%
    apply(2,function(col) approx(x=abs_data$wavelength,y=col*ifelse(unit[1] == "absorbance",100/cuvle*log(10),1),xout=as_vec)) %>%
    lapply(`[[`,"y") %>%
    bind_cols(wavelength=as_vec,.)

  res <- as %>%
    t() %>%
    data.frame() %>%
    setNames(paste0("a",.[1,])) %>%
    tibble::rownames_to_column(var="sample") %>%
    filter(sample != "wavelength") %>%
    mutate(E2_E3 = a250/a365, E4_E6 = a465/a665) %>%
    select(-one_of(paste0("a",c(250,365,465,665) %>% .[!(. %in% add_as)]))) %>%
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
#' data(absorbance)
#' abs_fit_slope(absorbance$wavelength,absorbance$sample1,lim=c(350,400),l_ref=350)
abs_fit_slope <- function(wl,abs,lim,l_ref = 350,control = drmc(errorm = FALSE, noMessage = TRUE),...){
  if(!all(is.na(wl)) & !all(is.na(abs))){
    if(lim[1] < min(wl) | lim[2] > max(wl)){
      warning("The absorbance wavelength is between ",min(wl)," and ",max(wl)," so slope calculations between ",lim[1]," and ",lim[2]," are not possible.")
      res <- NA
      attr(res, "lref") <- NA
    } else {
      al <- abs[wl > lim[1] & wl < lim[2]]
      l <- wl[wl > lim[1] & wl < lim[2]]
      dal <- diff(abs)/diff(wl)
      dal <- c(dal[1],dal)
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
    }
  } else {
    warning("Data columns are empty.")
    res <- NA
    attr(res, "lref") <- NA
  }
  res
}


#' Baseline correction for absorbance data
#'
#' @param abs_data data.frame containing samples in columns, the column containing wavelengths must be named "wavelength"
#' @param wlrange range of wavelengths that should be used for correction, absorbance mean in that range is subtracted from each value (sample-wise)
#'
#' @return data.frame
#' @export
#'
#' @import dplyr
#'
#' @examples
#' data(absorbance)
#' abs_data_cor <- abs_blcor(absorbance)
#'
#' abs_data_cor1 <- abs_blcor(absorbance[1:2])
#'
abs_blcor <- function(abs_data, wlrange = c(680,700)){
  ad <- abs_data %>%
    .[.$wavelength >= wlrange[1] & .$wavelength <= wlrange[2],] %>%
    apply(2,mean,na.rm=TRUE)
  ad2 <- sweep(data.matrix(abs_data),2,ad) %>%
    as.data.frame() %>%
    select(-1) %>%
    bind_cols(wavelength = abs_data$wavelength,.)
}

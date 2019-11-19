#' Runs a PARAFAC analysis on EEM data
#'
#' @description One or more PARAFAC models can be calculated depending on the number of components. The idea is to compare the different models to get the most suitable. B-mode is emmission wavelengths, C-mode is excitation wavelengths and, A-mode is the loadings of the samples. The calculation is done with \code{\link[multiway]{parafac}}, please see details there.
#'
#' @param eem_list object of class \code{\link[eemR]{eem}}
#' @param comps vector containing the desired numbers of components. For each of these numbers one model is calculated
#' @param maxit maximum iterations for PARAFAC algorithm
#' @param normalise state whether EEM data should be normalised in advance
#' @param nstart number of random starts
#' @param cores number of parallel calculations (e.g. number of physical cores in CPU)
#' @param ctol Convergence tolerance (R^2 change)
#' @param strictly_converging calculate nstart converging models and take the best. Please see details!
#' @param const constraints of PARAFAC analysis. Default is non-negative ("nonneg"), alternatively smooth and non-negative ("smonon") might be interesting for an EEM analysis.
#' @param output Output the \code{"best"} solution (default) only or additionally add \code{"all"} \code{nstart} solutions to the model as an element named \code{"models"}.
#' @param verbose print infos
#' @param ... additional parameters that are passed on to \code{\link[multiway]{parafac}}
#'
#' @details PARAFAC models are created based on multiple random starts. In some cases, a model does not converge and the resulting model is then based on less than nstart converging models. In case you want to have nstart converging models, set strictly_converging TRUE. This calculates models stepwise until the desired number is reached but it takes more calculation time. Increasing the number of models from the beginning is much more time efficient.
#'
#' @return object of class parafac
#' @export
#'
#' @import multiway
#' @import parallel
#' @import dplyr
#' @import tidyr
#'
#' @seealso \code{\link[multiway]{parafac}}
#'
#' @examples
#' \donttest{
#' data(eem_list)
#'
#' dim_min <- 3 # minimum number of components
#' dim_max <- 7 # maximum number of components
#' nstart <- 25 # random starts for PARAFAC analysis, models built simulanuously, best selected
#' cores <- parallel::detectCores(logical=FALSE) # use all cores but do not use all threads
#' maxit = 2500
#' ctol <- 10^-7 # tolerance for parafac
#'
#' pfres_comps <- eem_parafac(eem_list, comps = seq(dim_min, dim_max),
#'     normalise = TRUE, maxit = maxit, nstart = nstart, ctol = ctol, cores = cores)
#'
#' pfres_comps2 <- eem_parafac(eem_list, comps = seq(dim_min, dim_max),
#'     normalise = TRUE, maxit = maxit, nstart = nstart, ctol = ctol, cores = cores, output = "all")
#' }
eem_parafac <- function(eem_list, comps, maxit = 2500, normalise = TRUE, const = c("nonneg","nonneg","nonneg"), nstart = 30, ctol = 10^-8, strictly_converging = FALSE, cores = parallel::detectCores(logical=FALSE), verbose = FALSE, output = "best",...){
  eem_array <- eem2array(eem_list)
  if(normalise){
    if(verbose) cat("EEM matrices are normalised!",fill=TRUE)
    eem_array <- eem_array %>% norm_array()
  }
  if(verbose) cat(paste0(cores," cores are used for the calculation."),fill=TRUE)
  res <- lapply(comps,function(comp){
    #comp <- 6
    if(verbose) cat(paste0("calculating ",comp," components model..."),fill=TRUE)
    cl <- NULL
    if(cores > 1){
      cl <- makeCluster(min(cores,nstart), type="PSOCK")
      clusterExport(cl, c("eem_array","comp","maxit","const","ctol","cores"), envir=environment())
      clusterEvalQ(cl, library(multiway))
    }
    if(strictly_converging){
      cpresult <- parafac_conv(eem_array, nfac = comp, const = const, maxit = maxit, parallel = (cores > 1), cl = cl, ctol = ctol, nstart = nstart, output = "all", verbose = verbose, ...)
    } else {
      cpresult <- parafac(eem_array, nfac = comp, const = const, maxit = maxit, parallel = (cores > 1), cl = cl, ctol = ctol, nstart = nstart, output = "all",...)#, ...
    }
    Rsqs <- lapply(cpresult,`[[`,"Rsq") %>% unlist()
    cpresult1 <- cpresult[[which.max(Rsqs)]]
    if(cores > 1){
      stopCluster(cl)
    }
    cflags <- lapply(cpresult,`[[`,"cflag") %>% unlist()
    converged <- sum(cflags == 0)/nstart
    if(converged <= 0.5){
      warning("Calculating the ",comp," component",ifelse(comp > 1, "s","")," model, only ",sum(cflags == 0)," out of ",nstart," models converged! You might want to increase the number of initialisations (nstart) or iterations (maxit).")
    }
    if(cpresult1$cflag != 0){
      warning("The PARAFAC model with ",comp," component",ifelse(comp > 1, "s", "")," did not converge! Increasing the number of initialisations (nstart) or iterations (maxit) might solve the problem.")
    }
    if(output == "all"){
      cpresult1$models <- lapply(cpresult,.trans_parafac, em = eem_list[[1]]$em, ex = eem_list[[1]]$ex, samples = eem_list %>% eem_names(), comp = comp, const = const, norm_factors = attr(eem_array,"norm_factors"))
    }
    cpresult1 <- .trans_parafac(cpresult1, em = eem_list[[1]]$em, ex = eem_list[[1]]$ex, samples = eem_list %>% eem_names(), comp = comp, const = const, norm_factors = attr(eem_array,"norm_factors"))
    cpresult1$converged <- converged
    cpresult1
  })
  mostattributes(res) <- attributes(eem_array)
  return(res)
}

#' Add data of a PARAFAC model derived from multiway from EEMs
#'
#' @param parafac parafac model
#' @param em emission wavelengths
#' @param ex excitation wavelengths
#' @param samples sample names
#' @param comp number of components
#' @param const constraints
#' @param norm_factors factors to invert normalisation
#'
#' @return parafac model
.trans_parafac <- function(parafac, em, ex, samples, comp, const, norm_factors){
  attr(parafac,"norm_factors") <- norm_factors
  rownames(parafac$B) <- em #eem_list[[1]]$em
  rownames(parafac$C) <- ex #eem_list[[1]]$ex
  rownames(parafac$A) <- samples #eem_list %>% eem_names()
  labComp <- paste("Comp.",1:comp,sep="")
  colnames(parafac$A) <- labComp
  colnames(parafac$B) <- labComp
  colnames(parafac$C) <- labComp
  # small issue with multiway: slightly negative values are possible despite using nonnegative constraints
  non <- grepl("no|unsmpn",const)
  if(non[1]) parafac$A[parafac$A < 0] <- 0
  if(non[2]) parafac$B[parafac$B < 0] <- 0
  if(non[3]) parafac$C[parafac$C < 0] <- 0
  parafac
}

#' Calculate a PARAFAC model similar to and using \code{\link[multiway]{parafac}}.
#'
#' @description Please refer to \code{\link[multiway]{parafac}} for input parameters and details. This wrapper function ensures `nstart` converging models are calculated. On the contrary, parafac calculates `nstart` models regardless if they are converging.
#'
#' @param X array
#' @param nstart number of converging models to calculate
#' @param verbose logical, whether more information is supplied
#' @param output Output the best solution (default) or output all nstart solutions.
#' @param cl cluster to be used for parallel processing
#' @param ... arguments passed on to \code{\link[multiway]{parafac}}
#'
#' @return either a parafac model or a list of parafac models
#'
#' @seealso \code{\link[multiway]{parafac}}
#'
#' @export
#'
#' @examples
#' \donttest{
#' # sorry, no example provided yet
#' }
parafac_conv <- function(X, nstart, verbose = FALSE, output = c("best", "all"), cl = NULL, ...){
  nmod <- 0
  ntot <- 0
  cpresult_all <- list()
  while(nmod < nstart & ntot <= 10 * nstart){
    pred_factor <- ifelse(ntot == 0, 1, ifelse(nmod == 0, 3, ntot/nmod/2))
    if(verbose) cat("Due to previous model calculations, pred_factor was set", pred_factor, fill = TRUE)
    nmiss <- ceiling((nstart - nmod) * pred_factor / 8) * 8
    if(verbose) cat("start run with",nmiss,"models...",fill = TRUE)
    if(!is.null(cl)){
      clusterExport(cl, c("nmiss"), envir = environment())
    }
    cpresult <- parafac(X, nstart = nmiss, output = "all", cl = cl, ...) #, ...
    cpresult_conv <- cpresult[lapply(cpresult,`[[`,"cflag") %>% unlist() == 0]
    cpresult_all <- c(cpresult_all,cpresult_conv)
    nmod <- length(cpresult_all)
    ntot <- ntot + nmiss
    if(verbose) cat(length(cpresult_conv),"models converged successfully in this run!", fill = TRUE)
    if(verbose) cat(length(cpresult_all),"models calculated!", fill = TRUE)
  }
  if(verbose) cat(nmod,"out of",ntot,"models converged!",fill=TRUE)
  if(output[1] == "best"){
    sses <- cpresult_all %>% lapply(`[[`,"SSE") %>% unlist()
    cpresult_all <- cpresult_all[[which.min(sses)]]
  } else if (output[1] == "all"){
    sses <- cpresult_all %>% lapply(`[[`,"SSE") %>% unlist()
    cpresult_all <- cpresult_all[[sses[sort(order(-sses)[1:nstart])]]]
  }
  if(ntot >= 10 * nstart) warning("Maximum number of starts reached without generating the desired number of valid models.")
  cpresult_all
}

#' Rescale B and C modes of PARAFAC model
#'
#' @description B and C modes (emission and excitation wavelengths) are rescaled to RMS of value newscale. This is compensated in A mode (sample loadings).
#'
#' @param pfmodel object of class parafac
#' @param newscale If (default) newscale = "Fmax", each component will be scaled so the maximum of each component is 1. It is also possible to set a desired root mean-square for each column of the rescaled mode. Can input a scalar or a vector with length equal to the number of factors for the given mode.
#'
#' @return object of class parafac
#' @export
#'
#' @seealso \code{\link[multiway]{rescale}}
#'
#' @examples
#' data(pf_models)
#'
#' new_pf <- eempf_rescaleBC(pf4[[1]])
eempf_rescaleBC <- function(pfmodel,newscale = "Fmax"){
  nf <- attr(pfmodel,"norm_factors")
  comp <- ncol(pfmodel$A)
  if(newscale == "Fmax"){
    Bmax <- pfmodel$B %>% abs() %>% matrixStats::colMaxs()
    Cmax <- pfmodel$C %>% abs() %>% matrixStats::colMaxs()
    pfmodel$B <- pfmodel$B %*% diag(1/Bmax, nrow=comp) %>%
      `colnames<-`(colnames(pfmodel$B))
    pfmodel$C <- pfmodel$C %*% diag(1/Cmax, nrow=comp) %>%
      `colnames<-`(colnames(pfmodel$B))
    pfmodel$A <- pfmodel$A %*% diag(Bmax, nrow=comp) %*% diag(Cmax, nrow=comp) %>%
      `colnames<-`(colnames(pfmodel$B))
  } else {
    pfmodel <- rescale(pfmodel,mode = "C", newscale = newscale, absorb = "A")
    pfmodel <- rescale(pfmodel,mode = "B", newscale = newscale, absorb = "A")
  }
  attr(pfmodel,"norm_factors") <- nf
  labComp <- paste("Comp.",1:comp,sep="")
  colnames(pfmodel$A) <- labComp
  colnames(pfmodel$B) <- labComp
  colnames(pfmodel$C) <- labComp
  pfmodel
}

#' Extract names from PARAFAC model components
#'
#' @param pfmodel parafac model
#'
#' @return vector of names or list of vecters of names
#' @export
#'
#' @examples
#' data(pf_models)
#' eempf_comp_names(pf4)
#'
#' eempf_comp_names(pf4) <- c("A","B","C","D","E","F","G")
#'
#' value <- list(c("A1","B1","C1","D","E","F","G"),
#' c("A2","B2","C","D","E","F","G"),
#' c("A3","B3","C","D","E","F","G"),
#' c("A4","B4","C","D","E","F","G"),
#' c("A5","B5","C","D","E","F","G5")
#' )
#'
#' eempf_comp_names(pf4) <- value
#' eempf_comp_names(pf4)
#'
#' ggeem(pf4[[1]])
#'
eempf_comp_names <- function(pfmodel){
  if(class(pfmodel) == "parafac") {
    colnames(pfmodel$A)
  }else if(class(pfmodel) == "list" & class(pfmodel[[1]]) == "parafac"){
    lapply(pfmodel, function(pfm) colnames(pfm$A))
  } else{
    stop("pfmodel is not a parafac model or a list of parafac models!")
  }
}

#' Set names of PARAFAC components
#'
#' @param pfmodel model of class parafac
#' @param value character vector containing the new names for the components
#'
#' @return parafac model
#'
#' @import dplyr
#'
#' @export
#'
#' @examples
#' data(pf_models)
#'
#' eempf_comp_names(pf4) <- c("A","B","C","D","E","F","G")
`eempf_comp_names<-` <- function(pfmodel, value){
  if(class(pfmodel) == "parafac") {
    colnames(pfmodel$A) <- value
    colnames(pfmodel$B) <- value
    colnames(pfmodel$C) <- value
    pfmodel %>% `class<-`("parafac")
  }else if(class(pfmodel) == "list" & class(pfmodel[[1]]) == "parafac"){
    if(!is.list(value) | (length(value) == 1 & length(pfmodel) > 1)) value <- lapply(1:length(pfmodel), function(x) value)
    lapply(1:length(pfmodel), function(pfn){
      colnames(pfmodel[[pfn]]$A) <- value[[pfn]][1:ncol(pfmodel[[pfn]]$A)]
      colnames(pfmodel[[pfn]]$B) <- value[[pfn]][1:ncol(pfmodel[[pfn]]$A)]
      colnames(pfmodel[[pfn]]$C) <- value[[pfn]][1:ncol(pfmodel[[pfn]]$A)]
      pfmodel[[pfn]] %>% `class<-`("parafac")
    })
  } else{
    stop("pfmodel is not a parafac model or a list of parafac models!")
  }
}

#' Data from an eemlist is transformed into an array
#'
#' @description Data matrices from EEM are combined to an array that is needed for a PARAFAC analysis.
#'
#' @param eem_list object of class eemlist
#'
#' @return object of class array
#' @export
#'
#' @examples
#' data(eem_list)
#'
#' eem2array(eem_list)
eem2array <- function(eem_list){
  eem_matrices <- eem_list %>%
    sapply("[", "x")
  dv <- lapply(eem_matrices,dim) %>% bind_cols()
  if(all(dv[1,1] %>% unlist() == dv[1,]) & all(dv[2,1] %>% unlist() == dv[2,])) dim_eem <- c(dv[,1] %>% unlist(), eem_matrices %>% length()) else dim_eem <- NA
  if(is.na(dim_eem[1])) stop("dimensions mismatch!")

  eem_array <- array(eem_matrices %>% unlist, dim=dim_eem)
  eem_array <- eem_array %>% aperm(perm = c(3,1,2), resize = TRUE)
  attr(eem_array,"em") <- eem_list[[1]]$em
  attr(eem_array,"ex") <- eem_list[[1]]$ex
  attr(eem_array,"samples") <- eem_list %>% sapply("[","sample") %>% unlist()
  attr(eem_array,"mdim") <- dim(eem_array)
  eem_array
}

#' Normalise 3-dimensional array in first and second dimension
#'
#' @param eem_array 3-dimensional array
#'
#' @return array
#' @export
#'
#' @import dplyr
#' @import tidyr
#' @importFrom stats sd
#'
#' @examples
#' data(eem_list)
#'
#' a <- eem2array(eem_list)
#' an <- norm_array(a)
norm_array <- function(eem_array){
  norm_factors <- lapply(1:(dim(eem_array)[1]),function(s) {
    sd(eem_array[s,,], na.rm=TRUE)
  }) %>% unlist()
  eem_array <- eem_array / norm_factors
  attr(eem_array,"norm_factors") <- norm_factors
  eem_array
}


#' Extract EEM matrix for single components determined in the PARAFAC analysis
#'
#' @description The components of a PARAFAC analysis are extracted as a data frame
#'
#' @param pfmodel object of class parafac
#' @param gather logical value whether excitation wavelengths are a column, otherwise excitation wavelengths are column names
#'
#' @return a list of class data frames
#' @export
#'
#' @importFrom tibble rownames_to_column
#' @import dplyr
#' @import tidyr
#'
#' @examples
#' data(pf_models)
#'
#' eempf_comp_mat(pf4[[1]])
eempf_comp_mat <- function(pfmodel,gather=TRUE){
  mat <- lapply(seq(1:ncol(pfmodel$A)), function(comp){
    m <- matrix(pfmodel$B[,comp]) %*% t(matrix(pfmodel$C[,comp])) %>%
      data.frame()
    colnames(m) <- pfmodel$C %>% rownames()
    rownames(m) <- pfmodel$B %>% rownames()
    if(gather==TRUE) m <- m %>% tibble::rownames_to_column("em") %>% gather(ex,value,-em)
    m
  })
  names(mat) <- colnames(pfmodel$B)
  mat
}

#' Calculate the leverage of each emission and excitation wavelength and each sample from a single PARAFAC model
#'
#' @param pfmodel object of class parafac
#'
#' @return list of 3 named vectors (emission, excitation wavelengths and samples)
#' @export
#'
#' @importFrom pracma pinv
#' @importFrom stats setNames
#' @import dplyr
#'
#' @examples
#' data(pf_models)
#'
#' eempf_leverage(pf4[[1]])
eempf_leverage <- function(pfmodel){
  cpl <- lapply(pfmodel[c("A","B","C")], function(M) diag(M %*% pinv(t(M) %*% M) %*% t(M)))
  names <- list(rownames(pfmodel$A),rownames(pfmodel$B),rownames(pfmodel$C))
  cpl <- lapply(1:(cpl %>% length()),function(i){ cpl[[i]] %>% setNames(names[[i]])}) %>%
    setNames(c("A","B","C"))
}

#' Calculate the leverage of each emission and excitation wavelength and each sample from a list of PARAFAC models
#'
#' @param pfres_comps object of class parafac
#' @param ecdf logical, transforme leverages to according empirical quantiles (\code{\link[stats]{ecdf}})
#' @param stats logical, whether means and standard deviations are calculated from leverages
#'
#' @return data frame containing leverages of wavelengths and samples for each model
#' @export
#'
#' @import dplyr
#' @importFrom matrixStats rowSds
#' @importFrom stats ecdf
#'
#' @examples
#' data(pf_models)
#'
#' eempf_mleverage(pf3)
eempf_mleverage <- function(pfres_comps,ecdf = FALSE, stats = FALSE){
  cpls <- lapply(pfres_comps,eempf_leverage) %>%
    lapply(unlist) %>%
    lapply(function(ll) data.frame(parameter=names(ll),value=ll)) %>%
    list_join(by = "parameter") %>%
    `colnames<-`(c("parameter",paste0("comps",lapply(pfres_comps,function(cpout){ncol(cpout$A)}) %>% unlist())))
  if(ecdf){
    cpls <- cpls %>%
      mutate_if(is.numeric,function(col) ecdf(col)(col))
  }
  if(stats){
    cpls <- cpls %>%
      mutate(mean = rowMeans(select(., -parameter)), stdev=rowSds(select(., -parameter) %>% as.matrix()))
  }
  cpls
}


#' Combine leverages into one data frame and add optional labels.
#'
#' @param cpl leverage, outpout from \code{\link[staRdom]{eempf_leverage}}
#' @param qlabel optional, quantile of which labels are shown (1 = all, 0 = no labels)
#'
#' @return data frame
#' @export
#'
#' @importFrom tibble rownames_to_column
#' @import dplyr
#' @importFrom stats setNames
#' @importFrom stats quantile
#'
#' @examples
#' data(pf_models)
#'
#' leverage <- eempf_leverage(pf4[[1]])
#' lev_data <- eempf_leverage_data(leverage)
eempf_leverage_data <- function(cpl,qlabel=0.1){
  cpl <- cpl %>%
    lapply(. %>% data.frame() %>% rownames_to_column("x") %>% setNames(c("x","leverage")))
  mode_name <- c("sample","em","ex")
  cpl <- lapply(1:3,function(i){
    M <- cpl[[i]] %>% mutate(mode = mode_name[i])
  }) %>%
    bind_rows()
  pl <- cpl %>%
    group_by(mode) %>%
    mutate(q = quantile(leverage,1 - qlabel)) %>%
    mutate(label = ifelse(leverage > q,x,NA)) %>%
    ungroup()
}


#' Compensate for normalisation in C-modes
#'
#' @description Factors used for normalisation are saved separately in the PARAFAC models. With this function, the normalisation factors are combined with the A-modes of the model and removed as a separate vector. This means former normalisation is accounted for in the amount of each component in each sample. If no normalisation was done, the original model is returned without warning.
#'
#' @param pfmodel object of class parafac
#'
#' @return object of class parafac
#' @export
#'
#' @examples
#' data(pf_models)
#'
#' pf4[[1]] <- norm2A(pf4[[1]])
norm2A <- function(pfmodel){
  if(!is.null(attr(pfmodel,"norm"))){
    pfmodel$A <- pfmodel$A * attr(pfmodel,"norm_factors")
    attr(pfmodel,"norm_factors") <- NULL
  }
  pfmodel
}


#' Calculating correlations between the component loadings in all samples (C-Modes).
#'
#' @param pfmodel results from a PARAFAC analysis, class parafac
#' @param normalisation logical, whether normalisation is undone or not
#' @param method method of correlation, passed to \code{\link[stats]{cor}}
#' @param ... passed on to \code{\link[stats]{cor}}
#'
#' @return matrix
#' @export
#'
#' @importFrom stats cor
#'
#' @examples
#' data(pf_models)
#' eempf_cortable(pf4[[1]])
eempf_cortable <- function(pfmodel,normalisation = FALSE, method="pearson",...){
  if(normalisation) pfmodel <- norm2A(pfmodel)
  pfmodel %>%
    .$A %>%
    cor(method=method,...)
}

#' Extract data from emission and excitation wavelengths of the components of a PARAFAC model (scaled B- and C-modes)
#'
#' @description Data for each wavelengths is returned. For each component the lines intersecting at the component maxima are returned.
#'
#' @param pfmodel object of class parafac
#'
#' @return data frame
#' @export
#'
#' @import dplyr
#' @import tidyr
#'
#' @examples
#' data(pf_models)
#'
#' ml <- maxlines(pf4[[1]])
maxlines <- function(pfmodel){
  maxl <- lapply(colnames(pfmodel$C),function(comp){
    em = (pfmodel$B[,comp]*max(pfmodel$C[,comp])) %>%
      data.frame(e="em", wavelength = as.numeric(names(.)),value=.)
    ex = (pfmodel$C[,comp]*max(pfmodel$B[,comp])) %>%
      data.frame(e="ex", wavelength = as.numeric(names(.)),value=.)
    res <- bind_rows(em,ex) %>%
      setNames(c("e","wavelength",comp))
  }) %>%
    list_join(by=c("e","wavelength"))
}

#' Calculate residuals of EEM data according to a certain model
#'
#' @param pfmodel PARAFAC model of class parafac
#' @param eem_list eemlist containing EEM data
#' @param select character vector containing the names of the desired samples
#' @param cores number of cores to use for parallel processing
#'
#' @return data frame with EEM residuals
#' @export
#'
#' @import eemR
#' @import dplyr
#' @import tidyr
#' @import parallel
#'
#' @examples
#' \donttest{
#' data(eem_list)
#' data(pf_models)
#'
#' eempf_residuals(pf4[[1]],eem_list)
#' }
eempf_residuals <- function(pfmodel,eem_list,select=NULL, cores = parallel::detectCores(logical = FALSE)/2){
  pfmodel <- norm2A(pfmodel)
  if(!is.null(select)){
    eem_list <- eem_extract(eem_list,sample = select ,keep=TRUE,verbose = FALSE)
  }
  if(!all(eem_names(eem_list) %in% rownames(pfmodel$A)) | length(eem_list) == 0){
    pfmodel <- A_missing(eem_list,pfmodel,cores=cores)
  }
  what <- which(rownames(pfmodel$A) %in% (eem_list %>% eem_names()))
  pfmodel$A <- as.data.frame(pfmodel$A)[what,]
  res_data <- lapply(pfmodel$A %>% rownames(),function(sample){
    comps <- lapply(pfmodel$A %>% colnames(),function(component){
      pfmodel$B[,component] %*% t(pfmodel$C[,component]) * pfmodel$A[sample,component]
    })
    names(comps) <- pfmodel$C %>% colnames()
    fit <- comps %>%
      Reduce('+', .)
    eem <- eem_list[[which(eem_list %>% eem_names == sample)]]
    samp <- eem$x[eem$em %in% rownames(pfmodel$B),eem$ex %in% rownames(pfmodel$C)]
    res <- samp - fit

    comps <- lapply(pfmodel$A %>% colnames(),function(component){
      comps[[component]] %>% data.frame() %>% mutate(type = component, em = rownames(pfmodel$B)) %>%
        gather(ex,value,-em,-type) %>% mutate(ex = substr(ex,2,4))
    }) %>%
      bind_rows()
    colnames(samp) <- rownames(pfmodel$C)
    samp <- samp %>% data.frame() %>% mutate(type = "sample", em = rownames(pfmodel$B)) %>%
      gather(ex,value,-em,-type) %>% mutate(ex = substr(ex,2,4))
    rownames(res) <- rownames(pfmodel$B)
    colnames(res) <- rownames(pfmodel$C)
    res <- res %>% data.frame() %>% mutate(type = "residual", em = rownames(pfmodel$B)) %>%
      gather(ex,value,-em,-type) %>% mutate(ex = substr(ex,2,4))
    bind_rows(list(comps,samp,res)) %>%
      mutate(Sample = sample)
  }) %>%
    bind_rows()
}

#' Calculate the amount of each component for samples not involved in model building
#'
#' @description Samples from an eemlist that were not used in the modelling process are added as entries in the  A-modes. Values are calculated using fixed B and C modes in the PARAFAC algorithm. B and C modes can be provided via a previously calculated model or as matrices manually.
#'
#' @param eem_list object of class eemlist with sample data
#' @param pfmodel object of class parafac
#' @param cores number of cores to use for parallel processing
#' @param components optionally supply components to use manually, either as a variable of class parafac_components or as a list of variables of class parafac_components, if you do so,
#' @param const optional constraints for model, just used, when components are supplied
#' @param control optional constraint control parameters for model, just used, when components are supplied
#' @param ... additional arguments passed to eem_parafac
#'
#' @details This function can be used to calculate A modes (sample loadings) for samples that were previously excluded from the modelling process (e.g. outliers). Another way to use it would be a recombination of components from different models and calculating the according sample loadings. Expecially the later application is experimental and results have to be seen critically! Nevertheless, I decided to supply this function to stimulate some experiments on that and would be interested in your findings and feedback.
#'
#' @return object of class parafac
#' @export
#'
#' @import dplyr
#' @import parallel
#' @import foreach
#' @importFrom doParallel registerDoParallel
#'
#' @examples
#' \donttest{
#' data(eem_list)
#' data(pf_models)
#'
#' A_missing(eem_list,pf4[[1]])
#' }
A_missing <- function(eem_list,pfmodel = NULL,cores = parallel::detectCores(logical = FALSE),components = NULL, const = NULL, control = NULL, ...){
  eem_list <- eem_red2smallest(eem_list)
  if(is.null(pfmodel) & is.null(components)) stop("You must either specify a model or components as a base for the newly generated model!")

  exclude <- list("ex" = eem_list[[1]]$ex[!(eem_list[[1]]$ex %in% rownames(pfmodel$C))],
                  "em" = eem_list[[1]]$em[!(eem_list[[1]]$em %in% rownames(pfmodel$B))],
                  "sample" = c()
  )
  x <- eem_list %>%
    eem_exclude(exclude)
  if(!is.null(components)){
    if(!is.null(pfmodel)) warning("The base model is ignored since you provided components manually!")
    if(class(components[[1]]) == "parafac_components") components <- eempf_bindxc(components)
    if(class(components) == "parafac_components"){
      Bfixed <- components$B
      Cfixed <- components$C
      comps <- ncol(components$B)
      if(is.null(const)) const <- c("nonneg", "nonneg", "nonneg")
    } else {
      stop("The list of components you supplied is invalid!")
    }
  } else {
    Bfixed <- pfmodel$B
    Cfixed <- pfmodel$C
    comps <- pfmodel$A %>% ncol()
    normalise = (!is.null(attr(pfmodel,"norm_factors")))
    control = pfmodel$control
    const = pfmodel$const
  }
  missingAs <- eem_parafac(x,comps = comps,normalise = (!is.null(attr(pfmodel,"norm_factors"))),Bfixed = Bfixed, Cfixed = Cfixed,cores = cores,const = const, control = control, ...)
  missingAs[[1]]
}

#' Running a Split-Half analysis on a PARAFAC model
#'
#' @description The samples are split into four subsamples: A,B,C,D. Subsamples are then combined and compared: AB vs. CD, AC vs. BD, AD vs. BC. The results show graphs from the components of each of the 6 models.
#'
#' @param eem_list eemlist containing sample data
#' @param comps number of desired components
#' @param splits optional, list of 4 numerical vectors containing the sample numbers for A,B,C and D sample subsets
#' @param rand logical, splits are randomised
#' @param normalise state whether EEM data should be normalised in advance
#' @param nstart number of random starts
#' @param cores number of parallel calculations (e.g. number of physical cores in CPU)
#' @param maxit maximum iterations for PARAFAC algorithm
#' @param ctol Convergence tolerance (R^2 change)
#' @param rescale rescale splithalf models to Fmax, see \code{\link[staRdom]{eempf_rescaleBC}}
#' @param verbose states whether you want additional information during calculation
#' @param ... additional parameters that are passed on to \code{\link[multiway]{parafac}}
#'
#' @details Split data sets can be split suboptimal and cause low TCCs. Therefore, subsamples are recombined in 3 different ways and a TCC close to 1 in only one split combination per component is already a positive result. Check the split sets to check for sample independency.
#'
#' @return data frame containing components of the splithalf models
#' @export
#'
#' @import dplyr
#' @import tidyr
#' @importFrom utils combn
#'
#' @seealso \code{\link[staRdom]{splithalf_plot}}, \code{\link[staRdom]{splithalf_tcc}}
#'
#' @examples
#' \donttest{
#' data(eem_list)
#'
#' splithalf <- splithalf(eem_list, comps = 6)
#' splithalf_plot(splithalf)
#' }
splithalf <- function(eem_list, comps, splits = NA, rand = FALSE, normalise = TRUE, nstart = 20, cores = parallel::detectCores(logical = FALSE), maxit = 2500, ctol = 10^(-7), rescale = TRUE, verbose = FALSE, ...){
  a <- seq(1,eem_list %>% length())
  if(rand){
    a <- a %>% sample()
  }
  if(is.na(splits[1])) splits <- lapply(seq(1:4),function(sp) a[seq(sp,length(a),by=4)] %>% sort())
  names(splits) <- LETTERS[1:length(splits)]

  spl_eems <- lapply(combn(seq(1:length(splits)),2) %>% split(rep(1:ncol(.), each = nrow(.))), function(co){
    eem_list %>% eem_extract(splits[co] %>% unlist() %>% sort(),keep=TRUE, verbose=FALSE)
  })

  split_designations <- c("AB","AC","AD","BC","BD","CD")

  names(spl_eems) <- split_designations

  if(verbose){
    cat(paste0(cores," cores are used for the calculation."),fill=TRUE)
    if(normalise){
      if(verbose) cat("EEM matrices are normalised!",fill=TRUE)
    }
    cat(paste0("Calculating PARAFAC models with split-half data..."),fill=TRUE)
    pb <- txtProgressBar(max = length(spl_eems), style=3)
  }
  fits <- lapply(1:length(spl_eems),function(i){
    # i <- 1
    mod <- eem_parafac(spl_eems[[i]], comps = comps, normalise = normalise, maxit = maxit, nstart = nstart, cores = cores,ctol = ctol, verbose = FALSE)#,...
    if(rescale) mod <- lapply(mod,eempf_rescaleBC,newscale="Fmax")
    if(verbose) setTxtProgressBar(pb, i)
    mod
  }) #
  if(verbose) close(pb)

  sscs <- eempf_ssc(fits, tcc = TRUE)

  sscs2 <- sscs %>%
    lapply(lapply,ssc_max)

  C_sort <- sscs2 %>%
    .[grepl("1vs",names(.))] %>%
    lapply(`[[`,2) %>%
    lapply(attr, "order")

  fits <- lapply(1:length(spl_eems),function(sel){
    lapply(fits[[sel]],function(f){
      f$A <- f$A[,C_sort[[sel]]]
      f$B <- f$B[,C_sort[[sel]]]
      f$C <- f$C[,C_sort[[sel]]]
      f
    })
  })

  sel_comb <- lapply(1:(length(fits)/2), function(i){
    paste0(split_designations[i],"vs",split_designations[length(fits) + 1 - i]) %>%
      setNames(paste0(i,"vs",length(fits) + 1 - i))
  }) %>%
    unlist()

  attr(fits,"tcc_table") <- sscs2 %>%
    .[names(sel_comb)] %>%
    lapply(lapply,data.frame) %>%
    lapply(bind_cols) %>%
    lapply(setNames, c("tcc_ex", "tcc_em")) %>%
    lapply(mutate, component = paste0("Comp.",1:n())) %>%
    bind_rows(.id = "comb") %>%
    mutate(comb = sel_comb[comb]) %>%
    select(component,comb,tcc_em,tcc_ex) %>%
    arrange(component,comb)
  attr(fits,"splits") <- lapply(spl_eems,eem_names) %>% setNames(c("AB","AC","AD","BC","BD","CD"))

  fits
}

#' Reorders components of different PARAFAC models according to best fit (TCC)
#'
#' @description When running a splithalf analysis similar components are not necessarily on the same position. This function looks for best fits with Tucker's Congruence Coefficients and returns a list of models with reordered components.
#'
#' @param fits list of parafac models
#'
#' @return list of parafac models
#' @export
#'
#' @seealso \code{\link[staRdom]{splithalf}}
#'
#' @importFrom stats setNames
#' @import dplyr
#' @import tidyr
#'
#' @examples
#' \donttest{
#' data(eem_list)
#'
#' # function currently only used from within splithalf
#' splithalf(eem_list,6,nstart=2)
#' }
tcc_find_pairs <- function(fits){
  warning("This function is deprecated! Please use eempf_ssc and ssc_max.")
  sel <- 0
  problem <- FALSE
  table <- lapply(fits,function(fit){
    sel <<- sel + 1
    c <- fit %>% lapply(eempf_comp_mat)
    tab <- lapply(c,function(c1){
      nc1 <- length(c1)
      nc2 <- 0
      lapply(c1,function(c2){
        nc2 <<- nc2 + 1
        c2 <- c2 %>%
          mutate(comps = nc1, comp = paste0("Comp.",nc2))
      }) %>%
        bind_rows()
    }) %>%
      bind_rows() %>%
      mutate(selection = sel, ex = as.numeric(ex), em = as.numeric(em)) %>%
      group_by(comp) %>%
      filter(em == em[which.max(value)] | ex == ex[which.max(value)]) %>%
      mutate(em2 = em, em = ifelse(ex != ex[which.max(value)],NA,em),ex = ifelse(em2 != em2[which.max(value)],NA,ex)) %>%
      select(-em2) %>%
      ungroup()
  }) %>%
    bind_rows()

  comps <- max(table$comps)

  ct <- lapply(c("ex","em"),function(e){
    t <- table %>%
      mutate(o = paste0(selection,comp)) %>%
      filter_(.dots=paste0("!is.na(",e,")")) %>%
      select(one_of(!!e),value,o) %>%
      arrange_(.dots=e) %>%
      spread(o,value)
    comb_nam <- colnames(t) %>% .[. != e]
    t <- t %>%
      select(one_of(comb_nam)) %>%
      multiway::congru()
    rownames(t) <- comb_nam
    colnames(t) <- comb_nam
    attr(t,"e") <- e
    t
  }) %>%
    setNames(c("ex","em"))

  arr <- ct %>%
    .[["em"]] %>%
    .[1:comps,(comps+1):(length(fits)*comps)] %>%
    data.frame()

  a <- data.frame("comp1" = arr %>% rownames(),stringsAsFactors = FALSE)
  for(i in 2:6){
    b <- c()
    for(j in 1:(arr %>% nrow())){
      sel <- (arr %>% colnames() %>% substr(2,2) == i) & (!arr %>% colnames() %in% b)
      if(sum(sel) == 1){
        d <- arr %>% colnames() %>% .[sel]
      }else{
        d <- arr[j,which(arr %>% colnames() %>% substr(2,2) == i & !arr %>% colnames() %in% b)] %>% which.max() %>% names()
      }
      b <- c(b,d)
    }
    a <- a %>% bind_cols(data.frame(i=b,stringsAsFactors = FALSE))
  }

  arrange_emex <- a %>%
    mutate_all(as.character()) %>%
    gather(i,comp2,-comp1) %>%
    mutate(selection = substr(comp2,2,2), comp2 = substr(comp2,3,8))%>%
    select(comp1,selection,comp2) %>% #View()
    rowwise() %>%
    mutate(tcc_ex = ct[[1]][comp1,paste0(selection,comp2)],tcc_em = ct[[2]][comp1,paste0(selection,comp2)]) %>%
    ungroup()

  if(length(fits) == 6){
    tt <- arrange_emex %>%
      select(comp1) %>%
      distinct(comp1) %>%
      mutate(selection = "1", comp2 = substr(comp1,2,7)) %>%
      bind_rows(arrange_emex) %>%
      select(-tcc_em,-tcc_ex) %>%
      group_by(comp1) %>%
      spread(selection,comp2) %>%
      ungroup()

    split_designations <- c("AB","AC","AD","BC","BD","CD")

    ttt <- lapply(1:3,function(i){
      tt %>%
        select(comp1,one_of(paste0(i),paste0(7-i))) %>%
        mutate(i=i, comb = paste0(split_designations[i],"vs",split_designations[7-i])) %>%
        setNames(c("set","x","y","i","comb"))
    }) %>%
      bind_rows() %>%
      mutate(component = substr(set,2,7)) %>%
      rowwise() %>%
      mutate(tcc_ex = ct[[1]][paste0(i,x),paste0(7-i,y)],tcc_em = ct[[2]][paste0(i,x),paste0(7-i,y)]) %>%
      arrange(component) %>%
      select(component,comb,tcc_ex,tcc_em)
    attr(arrange_emex,"tcc_table") <- ttt
  }
  arrange_emex
}


#' Extracting TCC values from a splithalf analysis
#'
#' @param fits list of parafac models (from a splithalf analysis)
#'
#' @return data frame containing TCC values
#' @export
#'
#' @examples
#' data(sh)
#'
#' splithalf_tcc(sh)
splithalf_tcc <- function(fits){
  attr(fits,"tcc_table")
}

#' Extracting a list of sample names in each subsample from a splithalf analysis
#'
#' @param fits list of parafac models (from a splithalf analysis)
#'
#' @return data frame containing TCC values
#' @export
#'
#' @examples
#' data(sh)
#' splithalf_splits(sh)
splithalf_splits <- function(fits){
  attr(fits,"splits")
}


#' Caluclate Tucker's Congruence Coefficient of PARAFAC components
#'
#' @description Componets must be passed as modes, see \code{\link[staRdom]{maxlines}}
#'
#' @param maxl_table data frame containing the peak lines of components
#' @param na.action if "na.omit" NA are deleted from prior the test
#'
#' @return data.frame containing the TCCs
#' @export
#'
#' @importFrom stats setNames
#' @importFrom multiway congru
#' @import dplyr
#' @import tidyr
#' @importFrom stats na.omit
#'
#' @examples
#' data(pf_models)
#'
#' ml <- maxlines(pf4[[1]])
#'
#' tcc(ml)
tcc <- function(maxl_table,na.action="na.omit"){
  c <- lapply(c("em","ex"), function(E) {
    c2 <- maxl_table %>% filter(e == E) %>% select(-e) %>% arrange(wavelength)
    if(na.action == "na.omit") c2 <- c2 %>% na.omit()
    c2 <- c2 %>%
      select(-wavelength) %>%
      congru()
    na <- maxl_table %>% select(-e,-wavelength) %>% names()
    rownames(c2) <- na
    colnames(c2) <- na
    c2
  }) %>%
    setNames(c("em","ex"))
  c
}


#' Write out PARAFAC components to submit to openfluor.org.
#'
#' @description openfluor.org offers the possibility to compare your results to others, that were uploaded to the database. This functions writes out a txt containing the header lines and your components. Please open the file in an editor and fill in further information that cannot be covered by this function.
#'
#' @param pfmodel PARAFAC model
#' @param file string, path to outputfile. The directory must exist, the file will be created or overwritten if already present.
#' @param Fmax rescale modes so the A mode shows the maximum fluorescence. As openfluor does not accept values above 1, this is a way of scaling the B and C modes to a range between 0 and 1.
#'
#' @return txt file
#' @export
#'
#' @importFrom stringr str_replace_all
#'
#' @examples
#'   data(pf_models)
#'   eempf_openfluor(pf4[[1]],file.path(tempdir(),"openfluor_example.txt"))
eempf_openfluor <- function(pfmodel, file, Fmax = TRUE){
  if(!dir.exists(dirname(file.path(file)))){
    stop("The path to your file does not contain an existing directory. Please enter a correct path!")
  }
  factors <- rbind(pfmodel$C %>%
                     data.frame(mode = "Ex", wl = rownames(.),.),
                   pfmodel$B %>%
                     data.frame(mode = "Em", wl = rownames(.),.))
  template <- system.file("openfluor_template.txt",package="staRdom")
  template <- readLines(template)
  template <- stringr::str_replace_all(template,"toolbox","toolbox\tstaRdom")
  template <- stringr::str_replace_all(template,"nSample",paste0("nSample\t",pfmodel$A %>% nrow()))
  write(template,file)
  write.table(factors,file=file,append=TRUE,col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE)
  message("An openfluor file has been successfully written. Please fill in missing header fields manually!")
}


#' Create table of PARAFAC components and (optionally) EEM peaks and indices as well as absorbance slope parameters.
#'
#' @description Please refer to \code{\link[eemR]{eem_biological_index}}, \code{\link[eemR]{eem_coble_peaks}}, \code{\link[eemR]{eem_fluorescence_index}}, \code{\link[eemR]{eem_biological_index}} and \code{\link[staRdom]{abs_parms}} for details on the certain values
#'
#' @param pfmodel PARAFAC model where loadings of the components are extracted
#' @param eem_list optional eemlist used for peak and indices calculation
#' @param absorbance optional absorbance table used for absorbance slope parameter calculation
#' @param cuvl optional cuvette length of absorbance data in cm
#' @param n optional size of moving window in nm for data smoothing in advance of peak picking
#' @param export optional file path of csv or txt table where data is exported
#' @param ... additional parameters passed to \code{\link[utils]{write.table}}
#'
#' @return data frame
#' @export
#'
#' @import dplyr
#' @import tidyr
#' @import eemR
#' @importFrom tibble rownames_to_column
#'
#' @examples
#' \donttest{
#' data(eem_list)
#' data(pf_models)
#'
#' results <- eempf4analysis(pfmodel = pf4[[1]],
#'                           eem_list = eem_list,
#'                           cuvl = 5, n = 4)
#'                           }
eempf4analysis <- function(pfmodel,eem_list = NULL, absorbance = NULL, cuvl = NULL, n = 4, export = NULL,...){
  loadings <- pfmodel$A %>%
    norm2A() %>%
    data.frame() %>%
    tibble::rownames_to_column("sample")
  if(!is.null(eem_list)){
    eem4peaks <- eem_list %>% eem_smooth(n=n)
    indices_peaks <- eem4peaks %>% eem_biological_index() %>%
      full_join(eem4peaks %>% eem_coble_peaks(), by="sample")  %>%
      full_join(eem4peaks %>% eem_fluorescence_index(), by="sample") %>%
      full_join(eem4peaks %>% eem_humification_index(scale=TRUE), by="sample")
    loadings <- full_join(loadings,indices_peaks,by="sample")
  }
  if(!is.null(absorbance)){
    if(is.null(cuvl)){
      warning("Because of missing cuvette length, absorbance slope parameters were not calculated!")
    } else {
      abs_parameters <- abs_parms(absorbance %>%
                                    select(one_of(names(absorbance))), cuvl)
      loadings <- full_join(loadings,abs_parameters,by="sample")
    }
  }
  if(!is.null(export)){
    write.table(loadings, file=export, row.names=FALSE, ...)
  }
  loadings
}


#' Create one table containing the PARAFAC models factors and optionally exporting it to csv or txt
#'
#' @param pfmodel PARAFAC model
#' @param export file path to export table
#' @param Fmax rescale modes so the A mode shows the maximum fluorescence
#' @param ... additional parameters passed to \code{\link[utils]{write.table}}
#'
#' @return data frame
#' @export
#'
#' @importFrom tibble rownames_to_column
#' @import dplyr
#' @import tidyr
#' @importFrom utils write.table
#'
#' @examples
#' data(pf_models)
#'
#' factor_table <- eempf_export(pf4[[1]])
eempf_export<- function(pfmodel,export = NULL, Fmax = TRUE,...){
  pfmodel <- norm2A(pfmodel)
  if(Fmax) pfmodel <- eempf_rescaleBC(pfmodel)
  tabs <- list(pfmodel$A %>%
                 as.data.frame() %>%
                 tibble::rownames_to_column("sample") #%>%
               #`colnames<-`(colnames(.) %>% stringr::str_replace_all("Comp.","Fmax"))
               ,
               pfmodel$B %>%
                 as.data.frame() %>%
                 tibble::rownames_to_column("Em"),
               pfmodel$C %>%
                 as.data.frame() %>%
                 tibble::rownames_to_column("Ex"))
  rows <- lapply(tabs,nrow) %>% unlist() %>% max()
  tabs <- lapply(tabs,function(tab){
    if(nrow(tab) < rows){
      tab <- rbind(tab,as.data.frame(matrix(nrow = rows - nrow(tab), ncol = ncol(tab))) %>% `colnames<-`(colnames(tab)))
    }
    tab <- cbind(tab,as.data.frame(matrix(nrow = rows, ncol = 1)) %>% `colnames<-`(" "))
  })
  tabs <- do.call(cbind,tabs)
  if(!is.null(export)) write.table(tabs,file=export,row.names = FALSE,na="",...)
  tabs %>% invisible()
}


#' Calculate the core consistancy of an EEM PARAFAC model
#'
#' @description This is basically a wrapper for \code{\link[multiway]{corcondia}} that deals with the normalisation of the original data., Other than \code{\link[multiway]{corcondia}}, the default dicisor = "core".
#'
#'
#' @param pfmodel PARAFAC model
#' @param eem_list eemlist
#' @param divisor divisor, please refer to \code{\link[multiway]{corcondia}}
#'
#' @return numeric
#' @export
#'
#' @importFrom multiway corcondia
#' @import dplyr
#' @import tidyr
#'
#' @examples
#' \dontrun{
#' # due to data limitation in package, example does not work with that data!
#'
#' # eempf_corcondia(pfmodel,eem_list)
#'
#' }
eempf_corcondia <- function(pfmodel,eem_list,divisor="core"){
  arr <- eem_list %>% eem2array()
  if(!is.null(attr(eem_list,"norm_factors"))) arr <- arr %>% norm_array()
  corcondia(arr,pfmodel,divisor=divisor)
}

#' Calculating EEMqual which is an indicator of a PARAFAC model's quality
#'
#' @param pfmodel PARAFAC model
#' @param eem_list EEM data as eemlist
#' @param splithalf optionally, you can supplie available splithalf results from model to decrease computation time
#' @param ... additional arguments passed to splithalf
#'
#' @return data frame containing fit, corcondia, product of best TCCs from splithalf analysis, eemqual and splithalf models
#'
#' @export
#'
#' @importFrom multiway corcondia
#' @import dplyr
#'
#' @references Rasmus Bro, Maider Vidal, EEMizer: Automated modeling of fluorescence EEM data, Chemometrics and Intelligent Laboratory Systems, Volume 106, Issue 1, 2011, Pages 86-92, ISSN 0169-7439
#'
#' @examples
#' \donttest{
#' data(eem_list)
#' data(pf_models)
#'
#' pfmodel <- pf4[[1]]
#' eempf_eemqual(eem_list,pfmodel) # insuficient example data to run!
#' }
eempf_eemqual <- function(pfmodel,eem_list,splithalf = NULL, ...){
  corec <- eempf_corcondia(pfmodel,eem_list)
  fit <- pfmodel$Rsq
  if(is.null(splithalf)) sh <- splithalf(eem_list,comps = pfmodel$A %>% ncol(),...) else sh <- splithalf
  splth <- splithalf_tcc(sh) %>%
    group_by(component) %>%
    summarise_at(vars(contains("tcc")),max,na.rm=TRUE) %>%
    .[c(2,3)] %>%
    prod()
  tab <- data.frame(components = ncol(pfmodel$A) ,fit = fit, corec = corec, splithalf = splth, eemqual = fit*corec*splth)
  attr(tab,"shmodel") <- sh
  tab
}

#' Calculate the importance of each component.
#'
#' @param pfmodel model of class parafac
#' @param eem_list eemlist used to calculate that model
#' @param cores cores to be used for the calculation
#' @param ... other aruments passed to eem_parafac
#'
#' @details The importance of each variable is calculated by means of creating a model without a specific component and calculating the difference between the original R-squared and the one with the left out component. The derived values state the loss in model fit if one component is not used in the modeling process. For the creation of the new models, the exact components of the original model are used.
#'
#' @return numeric vector, values are in the same order of the components in the supplied model.
#' @export
#'
#' @import dplyr
#' @import tidyr
#'
#' @examples
#' \donttest{
#' data(pfmodel)
#' data(eem_list)
#'
#' eempf_varimp(pf4[[1]],eem_list)
#' }
eempf_varimp <- function(pfmodel, eem_list, cores = parallel::detectCores(logical=FALSE),...){

  exclude <- list("ex" = eem_list[[1]]$ex[!(eem_list[[1]]$ex %in% rownames(pfmodel$C))],
                  "em" = eem_list[[1]]$em[!(eem_list[[1]]$em %in% rownames(pfmodel$B))],
                  "sample" = c()
  )
  x <- eem_list %>%
    eem_red2smallest() %>%
    eem_exclude(exclude)
  mods <- lapply(1:(pfmodel$A %>% ncol()), function(c){
    eem_parafac(x, comps = ncol(pfmodel$A)-1,normalise = (!is.null(attr(pfmodel,"norm_factors"))),Bfixed = pfmodel$B[,-c], Cfixed = pfmodel$C[,-c],cores = cores,control = pfmodel$control, const = pfmodel$const, ...)
  })
  Rsq_al <- pfmodel$Rsq - (mods %>% lapply(function(mod) {mod[[1]]$Rsq}) %>% unlist())
}


#' Reorder PARAFAC components
#'
#' @param pfmodel model of class parafac
#' @param order vector containing desired new order or "em" or "ex" to reorder according to emission or excitation wavelengths of the peaks
#' @param decreasing logical, whether components are reordered according to peak wvalengths in a decreasing direction
#'
#' @return parafac model
#'
#' @import dplyr
#' @importFrom multiway reorder.parafac
#'
#' @export
#'
#' @examples
#' data(pf_models)
#' ggeem(pf4[[1]])
#'
#' pf4r <- eempf_reorder(pf4[[1]],"ex")
#' ggeem(pf4r)
eempf_reorder <- function(pfmodel,order,decreasing = FALSE){
  if(!(order[1] == "em" | order[1] == "ex" | is.vector(order))) stop("no valid data suppli ed for order!")
  if(order[1] == "em") order <- apply(pfmodel$B,2,which.max) %>% sort.list(decreasing = decreasing)
  if(order[1] == "ex") order <- apply(pfmodel$C,2,which.max) %>% sort.list(decreasing = decreasing)
  if(ncol(pfmodel$A) != length(order)) stop("the length of the order vector does not fit the number of components")

  mod <- try(reorder.parafac(pfmodel,neworder = order), silent=TRUE)
  if(class(mod) =="try-error") stop(mod) else {
    mod
  }
}


#' Extracting components of a PARAFAC model
#'
#' @param pfmodel parafac model
#' @param comps vector with numbers of components to extract
#'
#' @return list
#' @export
#'
#' @import dplyr
#'
#' @examples
#' data(pf_models)
#' pfmodel <- pf4[[1]]
#' comps <- eempf_excomp(pfmodel,c(1,3))
eempf_excomp <- function(pfmodel,comps){
  list(pfmodel$B[,comps], pfmodel$C[,comps]) %>%
    `names<-`(c("B","C")) %>%
    `class<-`("parafac_components")
}


#' Combining extracted components of PARAFAC models
#'
#' @param components list of parafac_components
#'
#' @return parafac_components
#' @export
#'
#' @import dplyr
#'
#' @examples
#'
#' data(pf_models)
#' pfmodel <- pf4[[1]]
#' comps <- eempf_excomp(pfmodel,c(1,3))
#' comps2 <- eempf_excomp(pfmodel,c(4,6))
#' comps3 <- eempf_bindxc(list(comps, comps2))
#'
eempf_bindxc <- function(components){
  B <- lapply(components, `[[`, "B") %>%
    do.call(cbind,.)
  C <- lapply(components, `[[`, "C") %>%
    do.call(cbind,.)
  list(B,C) %>%
    `names<-`(c("B","C")) %>%
    `class<-`("parafac_components")
}

#' Calculate the shift-and shape-sensitive congruence (SSC) between model components
#'
#' @param pfmodels list of either PARAFAC models or component matrices
#' @param tcc if set TRUE, TCC is returned instead
#' @param m logical, if TRUE, emission and excitation SSCs or TCCs are combined by calculating the geometric mean
#' @param cores number of CPU cores to be used
#'
#' @description The data variable pf_models can be supplied as list of PARAFAC models, output from a splithalf analysis or list of matrices
#' Please see details of calculation in:
#' U.J. Wünsch, R. Bro, C.A. Stedmon, P. Wenig, K.R. Murphy, Emerging patterns in the global distribution of dissolved matter fluorescence, Anal. Methods, 11 (2019), pp. 888-893
#'
#' @return (list of) tables containing SCCs between components
#' @export
#'
#' @import parallel
#' @import dplyr
#'
#' @examples
#' \donttest{
#' pf_models <- pf3[1:3]
#'
#' sscs <- eempf_ssc(pf_models, cores = 2)
#' sscs
#'
#' tcc <- eempf_ssc(pf_models, tcc = TRUE, cores = 2)
#' tcc
#' ## mixed tcc (combine em and ex)
#' mtcc <- eempf_ssc(pf_models, tcc = TRUE, m = TRUE, cores = 2)
#' mtcc
#'
#' ## compare results from splithalf analysis
#' sh_sscs <- eempf_ssc(sh, cores = 2)
#'
#' sh_sscs
#' ## view diagonals only (components with similar numbers only)
#' lapply(sh_sscs, lapply, diag)
#' }
eempf_ssc <- function(pfmodels, tcc = FALSE, m = FALSE, cores = parallel::detectCores(logical = FALSE)){
  classes <- unlist(lapply(unlist(pfmodels, recursive = FALSE),class))
  if(any(classes == "parafac") & !is.null(classes)){ ## Results from splithalf
    pfmodels %>%
      unlist(recursive = FALSE) %>%
      lapply(function(mod){
        list(B=mod$B,C=mod$C)
      }) %>%
      eempf_ssc(tcc = tcc, m = m, cores = cores)
  } else if(all(unlist(lapply(pfmodels,class)) == "parafac") & !is.null(unlist(lapply(pfmodels,class)))){ ## PARAFAC models
    pfmodels %>% lapply(function(mod){
      list(B=mod$B,C=mod$C)
    }) %>%
      eempf_ssc(pfmodels = ., tcc = tcc, m = m, cores = cores)
  } else if(all(classes == "matrix")){ ## matrices

    cl <- makePSOCKcluster(min(cores,length(pfmodels)))
    clusterExport(cl, c("pfmodels","tcc"), envir = environment())
    clusterEvalQ(cl,require(staRdom))

    SSCs <- parLapply(cl,1:length(pfmodels),function(k){
      lapply(k:length(pfmodels), function(l){
        B = ssc(pfmodels[[k]][[1]], pfmodels[[l]][[1]], tcc = tcc)
        C = ssc(pfmodels[[k]][[2]], pfmodels[[l]][[2]], tcc = tcc)
        list(B=B, C=C)
      }) %>%
        setNames(paste0(k,"vs",k:length(pfmodels)))
    }) %>% unlist(recursive = FALSE)
    stopCluster(cl)

    if(m){
      SSCs <- lapply(SSCs, function(mats){
        sqrt(mats[[1]]*mats[[2]])
      })
    }
    SSCs
  } else {
    stop("No suitable data supplied! Please refer to the eempf_ssc help.")
  }
}

#' Calculate the shift-and shape-sensitive congruence (SSC) between two matrices
#'
#' @param mat1 matrix
#' @param mat2 matrix
#' @param tcc if set TRUE, TCC is returned instead
#'
#' @description Please see details in:
#' U.J. Wünsch, R. Bro, C.A. Stedmon, P. Wenig, K.R. Murphy, Emerging patterns in the global distribution of dissolved matter fluorescence, Anal. Methods, 11 (2019), pp. 888-893
#'
#' @return table containing pairwise SCC of matrices columns
#' @export
#'
#' @examples
#' pf_models <- pf3
#' mat1 <- pf_models[[1]][[2]]
#' mat2 <- pf_models[[2]][[2]]
#'
#' ## calculate SSC
#' ssc(mat1,mat2)
#'
#' ## calculate TCC
#' ssc(mat1,mat2, tcc = TRUE)
#'
ssc <- function(mat1, mat2, tcc = FALSE){
  if(any(is.null(mat1),is.na(mat1),is.null(mat2), is.na(mat2))){
    a <- NA
  } else {
    a <- lapply(1:ncol(mat1),function(nc){
      col1 <- mat1[,nc]
      apply(mat2,2,function(col2){
        tcc_cal <- sum(col1*col2)/sqrt(sum(col1^2)*sum(col2^2))
        if(!tcc){
          wl <- as.numeric(names(col1))
          if(any(is.na(wl)) | pracma::isempty(wl)){
            stop("SSCs cannot be calculated. Please add wavelengths as rownames of the matrices!")
          }
          alpha <- abs((wl[which.max(col1)]-wl[which.max(col2)]) / diff(range(wl)))
          beta <- abs((sum(col1/max(col1)) - sum(col2/max(col2))) / diff(range(wl)))
          ssc <- tcc_cal -alpha - beta
        } else {
          tcc_cal
        }
      })
    }) %>% setNames(colnames(mat1)) %>%
      do.call(rbind,.)
  }
  attr(a,"method") <- ifelse(tcc, "TCC", "SSC")
  a
}


#' Check SSCs between different models or initialisations of one model
#'
#' @param pfmodels list of parafac models
#' @param best number of models with the highest R^2 to be used, default is all models
#' @param tcc logical, if TRUE, TCC instead of SSC is calculated
#' @param cores number of CPU cores to be used
#'
#' @return data.frame containing SSCs
#' @export
#'
#' @import dplyr
#' @importFrom stringr str_extract
#'
#' @examples
#' \donttest{
#' data(pf_models)
#'
#' eempf_ssccheck(pf3[1:2], cores = 2)
#'
#' # SSCs of split-half models, models need to be unlisted
#' data(sh)
#' eempf_ssccheck(unlist(sh, recursive = FALSE), cores = 2)
#' }
eempf_ssccheck <- function(pfmodels, best = length(pfmodels), tcc = FALSE, cores = parallel::detectCores(logical = FALSE)){
  #pfmodels <- pf3
  Rsqs <- lapply(pfmodels,`[[`,"Rsq") %>% unlist()
  checkmods <- pfmodels[order(Rsqs, decreasing = TRUE)[1:best]]
  not_conv <- checkmods %>%
    lapply(`[[`,"cflag") %>%
    unlist() %>%
    sapply(identical, 0) %>%
    sapply(`!`) %>%
    sum()
  if(not_conv){
    warning(paste0(not_conv," of the best ", best, " chosen models ",ifelse(not_conv == 1, "is","are")," not converging!"))
  }
  #Rsqs <- lapply(checkmods,`[[`,"Rsq") %>% unlist()
  sscs <- eempf_ssc(pfmodels = checkmods, tcc = tcc, cores = cores)

  cl <- makePSOCKcluster(min(cores, length(sscs)))
  clusterExport(cl, c("sscs"), envir=environment())
  clusterEvalQ(cl,require(staRdom))

  maxs <- parLapply(cl, sscs, lapply, ssc_max)

  stopCluster(cl)

  a <- names(maxs) %>%
    lapply(function(na){
      str_extract(na,"^[0-9]{1,2}") != str_extract(na,"[0-9]{1,2}$")
    }) %>%
    unlist() %>%
    maxs[.] %>%
    lapply(bind_rows)
  ssccheck <- a %>%
    lapply(mutate, comp = 1:n()) %>%
    bind_rows() %>%
    mutate(comparison = names(a)[cumsum(comp == 1)])
  attr(ssccheck,"method") <- ifelse(tcc, "TCC", "SSC")
  ssccheck
}

#' Calculate the combination of components giving the maximum of geometric mean of TCCs
#'
#' @param mat matrix
#'
#' @importFrom gtools permutations
#' @import dplyr
#'
#' @return vector with TCCs having the highest possible geometric mean
#' @export
#'
#' @examples
#' mat <- matrix(c(7,2,13,6,0,7,1,5,5), nrow = 3)
#' mat
#'
#' sscs <- ssc_max(mat)
#' sscs
#'
#' # order of components:
#' attr(sscs,"order")
ssc_max <- function(mat){
  n <- min(dim(mat))
  p <- permutations(n , n)
  combinations <- lapply(1:nrow(p),function(row){
    per <- p[row,]
    matrix(c(1:n,per), ncol = n,byrow = TRUE)
  })

  best_comb <- lapply(combinations, function(c){
    pair <- c[,2] %>% unlist()
    apply(c,2,function(pair){
      mat[pair[1],pair[2]] %>% `^`(2)
    }) %>%
      sum() %>%
      sqrt()
  }) %>%
    which.max() %>%
    combinations[[.]]

  res <- apply(best_comb,2,function(pair){
    mat[pair[1],pair[2]]
  })
  attr(res,"order") <- best_comb[2,]
  res
}

#' Extract modelling imformation from a PARAFAC model.
#'
#' @description The convergence behaviour of all initialisations in a PARAFAC model is shown by printing the numbers
#'
#' @param pfmodel PARAFAC model created with staRdom using output = "all"
#' @param print logical, whether you want console output or just a list with results
#'
#' @return list with numbers of converging models, cflags and SSEs
#' @export
#'
#' @examples
#' data("pf_models")
#'
#' pfmodel <- pf4[[1]]
#' conv_beh <- eempf_convergence(pfmodel)
eempf_convergence <- function(pfmodel, print = TRUE){
  if(!is.list(pfmodel$models)){
    stop("The supplied PARAFAC model does not contain the whole model set used in the calculation! Please rerun eem_parafac setting output = 'all'")
  } else {
    n <- length(pfmodel$models)
    sses <- lapply(pfmodel$models,`[[`,"SSE") %>% unlist()
    conv <- lapply(pfmodel$models,`[[`,"cflag") %>% unlist()
    res <- list(models = n, converging = sum(conv == 0), nc_itlim = sum(conv == 1), nc_other = sum(conv == 2), conv = conv, sses = sses)
    if(print){
    cat("Calculated models: ", n, fill = TRUE)
    cat("Converging models: ", sum(conv == 0), fill = TRUE)
    cat("Not converging Models, iteration limit reached: ", sum(conv == 1), fill = TRUE)
    cat("Not converging models, other reasons: ", sum(conv == 2), fill = TRUE)
    cat("Best SSE: ", min(sses), fill = TRUE)
    cat("Summary of SSEs of converging models:",fill = TRUE)
    summary(sses[conv == 0]) %>%
      print()
    }
    invisible(res)
  }
}

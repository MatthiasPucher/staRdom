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
#' @param const constraints of PARAFAC analysis
#' @param verbose print infos
#' @param ... additional parameters that are passed on to \code{\link[multiway]{parafac}}
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
#' nstart <- 10 # random starts for PARAFAC analysis, models built simulanuously, best selected
#' cores <- parallel::detectCores()/2 # use all cores but do not use all threads
#' maxit = 500
#' ctol <- 10^-4 # tolerance for parafac
#'
#' pfres_comps <- eem_parafac(eem_list,comps=seq(dim_min,dim_max),
#'     normalise = TRUE,maxit=10000,nstart=nstart,ctol=ctol,cores=cores)
#' }
eem_parafac <- function(eem_list,comps,maxit=500,normalise=TRUE,const=c(2,2,2),nstart = 10,ctol=10^-4,cores = parallel::detectCores()/2, verbose = FALSE, ...){
  #eem_list <- eem_list %>% eem_red2smallest()
  eem_array <- eem2array(eem_list)
  if(normalise){
    eem_array <- eem_array %>% norm_array()
    if(verbose) cat("EEM matrices were normalised!",fill=TRUE)
  }
  res <- lapply(comps,function(comp){
    if(verbose) cat(paste0(cores," cores are used for the calculation."),fill=TRUE)
    if(verbose) cat(paste0("calculating ",comp," components model..."),fill=TRUE)
    #comp <- 5
    #eem_array %>% dim()
    cl <- NULL
    if(cores > 1){
    cl <- makeCluster(cores, type="PSOCK")
    clusterExport(cl, c("eem_array","comp","maxit","nstart","const","ctol","cores"), envir=environment())
    clusterEvalQ(cl, library(multiway))
    }
    cpresult <- parafac(eem_array,nfac=comp,const = const,maxit = maxit,parallel = (cores > 1), cl=cl, ctol=ctol, nstart = nstart,...)
    if(cores > 1){
    stopCluster(cl)
    }
    attr(cpresult,"norm_factors") <- attr(eem_array,"norm_factors")
    rownames(cpresult$B) <- eem_list[[1]]$em
    rownames(cpresult$C) <- eem_list[[1]]$ex
    rownames(cpresult$A) <- eem_list %>% eem_names()
    labComp <- paste("Comp.",1:comp,sep="")
    colnames(cpresult$A) <- labComp
    colnames(cpresult$B) <- labComp
    colnames(cpresult$C) <- labComp
    return(cpresult)
  })
  mostattributes(res) <- attributes(eem_array)
  return(res)
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
  ## set NA to 0
  eem_array[is.na(eem_array) | is.nan(eem_array)] <- 0
  eem_array <- eem_array %>% aperm(perm = c(3,1,2), resize = TRUE)
  attr(eem_array,"em") <- eem_list[[1]]$em
  attr(eem_array,"ex") <- eem_list[[1]]$ex
  attr(eem_array,"samples") <- eem_list %>% sapply("[","sample") %>% unlist()
  attr(eem_array,"mdim") <- dim(eem_array)
  #attr(eem_array,"norm_factors") <- NULL
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
#'
#' @examples
#' data(eem_list)
#'
#' a <- eem2array(eem_list)
#' an <- norm_array(a)
norm_array <- function(eem_array){
  norm_factors <- lapply(1:(dim(eem_array)[1]),function(s) {
    apply(eem_array[s,,], c(1,2), function(x) x^2) %>% sum(na.rm=TRUE)
  }) %>% unlist() / (10^7)
  eem_array <- eem_array / norm_factors
  #eem_array <- lapply(1:(dim(eem_array)[1]),function(s) {
  #  eem_array[s,,]/norm_factors[s]
  #}) %>% unlist() %>%
  #  array(dim=dim(eem_array))
  attr(eem_array,"norm_factors") <- norm_factors
  eem_array
}


#' Extract EEM matrix for single components determined in the PARAFAC analysis
#'
#' @description The components of a PARAFAC analysis are extracted as a data frame
#'
#' @param cp_out object of class parafac
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
#' data(pfres_comps1)
#'
#' eempf_comp_mat(pfres_comps[[3]])
eempf_comp_mat <- function(cp_out,gather=TRUE){
  #b[[3]] -> cp_out
  #comp <- 1
  mat <- lapply(seq(1:ncol(cp_out$A)), function(comp){
    m <- matrix(cp_out$B[,comp]) %*% t(matrix(cp_out$C[,comp])) %>%
      data.frame()
    colnames(m) <- cp_out$C %>% rownames()
    rownames(m) <- cp_out$B %>% rownames()
    if(gather==TRUE) m <- m %>% rownames_to_column("em") %>% gather(ex,value,-em)
    m
  })
  names(mat) <- colnames(cp_out$B)
  mat
}

#' Calculate the leverage of each emission and excitation wavelength and each sample
#'
#' @param cp_out object of class parafac
#'
#' @return list of 3 named vectors (emission, excitation wavelengths and samples)
#' @export
#'
#' @importFrom pracma pinv
#' @importFrom stats setNames
#' @import dplyr
#'
#' @examples
#' data(pfres_comps1)
#'
#' eempf_leverage(pfres_comps[[2]])
eempf_leverage <- function(cp_out){
  cpl <- lapply(cp_out[c("A","B","C")], function(M) diag(M %*% pinv(t(M) %*% M) %*% t(M)))
  names <- list(rownames(cp_out$A),rownames(cp_out$B),rownames(cp_out$C))
  cpl <- lapply(1:(cpl %>% length()),function(i){ cpl[[i]] %>% setNames(names[[i]])}) %>%
    setNames(c("A","B","C"))
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
#' data(pfres_comps1)
#'
#' leverage <- eempf_leverage(pfres_comps[[2]])
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
    mutate(label = ifelse(leverage > q,x,NA))
}


#' Compensate for normalisation in C-modes
#'
#' @description Factors used for normalisation are saved separately in the PARAFAC models. With this function, the normalisation factors are combined with the A-modes of the model and removed as a separate vector. This means former normalisation is accounted for in the amount of each component in each sample. If no normalisation was done, the original model is returned without warning.
#'
#' @param cp_out object of class parafac
#'
#' @return object of class parafac
#' @export
#'
#' @examples
#' data(pfres_comps1)
#'
#' pfres_comps[[2]][[3]] <- norm2A(pfres_comps[[2]][[3]])
norm2A <- function(cp_out){
  if(!is.null(attr(cp_out,"norm"))){
    cp_out$A <- cp_out$A * attr(cp_out,"norm_factors")
    attr(cp_out,"norm_factors") <- NULL
  }
  cp_out
}


#' Calculating correlations between the component loadings in all samples (C-Modes).
#'
#' @param cp_out results from a PARAFAC analysis, class parafac
#' @param method method of correlation, passed to \code{\link[stats]{cor}}
#' @param ... passed on to \code{\link[stats]{cor}}
#'
#' @return matrix
#' @export
#'
#' @importFrom stats cor
#'
#' @examples
#' data(pfres_comps1)
#' eempf_cortable(pfres_comps[[2]])
eempf_cortable <- function(cp_out,method="pearson",...){
  cp_out %>%
    norm2A() %>%
    .$A %>%
    cor(method=method,...)
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
#' "sample" = c("sample3","sample5","sample14")
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
    #eem <- eem_list[[1]]
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

#' Extract data from emission and excitation wavelengths of the components of a PARAFAC model
#'
#' @description Data of wavelengths is returned. For each component the lines intersecting at the component maxima are returned.
#'
#' @param cp_out object of class parafac
#'
#' @return data frame
#' @export
#'
#' @import dplyr
#' @import tidyr
#'
#' @examples
#' data(pfres_comps1)
#'
#' ml <- maxlines(pfres_comps[[1]])
maxlines <- function(cp_out){
  #cp_out <- pfres_comps[[2]]
  maxl <- lapply(colnames(cp_out$C),function(comp){
    #comp <- colnames(cp_out$C)[1]
    em = (cp_out$B[,comp]*max(cp_out$C[,comp])) %>%
      data.frame(e="em", wavelength = as.numeric(names(.)),value=.)
    ex = (cp_out$C[,comp]*max(cp_out$B[,comp])) %>%
      data.frame(e="ex", wavelength = as.numeric(names(.)),value=.)
    res <- bind_rows(em,ex) %>%
      setNames(c("e","wavelength",comp))
  }) %>%
    list_join(by=c("e","wavelength"))
}


#' Calculate residuals of EEM data according to a certain model
#'
#' @param cp_out PARAFAC model of class parafac
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
#' data(pfres_comps2)
#'
#' eempf_residuals(pfres_comps2[[2]],eem_list)
#' }
eempf_residuals <- function(cp_out,eem_list,select=NULL, cores = parallel::detectCores(logical = FALSE)/2){
  # cp_out <- pfres_comps2[[which(comps==seq(dim_min,dim_max))]]
  cp_out <- norm2A(cp_out)
  if(!is.null(select)){
    eem_list <- eem_extract(eem_list,sample = select ,keep=TRUE,verbose = FALSE)
      #lapply(select,function(s) which(s == eem_list %>% eem_names())) %>%
      #unlist() %>%
      #unique() %>%
      #setdiff(seq(1,length(eem_list)),.)  %>% eem_extract(eem_list,.,verbose = FALSE)
  }
  if(!all(eem_names(eem_list) %in% rownames(cp_out$A)) | length(eem_list) == 0){
    cp_out <- A_missing(eem_list,cp_out,cores=cores)
  }
  what <- which(rownames(cp_out$A) %in% (eem_list %>% eem_names()))
  cp_out$A <- as.data.frame(cp_out$A)[what,]
  res_data <- lapply(cp_out$A %>% rownames(),function(sample){
    comps <- lapply(cp_out$A %>% colnames(),function(component){
      cp_out$B[,component] %*% t(cp_out$C[,component]) * cp_out$A[sample,component]
    })
    names(comps) <- cp_out$C %>% colnames()
    fit <- comps %>%
      Reduce('+', .)
    eem <- eem_list[[which(eem_list %>% eem_names == sample)]]
    samp <- eem$x[eem$em %in% rownames(cp_out$B),eem$ex %in% rownames(cp_out$C)]
    res <- samp - fit

    comps <- lapply(cp_out$A %>% colnames(),function(component){
      comps[[component]] %>% data.frame() %>% mutate(type = component, em = rownames(cp_out$B)) %>%
        gather(ex,value,-em,-type) %>% mutate(ex = substr(ex,2,4))
    }) %>%
      bind_rows()
    colnames(samp) <- rownames(cp_out$C)
    samp <- samp %>% data.frame() %>% mutate(type = "sample", em = rownames(cp_out$B)) %>%
      gather(ex,value,-em,-type) %>% mutate(ex = substr(ex,2,4))
    rownames(res) <- rownames(cp_out$B)
    colnames(res) <- rownames(cp_out$C)
    res <- res %>% data.frame() %>% mutate(type = "residual", em = rownames(cp_out$B)) %>%
      gather(ex,value,-em,-type) %>% mutate(ex = substr(ex,2,4))
    bind_rows(list(comps,samp,res)) %>%
      mutate(Sample = sample)
  }) %>%
    bind_rows()
}

#' Calculate the amount of each component for samples not involved in model building
#'
#' Samples from an eemlist that were not used in the modelling process are added as entries in the  A-modes. Values are calculated using fixed B and C modes in the PARAFAC algorithm.
#'
#' @param eem_list object of class eemlist with sample data
#' @param cp_out object of class parafac
#' @param cores number of cores to use for parallel processing
#' @param ... additional arguments passed to eem_parafac
#'
#' @return object of class parafac
#' @export
#'
#' @import parallel
#' @import foreach
#' @importFrom doParallel registerDoParallel
#'
#' @examples
#' \donttest{
#' data(eem_list)
#' data(pfres_comps2)
#'
#' A_missing(eem_list,pfres_comps2[[3]])
#' }
A_missing <- function(eem_list,cp_out,cores = parallel::detectCores(logical = FALSE)/2,...){
  #data("pfres_comps2")
  #data("eem_list")
  #cores=2
  #cp_out <- pfres_comps2[[2]]

  eem_list <- eem_red2smallest(eem_list)

  exclude <- list("ex" = eem_list[[1]]$ex[!(eem_list[[1]]$ex %in% rownames(cp_out$C))],
                  "em" = eem_list[[1]]$em[!(eem_list[[1]]$em %in% rownames(cp_out$B))],
                  "sample" = c()
  )

  x <- eem_list %>%
    eem_exclude(exclude)

  missingAs <- eem_parafac(x,comps = cp_out$A %>% ncol(),normalise = (!is.null(attr(cp_out,"norm_factors"))),Bfixed = cp_out$B, Cfixed = cp_out$C,cores = cores,control = cp_out$control, const = c(cp_out$const[1],0,0),...)

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
#' @param verbose states whether you want additional information during calculation
#' @param ... additional parameters that are passed on to \code{\link[multiway]{parafac}}
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
#' splithalf(eem_list,6,nstart=2)
#' }
splithalf <- function(eem_list,comps,splits=NA,rand=FALSE,normalise=TRUE,nstart=10,cores=parallel::detectCores()/2,maxit=500,ctol = 10^(-5),...,verbose =FALSE){
  a <- seq(1,eem_list %>% length())
  if(rand){
    a <- a %>% sample()
  }
  if(is.na(splits)) splits <- lapply(seq(1:4),function(sp) a[seq(sp,length(a),by=4)] %>% sort())
  spl_eems <- lapply(combn(seq(1:4),2) %>% split(rep(1:ncol(.), each = nrow(.))), function(co){

    eem_list %>% eem_extract(splits[co] %>% unlist() %>% sort(),keep=TRUE, verbose=FALSE)

  })
  if(verbose) cat(paste0(cores," cores are used for the calculation."),fill=TRUE)
  fits <- lapply(spl_eems,function(eem){eem_parafac(eem,comps=comps,normalise = normalise,maxit=maxit,nstart = nstart, cores = cores,ctol = ctol,verbose=verbose,...)}) #

  reallign <- tcc_find_pairs(fits)

  sel <- 4

  C_sort <- lapply(2:6,function(sel){
    fit <- fits[[sel]]
    reallign %>%
      filter(selection==sel) %>%
      arrange(comp1) %>%
      .$comp2 %>%
      substr(6,6) %>%
      as.numeric()
  })

  fits <- lapply(2:6,function(sel){
    lapply(fits[[sel]],function(f){
      f$A <- f$A[,C_sort[[sel - 1]]]
      f$B <- f$B[,C_sort[[sel - 1]]]
      f$C <- f$C[,C_sort[[sel - 1]]]
      f
    })
  })
  attr(fits,"tcc_table") <- attr(reallign,"tcc_table")
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
      congru()
    rownames(t) <- comb_nam
    colnames(t) <- comb_nam
    attr(t,"e") <- e
    t
  }) %>%
    setNames(c("ex","em"))

  arr <- ct %>%
    .[["em"]] %>%
    .[1:comps,(comps+1):(6*comps)] %>%
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


  split_designations <- c("AB","AC","AD","BC","BD","CD")

  tt <- arrange_emex %>%
    select(comp1) %>%
    distinct(comp1) %>%
    mutate(selection = "1", comp2 = substr(comp1,2,7)) %>%
    bind_rows(arrange_emex) %>%
    select(-tcc_em,-tcc_ex) %>%
    group_by(comp1) %>%
    spread(selection,comp2) %>%
    ungroup()

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
#' @description Componets must be passed as peak lines \code{\link[staRdom]{maxlines}}
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
#' data(pfres_comps1)
#'
#' ml <- maxlines(pfres_comps[[2]])
#'
#' tcc(ml)
tcc <- function(maxl_table,na.action="na.omit"){
  #maxl_table <- maxl %>% full_join(imp_table,by=c("e","wavelength"))
  #E="ex"
  c <- lapply(c("em","ex"), function(E) {
    c2 <- maxl_table %>% filter(e == E) %>% select(-e) %>% arrange(wavelength)
    if(na.action == "na.omit") c2 <- c2 %>% na.omit()
    #if(na.action == "approx") c2 <- c2 %>% zoo::na.approx(xout = .$wavelength, na.rm = FALSE) %>% data.frame() %>% na.omit()
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
#' @param cp_out PARAFAC model
#' @param file string, path to outputfile. The directory must exist, the file will be created or overwritten if already present.
#'
#' @return txt file
#' @export
#'
#' @importFrom stringr str_replace_all
#'
#' @examples
#'   data(pfres_comps2)
#'   eempf_openfluor(pfres_comps2[[2]],file.path(tempdir(),"openfluor_example.txt"))
eempf_openfluor <- function(cp_out,file){
  if(!dir.exists(dirname(file.path(file)))){
    stop("The path do your file does not contain an existing directory. Please enter a correct path!")
  }
  ml <- maxlines(cp_out) %>%
    arrange(desc(e)) %>%
    mutate(e = stringr::str_replace_all(e,"e","E"))
  template <- system.file("openfluor_template.txt",package="staRdom")
  template <- readLines(template)
  template <- stringr::str_replace_all(template,"toolbox","toolbox\tstaRdom")
  template <- stringr::str_replace_all(template,"nSample",paste0("nSample\t",cp_out$C %>% nrow()))
  write(template,file)
  write.table(ml,file=file,append=TRUE,col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE)
  message("An openfluor file has been successfully written. Please fill in missing header fields manually!")
}

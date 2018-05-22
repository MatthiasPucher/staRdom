#' Check for NAs in EEM data
#'
#' @param eem_list eemlist to check
#'
#' @return named character vector with sample names where EEM data contains NAs
#' @export
#'
#' @import dplyr
#' @import tidyr
#' @importFrom eemR eem_names
#' @importFrom stats setNames
#'
#' @examples
#' ### check
eem_is.na <- function(eem_list){
  lapply(eem_list,function(eem){
    any(is.na(eem$x))
  }) %>% unlist() %>%
    setNames(eem_names(eem_list))
}


#' Check for duplicate sample names
#'
#' @param data eemlist or data.frame containing absorbance data
#'
#' @return named character vector with duplicate sample names
#' @export
#'
#' @import dplyr
#' @import tidyr
#' @importFrom eemR eem_names
#'
#' @examples
#' ### check
eem_duplicates <- function(data) UseMethod("eem_duplicates")

#' @rdname eem_duplicates
#' @export
eem_duplicates.default <- function(data){
  stop("data is not of a suitable format!")
}

#' @rdname eem_duplicates
#' @export
#' @import dplyr
#' @import tidyr
#' @importFrom eemR eem_names
eem_duplicates.eemlist <- function(data){
  dupls <- eem_names(data) %>%
    table() %>%
    data.frame() %>%
    filter(Freq > 1) %>%
    .[,1] %>%
    as.character()
}

#' @rdname eem_duplicates
#' @export
#' @import dplyr
#' @import tidyr
#' @importFrom eemR eem_names
eem_duplicates.data.frame <- function(data){
  dupls <- colnames(data) %>%
    table() %>%
    data.frame() %>%
    filter(Freq > 1) %>%
    .[,1] %>%
    as.character()

  dots <- colnames(data) %>%
    grep("\\.[0-9]*$",.,value=TRUE) %>%
    gsub("(.)(\\.[0-9]*$){1}","\\1",.) %>%
    lapply(function(str){
      as <- colnames(data) %>%
        grep(paste0("^",str),.,value=TRUE)
      if(length(as > 1)) as else NA
    }) %>%
    unlist()

  dupls <- dupls %>%
    c(.,dots) %>% .[!is.na(.)]
}


#' Check your EEM, absorption and metadata before processing
#'
#' @description The function tries to lead you to possible problems in your data.
#'
#' @param eem_list eemlist continaing EEM data.
#' @param absorbance data.frame containing absorbance data.
#' @param metadata optional data.frame containing metadata.
#' @param metacolumns character vector of columns that are checkt for complete data sets
#' @param error logical, whether a problem should cause an error or not.
#'
#' @return logical, whether severe errors occured
#' @export
#'
#' @import dplyr
#' @import tidyr
#' @importFrom eemR eem_names eem_extract
#'
#' @examples
#' data(eem_list)
#' data(abs_data)
#'
#' eem_checkdata(eem_list,abs_data)
eem_checkdata <- function(eem_list,absorbance,metadata = NULL, metacolumns = NULL, error = TRUE){
  problem = FALSE
  nas <- which(eem_is.na(eem_list))
  if(length(nas) > 0){
    cat("NAs were found in the following samples: ",paste0(names(nas),collapse = ", ", sep=", "),fill=TRUE)
    problem = TRUE
  }

  eem_no_abs <- eem_names(eem_list)[!eem_names(eem_list) %in% colnames(absorbance)]
  if(length(eem_no_abs) > 0){
    cat(fill=TRUE)
    cat("EEM samples missing absorbance data:",fill=TRUE)
    lapply(eem_no_abs,function(dup){
      locs <- eem_list %>% eem_extract(dup, keep=TRUE, verbose=FALSE) %>%
        lapply(function(eem) eem$location)
      cat(dup,"in ",paste0(locs,collapse=", "), fill=TRUE)
    }) %>% invisible()
    problem = TRUE
  }

  abs_no_eem <- colnames(absorbance)[!colnames(absorbance) %in% eem_names(eem_list)] %>%
    .[. != "wavelength"]
  if(length(abs_no_eem) > 0){
    cat(fill=TRUE)
  cat("Absorbance data with missing EEM samples:",fill=TRUE)
  lapply(abs_no_eem,function(dup){
    locs <- attr(absorbance,"location")[which(colnames(absorbance) == dup) - 1]
    cat(dup,"in ",paste0(locs,collapse=", "), fill=TRUE)
  }) %>% invisible()
  cat("This can happen if you diluted for EEM and have additional undiluted absorbance samples.",fill=TRUE)
  }

  duplse <- eem_duplicates(eem_list)

  if(length(duplse) > 0) {
    cat(fill=TRUE)
    cat("Duplicate EEM sample names and according directories: ", fill=TRUE)
    lapply(duplse,function(dup){
      locs <- eem_list %>% eem_extract(dup, keep=TRUE, verbose=FALSE) %>%
        lapply(function(eem) eem$location)
      cat(dup,"in ",paste0(locs,collapse=", "), fill=TRUE)
    }) %>% invisible()
    problem = TRUE
  }

  duplsa <- eem_duplicates(absorbance)

  if(length(duplsa) > 0) {
    cat(fill=TRUE)
    cat("Duplicate absorbance sample names and according directories: ", fill=TRUE)
    lapply(duplsa,function(dup){
      locs <- attr(absorbance,"location")[which(colnames(absorbance) == dup) - 1]
      cat(dup,"in ",paste0(locs,collapse=", "), fill=TRUE)
    }) %>% invisible()
    cat("If sample names contain dots, please check if dots were in the original sample name or added by R due to duplicate sample names!")
    problem = TRUE
  }

  invalid_eem <- eem_names(eem_list)[!eem_names(eem_list) %in% make.names(eem_names(eem_list))]
  if(length(invalid_eem) > 0) {
    cat(fill=TRUE)
    cat("Invalid sample names in EEM data:",fill=TRUE)
    lapply(invalid_eem,function(inv){
      locs <- eem_list %>% eem_extract(inv, keep=TRUE, verbose=FALSE) %>%
        lapply(function(eem) eem$location)
      cat(inv,"in ",paste0(locs,collapse=", "), fill=TRUE)
    }) %>% invisible()
    problem = TRUE
  }

  invalid_abs <- colnames(absorbance)[!colnames(absorbance) %in% make.names(colnames(absorbance))]
  if(length(invalid_abs) > 0){
    problem = TRUE
    cat(fill=TRUE)
    cat("Invalid sample names in absorbance data:",fill=TRUE)
    lapply(invalid_abs,function(inv){
      locs <- attr(absorbance,"location")[which(colnames(absorbance) == inv) - 1]
      cat(inv,"in ",paste0(locs,collapse=", "), fill=TRUE)
    }) %>% invisible()
  }

  range_mismatch <- lapply(eem_list, function(eem){
    #eem <- eem_list[[1]]
    if(eem$sample %in% colnames(absorbance)){
      ar <- absorbance[c("wavelength",eem$sample)] %>%
        na.omit() %>%
        .$wavelength %>%
        range()
      er <- eem$em %>% range()
      if(er[1] < ar[1] | er[2] > ar[2]){
        cat(fill=TRUE)
        cat(eem$sample,": absorbance wavelength range is smaller than emission wavelength range, inner-filter effect correction is impossible!",fill=TRUE)
        problem <- TRUE
        return(eem$sample)
      }
    }
  }) %>% invisible()

  #metadata <- meta
  #metacolumns <- c("dilution")
  if(!is.null(metadata)){
    dupls <- rownames(metadata) %>%
      table() %>%
      data.frame() %>%
      filter(Freq > 1) %>%
      .[,1] %>%
      as.character()
    if(length(dupls) > 0) {
      cat(fill=TRUE)
      cat("The following sample names were duplicate in metadata:",paste0(dupls,collapse=", "),fill=TRUE)
    }
    if(!is.null(metacolumns)){
      problem <- lapply(metacolumns,function(col){
        problem <- FALSE
        if(col %in% colnames(metadata)){
          #col <- metacolumns[1]
          #col <- "dilution"
          #cat(col)
          #metadata[col] >= 0
          valid <- rownames(metadata)[metadata[col] >= 0]# %>%
          missing <- eem_names(eem_list)[!eem_names(eem_list) %in% valid]
          if(length(missing) > 0){
            cat(fill=TRUE)
            cat("Metadata column",col,"misses data for samples:",paste0(missing,collapse=", "),fill=TRUE)
            problem <- TRUE
          }
          missing <- valid[!valid %in% eem_names(eem_list)]
          if(length(missing) > 0){
            cat(fill=TRUE)
            cat("Metadata column",col,"contains additional data for samples with no EEM data present:",paste0(missing,collapse=", "),fill=TRUE)
          }
        } else cat("Column",col,"was not found in the table.", fill=TRUE)
        problem
      }) %>%
        unlist() %>%
        any() %>%
        `|`(problem)
    } else cat("No metadata was checked. No columns were supplied.",fill=TRUE)
  } else cat("No metadata was checked. Table was missing.",fill=TRUE)

  if(problem & error) stop("Please read the messages above carefully and correct the problems before continuing the analysis!")
  invisible(problem)
}

#' Create table that contains sample names and locations of files.
#'
#' @description You can use this table as an overview of your files and/or as a template for creating a metadata table.
#'
#' @param eem_list eemlist
#' @param absorbance data frame with absorbance data
#'
#' @return data frame
#' @export
#'
#' @import dplyr
#' @import tidyr
#'
#' @examples
#' data(eem_list)
#' data(abs_data)
#'
#' eem_metatemplate(eem_list,abs_data)
eem_metatemplate <- function(eem_list = NULL, absorbance = NULL){
  if(!is.null(eem_list)){
    t1 <- data.frame(sample = eem_names(eem_list), eem_location = lapply(eem_list,function(eem) eem$location) %>% unlist(),stringsAsFactors = FALSE)
  }
  if(!is.null(absorbance)){
    loc <- attr(absorbance,"location")
    t2 <- data.frame(sample = colnames(absorbance) %>% .[. != "wavelength"], abs_location = ifelse(is.null(loc),NA,loc),stringsAsFactors = FALSE)
    if(exists("t1")) t1 <- full_join(t1,t2,by="sample") else t1 <- t2
  }
  t1
}

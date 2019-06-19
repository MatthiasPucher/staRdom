#' Importer function for generic csv files to be used with eem_read().
#'
#' @description This function can be used to import generic csv files containing EEM data using \code{\link[eemR]{eem_read}}. Excitation wavelengths are assumed column-wise and emission wavelengths row-wise. If your data is arranged the other way round, please use \code{\link{eem_csv2}}
#'
#' @param file path to file passed from eem_read
#'
#' @return list with EEM data
#' @export
#'
#' @examples
#' eems <- system.file("extdata/EEMs",package="staRdom")
#' eem_list <- eem_read(eems, recursive = TRUE, import_function = eem_csv)
#'
#' eem_list
eem_csv <- function(file) {
  .eem_csv(file)
}

#' Importer function for generic csv files to be used with eem_read().
#'
#' @description This function can be used to import generic csv files containing EEM data using \code{\link[eemR]{eem_read}}. Excitation wavelengths are assumed row-wise and emission wavelengths column-wise If your data is arranged the other way round, please use \code{\link{eem_csv}}
#'
#' @param file path to file passed from eem_read
#'
#' @return list with EEM data
#' @export
#'
#' @examples
#' \donttest{
#' ## no example data provided with the package
#' ## below is an example how this could like like
#' eems <- "C:/some/path/to/eem.csv"
#' eem_list <- eem_read(eems, recursive = TRUE, import_function = eem_csv2)
#'
#' eem_list
#' }
eem_csv2 <- function(file) {
  .eem_csv(file, col ="em")
}

#' Import EEMs from generic csv files.
#'
#' @param file path to file
#' @param col either "ex" or "em", whatever wavelength is arranged in columns
#'
#' @return list with EEM data
#' @importFrom data.table fread
#'
.eem_csv <- function(file, col = "ex"){
    x <- fread(file, header = TRUE)
    if(col == "ex"){
      ex <- colnames(x)[-1] %>% as.numeric()
      em <- x[[1]]
    } else if(col == "em"){
      em <- colnames(x)[-1] %>% as.numeric()
      ex <- x[[1]]
    }
    x <- x[,-1] %>% as.matrix() %>% unname()
    x <- x[!is.na(em),!is.na(ex)]
    ex <- ex[!is.na(ex)]
    em <- em[!is.na(em)]
    if(col == "em") x <- t(x)

    l <- list(
      file = file,
      x = x,
      em = em,
      ex = ex
    )

    return(l)
}

#' Importer function for Hitachi F-7000 txt files to be used with eem_read().
#'
#' @description This function can be used to import txt files from Hitachi F-7000 containing EEM data using \code{\link[eemR]{eem_read}}.
#'
#' @param file path to file passed from eem_read
#'
#' @return list with EEM data
#' @export
#'
#' @importFrom readr read_lines
#' @importFrom stringr str_split
#'
#' @examples
#' \donttest{
#' ## no example data provided with the package
#' ## below is an example how this could like like
#' eems <- "C:/some/path/to/hitachi.TXT"
#' eem_list <- eem_read(eems, recursive = TRUE, import_function = eem_hitachi)
#'
#' eem_list
#' }
eem_hitachi <- function(file) {
  data <- read_lines(file)

  data <- str_split(data, "\t")

  data <- do.call(rbind.data.frame, data)

  header <- data[1:which(grepl("^Data Points",data[,1],ignore.case = TRUE)),1:2]

  data <- data[which(grepl("^Data Points",data[,1],ignore.case = TRUE))+1:nrow(data),]
  data <- as.data.frame(lapply(data,function(x) if(is.character(x)|is.factor(x)) gsub(",",".",x) else x))

  em <- as.vector(na.omit(data[2:nrow(data),1]))
  em <- as.numeric(em)

  ex <- as.vector(unlist(data[1,2:ncol(data)]))
  ex <- as.numeric(ex)

  data <- lapply(data, function(x) as.numeric(as.character(x)))

  data <- do.call(rbind.data.frame, data)

  data <- data.frame(data[2:nrow(data),2:(which(is.na(data[1,]))[2]-1)])

  eem <- matrix(unlist(t(data)), nrow = length(em), byrow = FALSE)

  l <- list(
    file = file,
    x = eem,
    em = em,
    ex = ex
  )

  return(l)
}

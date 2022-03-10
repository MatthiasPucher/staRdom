## ---- message=FALSE, warning=FALSE, include=FALSE-----------------------------
### BEWARE! ###
# This Rmd file was created to supply a tutrial for staRdom.
# Several steps in here are not part of the actual analysis but to create a nice looking and smoothly created
# tutorial file. Starting your PARAFAC analysis from this Rmd file is therefore not adviced. If you do so, you
# have to distinguish between the steps necessary for the analysis and those necessary for a nice tutorial.
# If you encounter troubles or you do not know, what this means please start with a fresh, plain R script, read
# the turial and use the code provided there!

library(knitr)
library(kableExtra)
cores <- 2

## ----message=FALSE, warning=TRUE----------------------------------------------
library(dplyr)
library(tidyr)

## ----eval=FALSE, include=TRUE-------------------------------------------------
#  install.packages("staRdom")

## ----eval=FALSE, include=TRUE-------------------------------------------------
#  install.packages("devtools") # Run this only, if devtools is not installed already.
#  devtools::install_github("MatthiasPucher/staRdom")

## ----eval=TRUE, message=FALSE, warning=FALSE, include=TRUE--------------------
library("staRdom")

## ----eval=FALSE, include=TRUE-------------------------------------------------
#  cores <- detectCores(logical = FALSE)

## ----eval=TRUE, include=TRUE--------------------------------------------------
folder <- system.file("extdata/EEMs/", package = "staRdom") # folder containing example EEMs
eem_list <- eem_read(folder, recursive = TRUE, import_function = eem_csv) # in case you use your own data, just replace folder by a path. e.g. "C:/folder/another folder" and change import_function according to instrument.

## ----eval=FALSE, include=TRUE-------------------------------------------------
#  eem_list <- eem_read(folder, import_function = "cary")

## ----eval=TRUE, fig.height=5, fig.width=6, message=FALSE, warning=FALSE, include=TRUE, paged.print=TRUE----
eem_overview_plot(eem_list, spp=9, contour = TRUE)

## ----eval=TRUE, include=TRUE--------------------------------------------------
absorbance_path = system.file("extdata/absorbance", package = "staRdom") # load example data, set a path without using system.file to use your own data e.g. "C:/folder/another folder"

## ----eval=TRUE, include=TRUE--------------------------------------------------
absorbance <- absorbance_read(absorbance_path, cores = cores) # load csv or txt tables in folder

## ----eval=TRUE, include=TRUE--------------------------------------------------
metatable <- system.file("extdata/metatable_dreem.csv",package = "staRdom") # path to example data, can be replaced by a path to your own data
meta <- read.table(metatable, header = TRUE, sep = ",", dec = ".", row.names = 1) # load data

## ----eval=FALSE, include=TRUE-------------------------------------------------
#  eem_metatemplate(eem_list, absorbance) %>%
#    write.csv(file="metatable.csv", row.names = FALSE)

## ----eval=TRUE, include=TRUE--------------------------------------------------
problem <- eem_checkdata(eem_list,absorbance,meta,metacolumns = c("dilution"),error=FALSE)

## ----eval=TRUE, include=TRUE--------------------------------------------------
eem_list <- eem_name_replace(eem_list,c("\\(FD3\\)"),c(""))

## ----eval=TRUE, include=TRUE--------------------------------------------------
absorbance <- abs_blcor(absorbance,wlrange = c(680,700))

## ----eval=TRUE, include=TRUE--------------------------------------------------
excorfile <- system.file("extdata/CorrectionFiles/xc06se06n.csv",package="staRdom")
Excor <- data.table::fread(excorfile)

emcorfile <- system.file("extdata/CorrectionFiles/mcorrs_4nm.csv",package="staRdom")
Emcor <- data.table::fread(emcorfile)

# adjust range of EEMs to cover correction vectors
eem_list <- eem_range(eem_list,ex = range(Excor[,1]), em = range(Emcor[,1]))

eem_list <- eem_spectral_cor(eem_list,Excor,Emcor)

## ----eval=TRUE, include=TRUE--------------------------------------------------
# extending and interpolation data
eem_list <- eem_extend2largest(eem_list, interpolation = 1, extend = FALSE, cores = cores)

# blank subtraction
eem_list <- eem_remove_blank(eem_list)

## ----eval=TRUE, fig.width=6, fig.height = 5, message=FALSE, warning=FALSE, include=TRUE, paged.print=TRUE----
eem_overview_plot(eem_list, spp=9, contour = TRUE)

## ----eval=TRUE, include=TRUE--------------------------------------------------
eem_list <- eem_ife_correction(eem_list,absorbance, cuvl = 5)

## ----eval=TRUE, fig.width=6, fig.height = 5, message=FALSE, warning=FALSE, include=TRUE----
eem_overview_plot(eem_list, spp=9, contour = TRUE)

## ----eval=TRUE, include=TRUE--------------------------------------------------
eem_list <- eem_raman_normalisation2(eem_list, blank = "blank")

## ----eval=TRUE, fig.width=6, fig.height = 5, message=FALSE, warning=FALSE, include=TRUE----
eem_overview_plot(eem_list, spp=9, contour = TRUE)

## ----eval=TRUE, include=TRUE--------------------------------------------------
eem_list <- eem_extract(eem_list, c("nano", "miliq", "milliq", "mq", "blank"),ignore_case = TRUE)

absorbance <- dplyr::select(absorbance, -matches("nano|miliq|milliq|mq|blank", ignore.case = TRUE))

## ----eval=TRUE, include=TRUE--------------------------------------------------
remove_scatter <- c(TRUE, TRUE, TRUE, TRUE)
remove_scatter_width <- c(15,15,15,15)

eem_list <- eem_rem_scat(eem_list, remove_scatter = remove_scatter, remove_scatter_width = remove_scatter_width)

## ----eval=TRUE, fig.width=6, fig.height = 5, message=FALSE, warning=FALSE, include=TRUE----
eem_overview_plot(eem_list, spp=9, contour = TRUE)

## ----eval=TRUE, include=TRUE--------------------------------------------------
eem_list <- eem_interp(eem_list, cores = cores, type = 1, extend = FALSE)

## ----eval=TRUE, fig.width=6, fig.height = 5, message=FALSE, warning=FALSE, include=TRUE----
eem_overview_plot(eem_list, spp=9, contour = TRUE)

## ----eval=TRUE, include=TRUE--------------------------------------------------
dil_data <- meta["dilution"]

eem_list <- eem_dilution(eem_list,dil_data)

## ----eval=FALSE, fig.width=7, message=FALSE, warning=FALSE, include=TRUE, paged.print=TRUE----
#  eem_overview_plot(eem_list, spp=9) # plot spared, due to no dilution it looks like the previous plot.

## ----eval=TRUE, include=TRUE--------------------------------------------------
eem4peaks <- eem_smooth(eem_list, n = 4, cores = cores)

## ----eval=FALSE, fig.width=7, message=FALSE, warning=FALSE, include=FALSE, paged.print=TRUE----
#  eem_overview_plot(eem4peaks, spp=6)

## -----------------------------------------------------------------------------
summary(eem_list)

## ----eval=FALSE, include=FALSE------------------------------------------------
#  eem_list %>%
#    summary() %>%
#    kable(format = "latex", booktabs = T) %>%
#    kable_styling(bootstrap_options = c("striped", "hover"), full_width = FALSE) %>%
#    #kable_styling(font_size = 5) %>%
#    row_spec(0, angle = -45) #%>%

## ----eval=TRUE, warning=FALSE, include=TRUE-----------------------------------

bix <- eem_biological_index(eem4peaks)
coble_peaks <- eem_coble_peaks(eem4peaks)
fi <- eem_fluorescence_index(eem4peaks)
hix <- eem_humification_index(eem4peaks, scale = TRUE)

indices_peaks <- bix %>%
  full_join(coble_peaks, by = "sample") %>%
  full_join(fi, by = "sample") %>%
  full_join(hix, by = "sample")

indices_peaks

## ----eval=FALSE, include=FALSE, echo=FALSE------------------------------------
#  kable(indices_peaks %>% mutate_if(is.numeric,prettyNum,digits = 2,format="fg"),format = "latex", booktabs = T) %>%
#    kable_styling(bootstrap_options = c("striped", "hover"), full_width = FALSE) %>%
#    #kable_styling(font_size = 5) %>%
#    row_spec(0, angle = 45) #%>%

## ----eval=TRUE, include=TRUE--------------------------------------------------
slope_parms <- abs_parms(absorbance, cuvl = 1, cores = cores)
slope_parms

## ----eval=FALSE, include=FALSE, echo=FALSE------------------------------------
#  kable(slope_parms %>% mutate_if(is.numeric,prettyNum,digits = 2,format="fg"),format = "latex", booktabs = T) %>%
#    #kable_styling(bootstrap_options = c("striped", "hover"), full_width = FALSE) %>%
#    #kable_styling(font_size = 5) %>%
#    row_spec(0, angle = 45) #%>%

## ----eval=FALSE, include=TRUE-------------------------------------------------
#  dreem_raw <- tempfile()
#  download.file("http://models.life.ku.dk/sites/default/files/drEEM_dataset.zip",dreem_raw)
#  dreem_data <- unz(dreem_raw, filename="Backup/PortSurveyData_corrected.mat", open = "rb") %>%
#    R.matlab::readMat()
#  unlink(dreem_raw)
#  
#  eem_list <- lapply(dreem_data$filelist.eem, function(file){
#    #file <- dreem_data$filelist.eem[1]
#    n <- which(dreem_data$filelist.eem == file)
#    file <- file %>%
#      gsub("^\\s+|\\s+$", "", .) %>% # trim white spaces in filenames
#      sub(pattern = "(.*)\\..*$", replacement = "\\1", .) # remove file extension from sample name
#    eem <- list(file = paste0("drEEM/dataset/",file),sample = file,x = dreem_data$XcRU[n,,] %>% as.matrix(),ex = dreem_data$Ex %>% as.vector(), em = dreem_data$Em.in %>% as.vector(), location = "drEEM/dataset/")
#    class(eem) <- "eem"
#    attr(eem, "is_blank_corrected") <- TRUE
#    attr(eem, "is_scatter_corrected") <- FALSE
#    attr(eem, "is_ife_corrected") <- TRUE
#    attr(eem, "is_raman_normalized") <- TRUE
#    attr(eem, "manufacturer") <- "unknown"
#    eem
#  }) %>%
#    `class<-`("eemlist")
#  
#  # add sample name suffix, R has sometimes troubles, when sample names start with a number.
#  eem_names(eem_list) <- paste0("d",eem_names(eem_list))

## ----eval=FALSE, include=TRUE-------------------------------------------------
#  ol <- function(x){x==("bl") | x == "0A"}
#  extract <- dreem_data$sites %>% unlist() %>% ol() %>% which()
#  eem_list <- eem_list %>% eem_extract(extract)

## ----message=FALSE, warning=FALSE, include=FALSE------------------------------
data(eem_list) # load example from staRdom package, this is just necessary for the actual tutorial creation. Remove this line, if you downloaded the drEEM dataset right above and want to use that.
eem_ex <- eem_extract(eem_list,sample ="^d667sf$",keep=TRUE)

## ----include = TRUE, eval = TRUE----------------------------------------------
eem_list <- eem_rem_scat(eem_list, remove_scatter = c(TRUE, TRUE, TRUE, TRUE), remove_scatter_width = c(15,15,18,19), interpolation = FALSE, cores = cores)

## ----include = FALSE, eval = TRUE---------------------------------------------
eem_ex <- eem_ex %>% 
  eem_extract(sample = "^d667sf$", keep = TRUE) %>%
  eem_rem_scat(eem_ex, remove_scatter = c(TRUE, TRUE, TRUE, TRUE), remove_scatter_width = c(15,15,18,19), interpolation = FALSE, cores = cores) %>%
  `eem_names<-`("d667sf_2_scatter") %>%
  eem_bind(eem_ex,.)

## ----include=TRUE,eval=FALSE--------------------------------------------------
#  eem_list <- eem_import_dir(dir)

## ----include=TRUE, eval = FALSE-----------------------------------------------
#  eem_list %>%
#    eem_extract(sample = "^d667sf$", keep = TRUE) %>%
#    ggeem(contour = TRUE)

## ----include=FALSE, eval = TRUE-----------------------------------------------
eem_ex %>% 
  ggeem(contour = TRUE)

## ----eval=TRUE----------------------------------------------------------------
eem_list <- eem_list %>% eem_range(ex = c(250,Inf), em = c(0,580))

## ----message=FALSE, warning=FALSE, include=FALSE------------------------------
eem_ex <- eem_ex %>% 
  eem_extract(sample = "^d667sf_2_scatter$", keep = TRUE) %>%
  eem_range(ex = c(250,Inf), em = c(0,580)) %>%
  `eem_names<-`("d667sf_3_cut") %>%
  eem_bind(eem_ex,.)

## ----eval = TRUE, message = FALSE, warning = FALSE, include = FALSE-----------
eem_ex <- eem_ex %>%
  eem_extract(sample = "^d667sf_3_cut$", keep = TRUE) %>%
  eem_setNA(sample = 1, ex = 345:350, interpolate = FALSE) %>% # sample 1 in staRdom data is sample 176 in drEEM data!
  eem_setNA(em = 560:576, ex = 280:295, interpolate = FALSE) %>%
  `eem_names<-`("d667sf_4_rem_noise") %>%
  eem_bind(eem_ex,.) 

## ----include=TRUE, eval=TRUE--------------------------------------------------
eem_list <- eem_list %>%
  eem_setNA(sample = 176, ex = 345:350, interpolate = FALSE) %>%
  eem_setNA(em = 560:576, ex = 280:295, interpolate = FALSE)

## ----eval = TRUE, include = TRUE----------------------------------------------
eem_list <- eem_interp(eem_list, type = 1, extend = FALSE, cores = cores)

## ----message=FALSE, warning=FALSE, include=FALSE------------------------------
eem_ex <- eem_ex %>% 
  eem_extract(sample = "^d667sf_4_rem_noise$", keep = TRUE) %>%
  eem_interp(type = 1, extend = FALSE, cores = cores) %>%
  `eem_names<-`("d667sf_5_interp") %>%
  eem_bind(eem_ex,.)

## ----echo=FALSE, message=FALSE, warning=FALSE, fig.width=6, fig.height = 4----
ggeem(eem_ex, contour = TRUE)

## ----include=TRUE, eval = TRUE------------------------------------------------
data(pf_models)

## ----include=TRUE-------------------------------------------------------------
# minimum and maximum of numbers of components
dim_min <- 3
dim_max <- 7

## ----eval=FALSE,include=TRUE--------------------------------------------------
#  nstart <- 25 # number of similar models from which best is chosen
#  maxit = 5000 # maximum number of iterations in PARAFAC analysis
#  ctol <- 10^-6 # tolerance in PARAFAC analysis
#  
#  # calculating PARAFAC models, one for each number of components
#  pf1 <- eem_parafac(eem_list, comps = seq(dim_min,dim_max), normalise = FALSE, const = c("uncons", "uncons", "uncons"), maxit = maxit, nstart = nstart, ctol = ctol, cores = cores)
#  
#  # same model but using non-negative constraints
#  pf1n <- eem_parafac(eem_list, comps = seq(dim_min,dim_max), normalise = FALSE, const = c("nonneg", "nonneg", "nonneg"), maxit = maxit, nstart = nstart, ctol = ctol, cores = cores)
#  
#  # rescale B and C modes to a maximum fluorescence of 1 for each component
#  pf1 <- lapply(pf1, eempf_rescaleBC, newscale = "Fmax")
#  pf1n <- lapply(pf1n, eempf_rescaleBC, newscale = "Fmax")

## ----eval=FALSE, include=TRUE, fig.width=6, fig.height=5----------------------
#  # This plot is not shown, because the components violate the assumptions for fluorescence peaks (negative fluorescence). Please try, if you are interested.
#  eempf_compare(pf1, contour = TRUE)

## ----eval=TRUE, include=TRUE, fig.width=6, fig.height=5-----------------------
eempf_compare(pf1n, contour = TRUE)

## ----eval=TRUE, include=TRUE, fig.width=7-------------------------------------
# check for correlation between components table
eempf_cortable(pf1n[[4]], normalisation = FALSE)

## ----eval=TRUE, include=TRUE, fig.width=6, fig.height=5-----------------------
eempf_corplot(pf1n[[4]], progress = FALSE, normalisation = FALSE)

## ----eval=FALSE,include=TRUE--------------------------------------------------
#  pf2 <- eem_parafac(eem_list, comps = seq(dim_min,dim_max), normalise = TRUE, const = c("nonneg", "nonneg", "nonneg"), maxit = maxit, nstart = nstart, ctol = ctol, cores = cores)
#  
#  # rescale B and C modes
#  pf2 <- lapply(pf2, eempf_rescaleBC, newscale = "Fmax")

## ----eval=TRUE, include=TRUE, fig.width=6, fig.height=5-----------------------
# eempf_compare(pf2, contour = TRUE) # use this to show the same plot as above
# for now, we are happy with just the components
eempf_plot_comps(pf2, contour = TRUE, type = 1)

## ----fig.width=6--------------------------------------------------------------
# calculate leverage
cpl <- eempf_leverage(pf2[[4]])

# plot leverage (nice plot)
eempf_leverage_plot(cpl,qlabel=0.1)

## ----eval=FALSE, include=TRUE, fig.width=7------------------------------------
#  # plot leverage, not so nice plot but interactive to select what to exclude
#  # saved in exclude, can be used to start over again with eem_list_ex <- eem_list %>% eem_exclude(exclude) above
#  exclude <- eempf_leverage_ident(cpl,qlabel=0.1)

## ----eval=FALSE, include=TRUE-------------------------------------------------
#  # samples, excitation and emission wavelengths to exclude, makes sense after calculation of leverage
#  exclude <- list("ex" = c(),
#                  "em" = c(),
#                  "sample" = c("dsfb676psp","dsgb447wt")
#  )
#  
#  # exclude outliers if neccessary. if so, restart analysis
#  eem_list_ex <- eem_exclude(eem_list, exclude)

## ----eval=FALSE, include=TRUE-------------------------------------------------
#  pf3 <- eem_parafac(eem_list_ex, comps = seq(dim_min,dim_max), normalise = TRUE, maxit = maxit, nstart = nstart, ctol = ctol, cores = cores)
#  pf3 <- lapply(pf3, eempf_rescaleBC, newscale = "Fmax")

## ----eval=TRUE, include=TRUE, fig.width=6, fig.height=5-----------------------
eempf_plot_comps(pf3, contour = TRUE, type = 1)

## ----eval=TRUE, include=TRUE, fig.width=7-------------------------------------
eempf_leverage_plot(eempf_leverage(pf3[[4]]),qlabel=0.1)

## ----eval=FALSE, include=TRUE, fig.width=6, fig.height=5----------------------
#  eempf_residuals_plot(pf3[[4]], eem_list, residuals_only = TRUE, select = c("d0680sfK", "d1266sf", "d1268sfK", "d1543sfK", "dsfb676psp", "dsgb447wt"), spp = 6, cores = cores, contour = TRUE)

## ----eval=TRUE, echo=FALSE, message=FALSE, warning=FALSE, fig.width=6, fig.height=5----
data("eem_list_outliers")
eem_list %>%
  eem_extract(5:15) %>%
  eem_bind(eem_list_outliers) %>%
  eem_red2smallest() %>%
  eempf_residuals_plot(pf3[[4]], ., residuals_only = TRUE, spp = 6, cores = cores, contour = TRUE) %>%
  invisible()

## ----eval=FALSE, include=TRUE-------------------------------------------------
#  ctol <- 10^-8 # decrease tolerance in PARAFAC analysis
#  nstart = 20 # number of random starts
#  maxit = 10000 # increase number of maximum interations
#  
#  pf4 <- eem_parafac(eem_list_ex, comps = 6, normalise = TRUE, const = c("nonneg", "nonneg", "nonneg"), maxit = maxit, nstart = nstart, ctol = ctol, output = "all", cores = cores, strictly_converging = TRUE)
#  
#  pf4 <- lapply(pf4, eempf_rescaleBC, newscale = "Fmax")

## ----eval=TRUE, include=TRUE, fig.width=4, fig.height=6-----------------------
eempf_convergence(pf4[[1]])

## ----eval=FALSE, include=TRUE, fig.width=4, fig.height=6----------------------
#  # just one model, not really a need to compare
#  eempf_compare(pf4, contour = TRUE)

## ----eval=TRUE, include=TRUE, fig.width=6-------------------------------------
eempf_leverage_plot(eempf_leverage(pf4[[1]])) # [[1]] means the 4th model in the list, 6 component model in that case

## ----eval=TRUE, include=TRUE, fig.width=6, fig.height=5-----------------------
eempf_corplot(pf4[[1]], progress = FALSE)

## ----eval=FALSE, include=TRUE-------------------------------------------------
#  # calculating a rough model, nstart is high (100) but ctol is 2 magnitudes larger or at least 0.01
#  pf5 <- eem_parafac(eem_list, comps = 6, normalise = TRUE, const = c("nonneg", "nonneg", "nonneg"), maxit = maxit, nstart = 100, ctol = min(ctol*100,0.01), cores = cores)
#  
#  # plot is not shown
#  ggeem(pf5[[1]], contour = TRUE)
#  
#  nstart <- 5
#  pf4 <- eem_parafac(eem_list_ex, comps = 6, normalise = TRUE, const = c("nonneg", "nonneg", "nonneg"), maxit = maxit, nstart = nstart, ctol = ctol, cores = cores, Bstart = pf5[[1]]$B, Cstart = pf5[[1]]$C)
#  
#  pf4 <- lapply(pf4, eempf_rescaleBC, newscale = "Fmax")
#  
#  # plot is not shown
#  ggeem(pf4[[1]], contour = TRUE)

## ----eval=TRUE, include=TRUE, fig.width=6, fig.height=5-----------------------
eempf_comp_load_plot(pf4[[1]], contour = TRUE)

# eempf_plot_comps(pf4[1], type = 2) # this function can be used to view the B- and C-modes

## ----eval=TRUE, fig.height=7, fig.width=6, warning=FALSE, include=TRUE--------
# plot components in each sample, residual and whole sample
eempf_residuals_plot(pf4[[1]], eem_list, select = eem_names(eem_list)[10:14], cores = cores, contour = TRUE)

## ----eval=FALSE, include=TRUE-------------------------------------------------
#  #calculate split_half analysis
#  sh <- splithalf(eem_list_ex, 6, normalise = TRUE, rand = FALSE, cores = cores, nstart = nstart, strictly_converging = TRUE, maxit = maxit, ctol = ctol)

## ----eval=TRUE, include=TRUE, fig.width=7-------------------------------------
data(sh)

## ----eval=TRUE, include=TRUE, fig.width=6-------------------------------------
splithalf_plot(sh)

# you can also use the already known plots from eempf_compare
# sh %>% unlist(recursive = FALSE) %>% eempf_compare()

## ----eval=FALSE, include=TRUE-------------------------------------------------
#  sh_r <- splithalf(eem_list_ex, 6, normalise = TRUE, rand = TRUE, cores = cores, nstart = nstart, maxit = maxit, ctol = ctol)

## ----eval=TRUE, include=TRUE, fig.width=7-------------------------------------
tcc_sh_table <- splithalf_tcc(sh)

tcc_sh_table

## ----include = TRUE, eval = TRUE----------------------------------------------
pf4_wOutliers <- A_missing(eem_list, pfmodel = pf4[[1]], cores = cores)

## ----eval=FALSE, include=TRUE, fig.width=7------------------------------------
#  corcondia <- eempf_corcondia(pf4[[1]], eem_list_ex)

## ----eval=FALSE, include=TRUE, fig.width=7------------------------------------
#  eemqual <- eempf_eemqual(pf4[[1]], eem_list_ex, sh, cores = cores)

## ----eval=FALSE, include=TRUE, fig.width=7------------------------------------
#  varimp <- eempf_varimp(pf4[[1]], eem_list_ex, cores = cores)

## ----fig.width=6--------------------------------------------------------------
# get current model names (none set so far!)
names(pf3)
# set new model names, number of models must be equal to number of names
names(pf3) <- c("3 components", "4 components xy","5 components no outliers","6 components","7 components")
names(pf3)

# get current component names
eempf_comp_names(pf4)
# set new model names, number of models must be equal to number of names
eempf_comp_names(pf4) <- c("A4","B4","C4","D4","E4","F4")

# in case of more than one model(e.g. pf3):
eempf_comp_names(pf3) <- list(c("A1","B1","C1"), # names for 1st model
                                       c("humic","T2","whatever","peak"),
                                       c("rose","peter","frank","dwight","susan"),
                                       c("A4","B4","C4","D4","E4","F4"),
                                       c("A5","B5","C5","D5","E5","F5","G5") # names for 5th model
)

pf4[[1]] %>%
  ggeem(contour = TRUE)

## ----eval=FALSE, include = TRUE-----------------------------------------------
#  eempf_openfluor(pf4[[1]], file = "my_model_openfluor.txt")

## ----eval=FALSE, include = TRUE-----------------------------------------------
#  eempf_report(pf4[[1]], export = "parafac_report.html", eem_list = eem_list, shmodel = sh, performance = TRUE)

## ----eval=FALSE, include=TRUE, fig.width=7------------------------------------
#  pf4 <- eem_parafac(eem_list_ex, comps = 6, normalise = TRUE, const = c("nonneg", "nonneg", "nonneg"), maxit = maxit, nstart = nstart, ctol = ctol, output = "all", cores = cores)
#  
#  ssccheck <- eempf_ssccheck(pf4[[1]]$models, best = 3, cores = cores) # best 3 models are shown
#  
#  eempf_plot_ssccheck(ssccheck)


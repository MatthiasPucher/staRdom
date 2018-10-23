## ---- message=FALSE, warning=FALSE, include=FALSE------------------------
library(knitcitations)
library(dplyr)
library(tidyr)
cleanbib()
options("citation_format" = "pandoc")
bibliography() #style="apalike"
library(staRdom)
library(knitr)
library(kableExtra)
cores <- 2

## ----eval=FALSE, include=TRUE--------------------------------------------
#  cores <- parallel::detectCores()/2

## ----eval=FALSE, include=FALSE-------------------------------------------
#  data(eem_list) # load example data

## ----eval=TRUE, include=TRUE---------------------------------------------
folder <- system.file("extdata/cary/scans_day_1", package = "eemR") # load example data
eem_list <- eem_read(folder)

## ----eval=TRUE, fig.width=7, message=FALSE, warning=FALSE, include=TRUE, paged.print=TRUE----
ggeem(eem_list)

## ----eval=TRUE, include=TRUE---------------------------------------------
data(absorbance) # load example data

## ----eval=TRUE, include=TRUE---------------------------------------------
meta <- read.table(system.file("extdata/metatable_eemR.csv",package = "staRdom"), header = TRUE, sep = " ", dec = ".", row.names = 1) # load example data

## ----eval=TRUE, include=TRUE---------------------------------------------
problem <- eem_checkdata(eem_list,absorbance,meta,metacolumns = c("dilution"),error=FALSE)

## ----eval=TRUE, include=TRUE---------------------------------------------
eem_list <- eem_name_replace(eem_list,c("\\(FD3\\)"),c(""))

## ----eval=TRUE, include=TRUE---------------------------------------------
absorbance <- abs_blcor(absorbance,wlrange = c(680,700))

## ----eval=TRUE, include=TRUE---------------------------------------------
eem_list <- eem_remove_blank(eem_list)

## ----eval=TRUE, fig.width=7, message=FALSE, warning=FALSE, include=TRUE, paged.print=TRUE----
ggeem(eem_list)

## ----eval=TRUE, include=TRUE---------------------------------------------
eem_list <- eem_ife_correction(eem_list,absorbance, cuvl = 5)

## ----eval=TRUE, fig.width=7, message=FALSE, warning=FALSE, include=TRUE, paged.print=TRUE----
ggeem(eem_list)

## ----eval=TRUE, include=TRUE---------------------------------------------
eem_list <- eem_raman_normalisation2(eem_list, blank = "blank")

## ----eval=TRUE, fig.width=7, message=FALSE, warning=FALSE, include=TRUE, paged.print=TRUE----
ggeem(eem_list)

## ----eval=TRUE, include=TRUE---------------------------------------------
eem_list <- eem_extract(eem_list, c("nano", "miliq", "milliq", "mq", "blank"),ignore_case = TRUE)

absorbance <- select(absorbance, -matches("nano|miliq|milliq|mq|blank", ignore.case = TRUE))

## ----eval=TRUE, include=TRUE---------------------------------------------
remove_scatter <- c(TRUE, TRUE, TRUE, TRUE)
remove_scatter_width <- c(15,15,15,15)

eem_list <- eem_rem_scat(eem_list, remove_scatter = remove_scatter, remove_scatter_width = remove_scatter_width)

## ----eval=TRUE, fig.width=7, message=FALSE, warning=FALSE, include=TRUE, paged.print=TRUE----
ggeem(eem_list)

## ----eval=TRUE, include=TRUE---------------------------------------------
eem_list <- eem_interp(eem_list, cores = cores, type = 3)

## ----eval=TRUE, fig.width=7, message=FALSE, warning=FALSE, include=TRUE, paged.print=TRUE----
ggeem(eem_list)

## ----eval=TRUE, include=TRUE---------------------------------------------
dil_data <- meta["dilution"]

eem_list <- eem_dilution(eem_list,dil_data)

## ----eval=TRUE, fig.width=7, message=FALSE, warning=FALSE, include=TRUE, paged.print=TRUE----
ggeem(eem_list)

## ----eval=TRUE, include=TRUE---------------------------------------------
eem4peaks <- eem_smooth(eem_list, n = 4)

## ----eval=TRUE, fig.width=7, message=FALSE, warning=FALSE, include=TRUE, paged.print=TRUE----
ggeem(eem4peaks)

## ------------------------------------------------------------------------
summary(eem_list)

## ----eval=FALSE, include=FALSE-------------------------------------------
#  eem_list %>%
#      summary() %>%
#      kable(format = "latex", booktabs = T) %>%
#      kable_styling(bootstrap_options = c("striped", "hover"), full_width = FALSE) %>%
#      #kable_styling(font_size = 5) %>%
#      row_spec(0, angle = -45) #%>%

## ----eval=TRUE, warning=FALSE, include=TRUE------------------------------

bix <- eem_biological_index(eem4peaks)
coble_peaks <- eem_coble_peaks(eem4peaks)
fi <- eem_fluorescence_index(eem4peaks)
hix <- eem_humification_index(eem4peaks, scale = TRUE)

indices_peaks <- bix %>%
  full_join(coble_peaks, by = "sample") %>%
  full_join(fi, by = "sample") %>%
  full_join(hix, by = "sample")

indices_peaks

## ----eval=FALSE, include=FALSE, echo=FALSE-------------------------------
#    kable(indices_peaks %>% mutate_if(is.numeric,prettyNum,digits = 2,format="fg"),format = "latex", booktabs = T) %>%
#      kable_styling(bootstrap_options = c("striped", "hover"), full_width = FALSE) %>%
#      #kable_styling(font_size = 5) %>%
#      row_spec(0, angle = 45) #%>%

## ----eval=TRUE, include=TRUE---------------------------------------------
slope_parms <- abs_parms(absorbance, cuvl = 1, cores = cores)
slope_parms

## ----eval=FALSE, include=FALSE, echo=FALSE-------------------------------
#    kable(slope_parms %>% mutate_if(is.numeric,prettyNum,digits = 2,format="fg"),format = "latex", booktabs = T) %>%
#      #kable_styling(bootstrap_options = c("striped", "hover"), full_width = FALSE) %>%
#      #kable_styling(font_size = 5) %>%
#      row_spec(0, angle = 45) #%>%

## ----include=TRUE--------------------------------------------------------
data(eem_list) # load example data

## ----include=TRUE--------------------------------------------------------
eem_list <- eem_rem_scat(eem_list, remove_scatter = c(TRUE, TRUE, TRUE, TRUE), remove_scatter_width = c(15,15,18,19), interpolation = FALSE, cores = cores)

## ----include=TRUE,eval=FALSE---------------------------------------------
#  eem_list <- eem_import_dir(dir)

## ----include=TRUE--------------------------------------------------------
eem_list %>% 
  eem_extract(sample = "^667sf$", keep = TRUE) %>%
  ggeem()

## ----eval=FALSE, message=FALSE, warning=FALSE, include=TRUE--------------
#  eem_list <- eem_range(eem_list, ex = c(250,Inf), em = c(0,580))

## ----message=FALSE, warning=FALSE, include=FALSE-------------------------
eem_example <- eem_list %>% 
  eem_extract(sample = "^667sf$", keep = TRUE) %>%
  `eem_names<-`("667sf_2_cut") %>%
  eem_bind(eem_example,.)

## ------------------------------------------------------------------------
eem_list <- eem_list %>%
  eem_setNA(sample = 176, ex = 345:350, interpolate = FALSE) %>%
  eem_setNA(em = 560:576, ex = 280:295, interpolate = FALSE)

## ----message=FALSE, warning=FALSE, include=FALSE-------------------------
eem_example <- eem_list %>% 
  eem_extract(sample = "^667sf$", keep = TRUE) %>%
  `eem_names<-`("667sf_3_rem_noise") %>%
  eem_bind(eem_example,.)

## ------------------------------------------------------------------------
eem_list <- eem_interp(eem_list, type = 3, cores = cores)

## ----message=FALSE, warning=FALSE, include=FALSE-------------------------
eem_example <- eem_list %>% 
  eem_extract(sample = "^667sf$", keep = TRUE) %>%
  `eem_names<-`("667sf_4_interp") %>%
  eem_bind(eem_example,.)

## ----echo=FALSE, message=FALSE, warning=FALSE, fig.width=7---------------
ggeem(eem_example)

## ----include=TRUE--------------------------------------------------------
# minimum and maximum of numbers of components
dim_min <- 3
dim_max <- 7

## ----eval=FALSE,include=TRUE---------------------------------------------
#  nstart <- 16 # number of similar models from which best is chosen
#  maxit = 1000 # maximum number of iterations in PARAFAC analysis
#  ctol <- 10^-5 # tolerance in PARAFAC analysis
#  
#  # calculating PARAFAC models, one for each number of components
#  pf1 <- eem_parafac(eem_list, comps = seq(dim_min,dim_max), normalise = FALSE, const = c("uncons", "uncons", "uncons"), maxit = maxit, nstart = nstart, ctol = ctol, cores = cores)
#  
#  pf1n <- eem_parafac(eem_list, comps = seq(dim_min,dim_max), normalise = FALSE, const = c("nonneg", "nonneg", "nonneg"), maxit = maxit, nstart = nstart, ctol = ctol, cores = cores)
#  
#  # rescale B and C modes
#  pf1 <- lapply(pf1, eempf_rescaleBC, newscale = "Fmax")
#  pf1n <- lapply(pf1n, eempf_rescaleBC, newscale = "Fmax")

## ----include=TRUE--------------------------------------------------------
data(pf_models)

## ----eval=TRUE, include=TRUE, fig.width=7, fig.height=6------------------
eempf_compare(pf1)

## ----eval=TRUE, include=TRUE, fig.width=7, fig.height=6------------------
eempf_compare(pf1n)

## ----eval=TRUE, include=TRUE, fig.width=7--------------------------------
# check for correlation between components table
eempf_cortable(pf1n[[4]])

## ----eval=TRUE, include=TRUE, fig.width=7, fig.height=6------------------
eempf_corplot(pf1n[[4]], progress = FALSE)

## ----eval=FALSE,include=TRUE---------------------------------------------
#  pf2 <- eem_parafac(eem_list, comps = seq(dim_min,dim_max), normalise = TRUE, const = c("nonneg", "nonneg", "nonneg"), maxit = maxit, nstart = nstart, ctol = ctol, cores = cores)
#  
#  # rescale B and C modes
#  pf2 <- lapply(pf2, eempf_rescaleBC, newscale = "Fmax")

## ----eval=TRUE, include=TRUE, fig.width=7, fig.height=6------------------
eempf_compare(pf2)

## ----fig.width=7---------------------------------------------------------
# calculate leverage
cpl <- eempf_leverage(pf2[[4]])
# plot leverage (nice plot)
eempf_leverage_plot(cpl,qlabel=0.1)
# plot leverage, not so nice plot but interactive to select what to exclude
# saved in exclude, can be used to start over again with eem_list_ex <- eem_list %>% eem_exclude(exclude) above
exclude <- eempf_leverage_ident(cpl,qlabel=0.1)

## ----eval=TRUE, include=TRUE---------------------------------------------
# samples, excitation and emission wavelengths to exclude, makes sense after calculation of leverage
exclude <- list("ex" = c(),
                "em" = c(),
                "sample" = c("sfb676psp","sgb447wt")
)

# exclude outliers if neccessary. if so, restart analysis
eem_list_ex <- eem_exclude(eem_list, exclude)

## ----eval=FALSE, include=TRUE--------------------------------------------
#  pf3 <- eem_parafac(eem_list_ex, comps = seq(dim_min,dim_max), normalise = TRUE, maxit = maxit, nstart = nstart, ctol = ctol, cores = cores)
#  pf3 <- lapply(pf3, eempf_rescaleBC, newscale = "Fmax")

## ----eval=TRUE, include=TRUE, fig.width=7, fig.height=6------------------
eempf_compare(pf3)

eempf_leverage_plot(eempf_leverage(pf3[[4]]),qlabel=0.1)

## ----eval=TRUE, include=TRUE, fig.width=7, fig.height=6------------------
eempf_residuals_plot(pf3[[4]], eem_list, residuals_only = TRUE, select = eem_list %>% eem_names() %>% .[c(1:16,205,208)], spp = 9, cores = cores)

## ----eval=FALSE, include=TRUE--------------------------------------------
#  ctol <- 10^-8 # tolerance in PARAFAC analysis
#  
#  pf4 <- eem_parafac(eem_list_ex, comps = seq(dim_min,dim_max), normalise = TRUE, const = c("nonneg", "nonneg", "nonneg"), maxit = maxit, nstart = nstart, ctol = ctol, cores = cores)
#  
#  pf4 <- lapply(pf4, eempf_rescaleBC, newscale = "Fmax")

## ----eval=TRUE, include=TRUE, fig.width=7, fig.height=6------------------
eempf_compare(pf4)

eempf_leverage_plot(eempf_leverage(pf4[[4]]))

eempf_corplot(pf4[[4]], progress = FALSE)

## ----eval=TRUE, include=TRUE, fig.width=7, fig.height=6------------------
eempf_comp_load_plot(pf4[[4]])
#eempf_plot_comps(pf4[4], type = 2)

## ----eval=TRUE, include=TRUE, fig.width=7, fig.height=8------------------
# plot components in each sample, residual and whole sample
eempf_residuals_plot(pf4[[4]], eem_list, select = eem_names(eem_list)[10:14], cores = cores)


## ----eval=FALSE, include=TRUE--------------------------------------------
#  #calculate split_half analysis
#  sh <- splithalf(eem_list_ex, 6, normalise = TRUE, rand = FALSE, cores = cores, nstart = nstart, maxit = maxit, ctol = ctol)

## ----eval=TRUE, include=TRUE, fig.width=7--------------------------------
data(sh)

## ----eval=TRUE, include=TRUE, fig.width=7--------------------------------
splithalf_plot(sh)

## ----eval=FALSE, include=TRUE--------------------------------------------
#  sh_r <- splithalf(eem_list_ex, 6, normalise = TRUE, rand = TRUE, cores = cores, nstart = nstart, maxit = maxit, ctol = ctol)

## ----eval=TRUE, include=TRUE, fig.width=7--------------------------------
splithalf_plot(sh_r)

## ----eval=TRUE, include=TRUE, fig.width=7--------------------------------
tcc_sh_table <- splithalf_tcc(sh_r)

tcc_sh_table

## ----eval=TRUE, include=TRUE, fig.width=7--------------------------------
corcondia <- eempf_corcondia(pf4[[4]], eem_list_ex)

corcondia

## ----eval=TRUE, include=TRUE, fig.width=7--------------------------------
eemqual <- eempf_eemqual(pf4[[4]], eem_list_ex, sh_r, cores = cores)

eemqual

## ----eval=FALSE, include = TRUE------------------------------------------
#  eempf_openfluor(pf4[[4]], file = "my_model_openfluor.txt")

## ----eval=FALSE, include = TRUE------------------------------------------
#  eempf_report(pf4[[4]], export = "my_model_openfluor.txt", eem_list = eem_list, shmodel = sh, performanec = TRUE)

## ---- message=FALSE, warning=FALSE, include=FALSE------------------------
write.bibtex(file="references2.bib")


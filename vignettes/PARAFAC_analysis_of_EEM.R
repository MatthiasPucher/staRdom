## ---- message=FALSE, warning=FALSE, include=FALSE------------------------
library(knitcitations)
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
#  data(eem_list)

## ----eval=TRUE, include=TRUE---------------------------------------------
folder <- system.file("extdata/cary/scans_day_1", package = "eemR")
eem_list <- eem_read(folder)

## ----eval=TRUE, fig.width=7, message=FALSE, warning=FALSE, include=TRUE, paged.print=TRUE----
ggeem(eem_list)

## ----eval=TRUE, include=TRUE---------------------------------------------
data("absorbance")

## ----eval=TRUE, include=TRUE---------------------------------------------
meta <- read.table(system.file("extdata/metatable_eemR.csv",package = "staRdom"), header = TRUE, sep = " ", dec = ".", row.names = 1)

## ----eval=TRUE, include=TRUE---------------------------------------------
eem_list <- eem_name_replace(eem_list,c("\\(FD3\\)"),c(""))

## ----eval=TRUE, include=TRUE---------------------------------------------
eem_list <- eem_list <- eem_remove_blank(eem_list)

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
remove_scatter <- c("raman1" = TRUE, "raman2" = TRUE, "rayleigh1" = TRUE, "rayleigh2" = TRUE)
remove_scatter_width <- c(15,15,15,15)

eem_list <- eem_rem_scat(eem_list, remove_scatter = remove_scatter, remove_scatter_width = remove_scatter_width)

## ----eval=TRUE, fig.width=7, message=FALSE, warning=FALSE, include=TRUE, paged.print=TRUE----
ggeem(eem_list)

## ----eval=TRUE, include=TRUE---------------------------------------------
eem_list <- eem_interp(eem_list, cores = cores)

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
data(eem_list)

## ----include=TRUE,eval=FALSE---------------------------------------------
#  eem_list <- eem_import_dir(dir)

## ----include=TRUE--------------------------------------------------------
# minimum and maximum of numbers of components
dim_min <- 5
dim_max <- 8

## ----eval=FALSE,include=TRUE---------------------------------------------
#  nstart <- 10 # number of similar models from which best is chosen
#  cores <- parallel::detectCores()/2 # use all cores but do not use all threads
#  maxit = 500 # maximum number of iterations in PARAFAC analysis
#  ctol <- 10^-5 # tolerance in PARAFAC analysis
#  
#  pfres_comps <- eem_parafac(eem_list, comps = seq(dim_min,dim_max), normalise = TRUE, maxit = maxit, nstart = nstart, ctol = ctol, cores = cores)

## ----include=TRUE--------------------------------------------------------
data(pfres_comps1)

## ----eval=TRUE, include=TRUE, fig.width=7, fig.height=6------------------
eempf_compare(pfres_comps)

## ----eval=TRUE, include=TRUE---------------------------------------------
comps <- 6

cp_out <- pfres_comps[[which(comps==seq(dim_min,dim_max))]]

## ----fig.width=7---------------------------------------------------------
# calculate leverage
cpl <- eempf_leverage(cp_out)
# plot leverage (nice plot)
eempf_leverage_plot(cpl,qlabel=0.1)
# plot leverage, not so nice plot but interactive to select what to exclude
# saved in exclude, can be used to start over again with eem_list_ex <- eem_list %>% eem_exclude(exclude) above
exclude <- eempf_leverage_ident(cpl,qlabel=0.1)

## ----eval=TRUE, include=TRUE---------------------------------------------
# samples, excitation and emission wavelengths to exclude, makes sense after calculation of leverage
exclude <- list("ex" = c(200,205,210,215,220,225,230,235,240,245,250,255,260,265,270,275,280,285,290,295,300),
                "em" = c(534,536,538,540,542,544,546,548,550,552,554,556,558,560,562,564,566,568,570,572,574,576,578,580,582,584,586,588,590,592,594,596,598,600),
                "sample" = c("sample87","sample78","sample95","sample12","sample17","sample51")
)

# exclude outliers if neccessary. if so, restart analysis
eem_list_ex <- eem_exclude(eem_list, exclude)

## ----eval=FALSE, include=TRUE--------------------------------------------
#  pfres_comps2 <- eem_parafac(eem_list_ex, comps = seq(dim_min,dim_max), normalise = TRUE, maxit = maxit, nstart = nstart, ctol = ctol, cores = cores)

## ----include=TRUE--------------------------------------------------------
data(pfres_comps2)

## ----eval=TRUE, include=TRUE, fig.width=7, fig.height=6------------------
eempf_compare(pfres_comps2)

## ----eval=TRUE, include=TRUE, fig.width=7, fig.height=6------------------
comps <- 6

cp_out <- pfres_comps2[[which(comps==seq(dim_min,dim_max))]]

## ----eval=TRUE, include=TRUE, fig.width=7, fig.height=6------------------
eempf_comp_load_plot(cp_out)

## ----eval=FALSE, include=TRUE, fig.width=7-------------------------------
#  eempf_comps3D(cp_out)

## ----eval=TRUE, include=TRUE, fig.width=7--------------------------------
# check for correlation between components table
# high correlations should be avoided
# try to normalise data or remove outliers as first step
eempf_cortable(cp_out)
# plot correlations
eempf_corplot(cp_out)

## ----eval=TRUE, include=TRUE, fig.width=7, fig.height=8------------------
# plot components in each sample, residual and whole sample
eempf_residuals_plot(cp_out, eem_list, select = eem_names(eem_list)[10:14], cores = cores)


## ----eval=TRUE, include=TRUE, fig.width=7, fig.height=8------------------
# plot components in each sample, residual and whole sample
eempf_residuals_plot(cp_out, eem_list, select = eem_names(eem_list)[c(10,11,13:16)], residuals_only = TRUE, cores = cores, spp = 6)


## ----eval=TRUE, include=TRUE, fig.width=7, fig.height=8------------------
# plot components in each sample, residual and whole sample
eempf_residuals_plot(cp_out, eem_list, select = c("sample12","sample17"), residuals_only = TRUE, cores = cores, spp = 6)


## ----eval=FALSE, include=TRUE--------------------------------------------
#  #calculate split_half analysis
#  sh <- splithalf(eem_list_ex, comps, normalise = TRUE, rand = FALSE, cores = cores)

## ----eval=TRUE, include=TRUE, fig.width=7--------------------------------
data(sh)

## ----eval=TRUE, include=TRUE, fig.width=7--------------------------------
splithalf_plot(sh)

## ----eval=TRUE, include=TRUE, fig.width=7--------------------------------
tcc_sh_table <- splithalf_tcc(sh)

tcc_sh_table

## ---- message=FALSE, warning=FALSE, include=FALSE------------------------
write.bibtex(file="references2.bib")


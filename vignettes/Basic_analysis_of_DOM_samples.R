## ----message=FALSE, warning=FALSE, include=FALSE------------------------------
library(dplyr)
library(tidyr)

## ----eval=FALSE, message=FALSE, warning=FALSE, include=TRUE, paged.print=FALSE----
#  title: "DOM analysis"
#  subtitle: "EEM peak picking, absorbance slope parameters"
#  date: "`r format(Sys.time(), '%B %e %Y')`"
#  author: "WCL"
#  knit: (function(inputFile, encoding) {
#        output_dir = 'C:/some_folder/output/';
#        rmarkdown::render(inputFile,
#                          encoding=encoding,
#                          output_file=file.path(output_dir,paste0('DOM_EEM_analysis_report_',format(Sys.time(), "%Y%m%d_%H%M%S"),'.htm')))
#                          })

## ----eval=FALSE---------------------------------------------------------------
#  # Set the directory where all output files are put in.
#  # The directory is automatically created if it does not exist.
#  # Folder delimiters can be / or \\. \, as it is usually used in Windows, will not work!
#  output_dir = "C:/some_folder/another_folder" # e.g. output_dir = "C:/some_folder/output/"

## ----eval=FALSE---------------------------------------------------------------
#  # Set the directory with your sample files. Please see eem_read() help for details on file formats.
#  # Sub folders are read in and are considered different sample sets.
#  # Import is done with eem_read() (package eemR), please see details there.
#  # The template refers to data coming with the package. Please use your data
#  # by setting the path to your files!
#  sample_dir = "C:/some_folder/input/fluor/" # e.g. sample_dir = "C:/some_folder/input/fluor/", system.file() accesses the example data coming with the package!
#  # Set the used instrument (with hyphens!):
#  # Cary Eclipse: "cary"
#  # Aqualog: "aqualog"
#  # Shimadzu: "shimadzu"
#  # Fluoromax-4: "fluoromax4"
#  # And furthermore, without hyphens:
#  # generic csv, excitation column-wise: eem_csv
#  # generic csv, emission column-wise: eem_csv2
#  # Hitachi F-7000: eem_hitachi
#  fluorometer = eem_csv

## ----eval=FALSE---------------------------------------------------------------
#  ### Absorbance data ###
#  #~~~~~~~~~~~~~~~~~~~~~#
#  # Absorbance data is read from *.TXT or *.CSV files.
#  # Either a directory containing one or more files can be named or a single file containing all samples.
#  # Absorbance data is used for inner-filter-effect correction and calculation of the slope parameters.
#  # Those steps can be skipped but keep in mind it is important for a profound analysis!
#  #
#  # path of adsorbance data as directory or single file, sub folders are not read:
#  absorbance_dir = "C:/some_folder/input/absorbance/" # e.g. absorbance_dir = "C:/some_folder/input/absorbance/", system.file() accesses the exmaple data coming with the package!
#  
#  # Path length of absorbance measurement in cm that was used in absorbance measurement.
#  # If it is set to "meta" data from the metadata table is used (details see below).
#  absorbance_path = 5 # e.g. absorbance_path = 5

## ----eval=FALSE---------------------------------------------------------------
#  ### Meta data ###
#  #~~~~~~~~~~~~~~~#
#  # Adding a table with meta data is OPTIONAL!
#  # The table can contain dilution factors, path lengths of
#  # the photometer and raman areas and is intended
#  # for cases where different values should be used for different
#  # samples. Each column can be used optionally.
#  
#  # read table with metadata as *.TXT or *.CSV
#  # either a path or FALSE if no metadata file is used.
#  metadata = system.file("extdata/metatable.csv",package = "staRdom") # e.g. metadata = "C:/some_folder/input/metatable.csv"", system.file() accesses the exmaple data coming with the package!
#  
#  ### Meta data: names of columns ###
#  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#  # designation of column with sample names
#  col_samples = "sample"
#  
#  # if you want to use dilution factors (e.g. 10 if 1 part sample and 9 parts solvent) from the meta data table, state the name
#  # of the column containing the dilution data and set dilution = "meta" (below)
#  col_dilution = "dilution"
#  
#  # if you want to use the cuvette length (in cm) for the absorbance from the meta data table,
#  # state the name of the column containing the cuvette lengths and set absorbance_path = "meta" (below)
#  col_cuv_len = "cuv_len"
#  
#  # if you want to use the raman area (under the curve) data from the meta data table, state the name
#  # of the column containing the raman areas and set raman_normalisation = "meta" (below)
#  col_raman_area = "raman"

## ----eval=FALSE---------------------------------------------------------------
#  #### Spectral correction of EEMs ####
#  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#  # Some instruments, but not all need a spectral correction to compensate
#  # for specific deviations in the measurements. A vector for emission and
#  # excitation is used each. EEMs are cut to match the wavelength range of
#  # the vectors if used. Please provide paths to csv tables containing
#  # wavelengths in the first column and correction values in the second. If
#  # you do not want spectral correction to be done, setting these two input
#  # files it not necessary.
#  # Emission correction vector
#  Emcor <- system.file("extdata/CorrectionFiles/mcorrs_4nm.csv",package="staRdom") # e.g. "C:\folder\emcor.csv", FALSE
#  # Excitation correction vector
#  Excor <- system.file("extdata/CorrectionFiles/xc06se06n.csv",package="staRdom")

## ----eval=FALSE---------------------------------------------------------------
#  ### Table output ###
#  #~~~~~~~~~~~~~~~~~~#
#  # Write a table with peaks and slope parameters.
#  # Written as xls or, in case of missing java environment or the package xlsx as csv.
#  output_xls = TRUE # e.g. TRUE
#  
#  # In case of a csv export you can define the separator and the decimal point here.
#  out_sep_dec = c("\t",".") # e.g. out_sep_dec = c("\t",".")

## ----eval=FALSE---------------------------------------------------------------
#  ### Plot settings PNG ###
#  #~~~~~~~~~~~~~~~~~~~~~~~#
#  # State whether you want pngs of the single EEM spectra written in your output directory
#  output_single_png = FALSE # e.g. TRUE
#  
#  ## State whether you want pngs of multiple EEM spectra written in your output directory
#  output_overview_png = FALSE # e.g. TRUE
#  
#  ## number of EEM spectra plottet in each overview image
#  overview_number = 6 # e.g. 6
#  
#  # The scaling of the different sample plots can be chosen.
#  # Either all samples are coloured according to the range of the
#  # complete sample set (TRUE) or each plot is scaled separately (FALSE).
#  scale_col = FALSE # e.g. TRUE
#  
#  ### Plot settings report ###
#  #~~~~~~~~~~~~~~~~~~~~~~~~~~#
#  # This block defines which plots are included in the report.
#  #
#  # Add plots with several EEM samples per plot.
#  # The number per plot is defined by overview_number above.
#  overview = TRUE # e.g. TRUE
#  
#  # State whether you want plots from single EEM spectra in the report.
#  single_plots = FALSE # e.g. TRUE

## ----eval=FALSE---------------------------------------------------------------
#  #### Save data for further analysis in R ####
#  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#  # File name where data is stored in RData format in the output directory.
#  # Set to FALSE if you dont want your eem data saved.
#  # Date, time and file extension is added automatically so you do not overwrite previous saved data.
#  data_file = "eem_data" # e.g. "eem_data"" or FALSE
#  
#  # Desired name for the variable containing the eem data.
#  eem_name = "eem_list"

## ----eval=FALSE---------------------------------------------------------------
#  #### Normalising absorbance data to baseline ####
#  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#  # Absorbance data can be corrected by subtracting a baseline value from each sample.
#  # In high wavelength ranges (default 680-700 nm), the absorbance is assumed to be 0.
#  # The average value of that range (or any other range) is subtracted from the whole spectrum.
#  # abs_norm can be set TRUE to use the default range, you can specify the desired range by a vector of length 2 and you can set it FALSE to skip this correction.
#  abs_norm = TRUE # e.g. TRUE, c(700,800)

## ----eval=FALSE---------------------------------------------------------------
#  ### Correction of diluted samples ###
#  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#  # Set a dilution factor if your samples were diluted.
#  # All samples are multiplied with this factor.
#  # Please use a meta table (above) if your dilutions are differing
#  # 1 for no dilution, 10 for dilution 1:10 (1 part sample and 9
#  # parts ultrapure water), "meta" for data from meta table
#  dilution = "meta" # e.g. 1 for undiluted samples

## ----eval=FALSE---------------------------------------------------------------
#  # In case of diluted samples, two absorbance measurements of the
#  # same sample in different dilutions might be present. If this is
#  # the case, EEMs are renamed to the undiluted sample, absorbance
#  # data might be multiplied by the dilution factor if it is only
#  # presentas diluted sample. This can be done automatically. In the
#  # final protocol a table shows, what has been done to the samples.
#  # Please check this table and see, if the output is what you
#  # wanted it to be!
#  dil_sample_name_correction = FALSE

## ----eval=FALSE---------------------------------------------------------------
#  #### Spectral correction of EEMs ####
#  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#  # Some instruments, but not all need a spectral correction to compensate
#  # for specific deviations in the measurements. Please sepecify, if you want
#  # spectral correection to be done.
#  spectral_cor = TRUE # e.g. TRUE, set to FALSE, if your instrument already provided EEMs with spectral correction

## ----eval=FALSE---------------------------------------------------------------
#  ### Cut data to certain range ###
#  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#  # Set a vector with range of wavelengths to be plotted and saved.
#  # Peak picking is done before range reduction.
#  # Emission wavelength:
#  em_range = c(0,Inf) # e.g. c(300,500), c(0,Inf) to use everything
#  
#  # Excitation wavelength:
#  ex_range = c(0,Inf) # e.g. c(300,500), c(0,Inf) to use everything
#  
#  # Cut all samples to fit largest range available in all samples
#  cut_range_to_smallest = FALSE # e.g. FALSE

## ----eval=FALSE---------------------------------------------------------------
#  ### Blank correction ###
#  #~~~~~~~~~~~~~~~~~~~~~~#
#  # A blank sample is subtracted from each sample. Blank samples have to be
#  # in the same (sub)folder as the according EEM samples. So different blanks are used
#  # for different subsets. The file names of the blanks have to contain nano,
#  # miliq, milliq, mq or blank (cases are ignored). Other samples must not
#  # contain these words in their names respectively!
#  blank_correction = FALSE # e.g. FALSE

## ----eval=FALSE---------------------------------------------------------------
#  ### Inner filter effect correction ###
#  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#  # Inner filter effects are corrected. Absorbance data is needed. File or column designations
#  # of the absorbance data have to resamble file names of the EEM data.
#  ife_correction = TRUE # e.g. FALSE

## ----eval=FALSE---------------------------------------------------------------
#  ### Remove scattering and interpolate missing data ###
#  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#  # Scattering is removed from the EEM spectra.
#  remove_scatter <- c(TRUE, TRUE, TRUE, TRUE) # logical values, ordered by raman1,raman2,rayleigh1,rayleigh2
#  
#  # Set the width of removed scatter slot (usually 10 to 20).
#  # If you can still see traces of scattering after interpolation,
#  # this value should be increased. You can specify a vector containing
#  # separate widths for each scatter c(15,16,16,14), ordered by raman1,raman2,rayleigh1,rayleigh2
#  # In case one or more scatter peaks are skipped, this vector must remain of length 4 and positions of the certain widths must be kept.
#  remove_scatter_width = c(15,15,15,15) # e.g. 15 or c(15,15,15,15)
#  
#  # state whether removed scattering should be interpolated
#  interpolation <- TRUE # e.g. TRUE

## ----eval=FALSE---------------------------------------------------------------
#  ### Raman normailsation ###
#  #~~~~~~~~~~~~~~~~~~~~~~~~~#
#  # State whether a Raman normalisation should be performed
#  # Either "blank" if a blank is present in each (sub)folder of the EEM data.
#  # Blank samples have to be in the same (sub)folder as the EEM samples. So
#  # different blanks are used for different subsets. The file names of the
#  # blanks have to contain nano, miliq, milliq, mq or blank (cases are ignored).
#  # Other samples must not contain these words in their names respectively!
#  # Normalisation is then calculated with this blank, the raman area as a number
#  # or "meta" if the raman areas should be taken from the meta data table.
#  raman_normalisation = "blank" # e.g. "blank", FALSE, 160, "meta"

## ----eval=FALSE---------------------------------------------------------------
#  ### Smooth data for peak picking ###
#  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#  # Moving window size for smoothing data along excitation wavelengths.
#  # Data must be interpolated if you want to use smoothing.
#  # This is used for peak picking but not saved.
#  smooth = 4 # e.g. FALSE, 4

## ----eval=FALSE---------------------------------------------------------------
#  #############################################
#  #                                           #
#  #       THERE ARE NO SETTINGS BELOW.        #
#  #  YOU CAN KLICK "KNIT" AT THE MENU BAR.    #
#  #  In case of errors, chunk-wise execution  #
#  #     of the code can reveal problems!      #
#  #                                           #
#  #     Please read the help of the used      #
#  #        functions if you encounter         #
#  #               any problems:               #
#  #   Press F1 while cursor in function or    #
#  #   type help(function) in command line!    #
#  #                                           #
#  #             Please read the               #
#  #        error messages carefully!          #
#  #    Naming of the input files and table    #
#  #     column and row names is crucial!      #
#  #                                           #
#  #############################################


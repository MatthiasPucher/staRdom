## ----message=FALSE, warning=FALSE, include=FALSE-------------------------
library(knitcitations)
cleanbib()
options("citation_format" = "pandoc")
bibliography() #style="apalike"

## ----eval=FALSE, message=FALSE, warning=FALSE, include=FALSE, paged.print=FALSE----
#  title: "DOM analysis"
#  subtitle: "EEM peak picking, absorbance slope parameters"
#  date: "`r format(Sys.time(), '%B %e %Y')`"
#  author: "WCL"
#  knit: (function(inputFile, encoding) {
#        output_dir = 'C:/some_folder/another_folder';
#        rmarkdown::render(inputFile,
#                          encoding=encoding,
#                          output_file=file.path(output_dir,paste0('DOM_EEM_analysis_report_',format(Sys.time(), "%Y%m%d_%H%M%S"),'.htm')))
#                          })

## ----eval=FALSE----------------------------------------------------------
#  # Set the directory where all output files are put in.
#  # The directory is automatically created if it does not exist.
#  output_dir = "C:/some_folder/another_folder" # e.g. output_dir = "C:/some_folder/output/"

## ----eval=FALSE----------------------------------------------------------
#  ### Directory containing EEM data ###
#  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#  # Set the directory with your sample files. Please see eem_read() help for details on file formats.
#  # Sub folders are read in and are considered different sample sets.
#  # Import is done with eem.read() (package eemR), please see details there.
#  # The template refers to data coming with the package. Please use your data
#  # by setting the path to your files!
#  sample_dir = "C:/some_folder/input/fluor/" # e.g. sample_dir = "C:/some_folder/input/fluor/"

## ----eval=FALSE----------------------------------------------------------
#  ### Absorbance data ###
#  #~~~~~~~~~~~~~~~~~~~~~#
#  # Absorbance data is read from *.TXT or *.CSV files.
#  # Either a directory containing one or more files can be named or a single file containing all samples.
#  # Absorbance data is used for inner-filter-effect correction and calculation of the slope parameters.
#  # Those steps can be skipped but keep in mind it is important for a profound analysis!
#  #
#  # path of adsorbance data as directory or single file, sub folders are not read:
#  absorbance_dir = "C:/some_folder/input/absorbance/" # e.g. absorbance_dir = "C:/some_folder/input/absorbance/"
#  
#  # Cuvette length in cm that was used in absorbance measurement.
#  # If it is set to "meta" data from the metadata table is used (details see below).
#  absorbance_cuv_len = 5 # e.g. absorbance_cuv_len = 5

## ----eval=FALSE----------------------------------------------------------
#  ### Meta data ###
#  #~~~~~~~~~~~~~~~#
#  # Adding a table with meta data is OPTIONAL!
#  # The table can contain dilution factors, cuvette lengths of
#  # the photometer and absorbance cuvette lengths and is intended
#  # for cases where different values should be used for different
#  # samples. Each column can be used optionally.
#  
#  # read table with metadata as *.TXT or *.CSV
#  # either a path or FALSE if no metadata file is used.
#  metadata = system.file("extdata/metatable.csv",package = "staRdom") # e.g. metadata = "C:/some_folder/input/metatable.csv"
#  
#  ### Meta data: names of columns ###
#  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#  # column with sample names
#  col_samples = "sample"
#  
#  # if you want to use dilution factors (e.g. 10 if 1 part sample and 9 parts solvent) from the meta data table, state the name
#  # of the column containing the dilution data and set dilution = "meta" (below)
#  col_dilution = "dilution"
#  
#  # if you want to use the cuvette length (in cm) for the absorbance from the meta data table,
#  # state the name of the column containing the cuvette lengths and set absorbance_cuv_len = "meta" (below)
#  col_cuv_len = "cuv_len"
#  
#  # if you want to use the raman area (under the curve) data from the meta data table, state the name
#  # of the column containing the raman areas and set raman_normalisation = "meta" (below)
#  col_raman_area = "raman"

## ----eval=FALSE----------------------------------------------------------
#  ### Table output ###
#  #~~~~~~~~~~~~~~~~~~#
#  # Write a table with peaks and slope parameters.
#  # Written as xls or, in case of missing java environment as csv.
#  output_xls = TRUE # e.g. TRUE
#  
#  # In case of csv export you can define the separator and the decimal point here.
#  out_sep_dec = c("\t",".") # e.g. out_sep_dec = c("\t",".")

## ----eval=FALSE----------------------------------------------------------
#  ### Plot settings PNG ###
#  #~~~~~~~~~~~~~~~~~~~~~~~#
#  # The scaling of the different sample plots can be chosen.
#  # Either all samples are coloured according to the range of the
#  # complete sample set (TRUE) or each plot is scaled separately (FALSE).
#  scale_col = FALSE # e.g. TRUE
#  
#  # State whether you want pngs of the single EEM spectra written in your output directory
#  output_single_png = FALSE # e.g. TRUE
#  
#  ## State whether you want pngs of multiple EEM spectra written in your output directory
#  output_overview_png = FALSE # e.g. TRUE
#  
#  ## number of EEM spectra plottet in each overview image
#  overview_number = 6 # e.g. 6
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

## ----eval=FALSE----------------------------------------------------------
#  ### Correction of diluted samples ###
#  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#  # Set a dilution factor if your sample was diluted.
#  # All samples are multiplied with this factor.
#  # Please use a meta table (above) if your dilutions are differing
#  # 1 for no dilution, 10 for dilution 1:10 (1 part sample and 9
#  # parts milliq), "meta" for data from meta table
#  dilution = "meta" # e.g. 1

## ----eval=FALSE----------------------------------------------------------
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

## ----eval=FALSE----------------------------------------------------------
#  ### Blank correction ###
#  #~~~~~~~~~~~~~~~~~~~~~~#
#  # A blank sample is substracted from each sample. Blank samples have to be
#  # in the same (sub)folder as the EEM samples. So different blanks are used
#  # for different subsets. The file names of the blanks have to contain nano,
#  # miliq, milliq, mq or blank (cases are ignored). Other samples must not
#  # contain these words in their names respectively!
#  blank_correction = FALSE # e.g. FALSE

## ----eval=FALSE----------------------------------------------------------
#  ### Inner filter effect correction ###
#  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#  # Inner filter effects are corrected. Absorbance data is needed. File or column names
#  # of the absorbance data have to resamble file names of the EEM data.
#  ife_correction = TRUE # e.g. FALSE

## ----eval=FALSE----------------------------------------------------------
#  ### Remove scattering and interpolate missing data ###
#  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#  # Scattering is removed from the EEM spectra.
#  remove_scatter <- c() # Please leave that as it is ;-)
#  
#  # remove Raman scatter order 1
#  remove_scatter["raman1"] = TRUE # e.g. TRUE
#  
#  # remove Raman scatter order 2
#  remove_scatter["raman2"] = TRUE # e.g. TRUE
#  
#  # remove Rayleigh scatter order 1
#  remove_scatter["rayleigh1"] = TRUE # e.g. TRUE
#  
#  # remove Rayleigh scatter order 2
#  remove_scatter["rayleigh2"] = TRUE # e.g. TRUE
#  
#  # Set the width of removed scatter slot (usually 10 to 16).
#  # If you can still see traces of scattering after interpolation,
#  # this value should be increased. You can specify a vector containing
#  # separate widths for each scatter c(15,16,16,14), ordered by raman1,raman2,rayleigh1,rayleigh2.
#  remove_scatter_width = c(15,15,15,15) # e.g. 15 or c(15,15,15,15)
#  
#  # state whether removed scattering should be interpolated
#  interpolation <- TRUE # e.g. TRUE

## ----eval=FALSE----------------------------------------------------------
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

## ----eval=FALSE----------------------------------------------------------
#  ### Smooth data for peak picking ###
#  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#  # Moving window size for smoothing data along excitation wavelengths.
#  # Data must be interpolated if you want to use smoothing.
#  # This is used for peak picking but not saved.
#  smooth = 4 # e.g. FALSE, 4

## ----eval=FALSE----------------------------------------------------------
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
#  #    Naming of the imput files and table    #
#  #     column and row names is crucial!      #
#  #                                           #
#  #############################################

## ----message=FALSE, warning=FALSE, include=FALSE-------------------------
write.bibtex(file="references.bib")


##############################
### Allium metabolomics script 
### by Jakob Franke, Leibniz University Hannover, Germany
### jakob.franke@botanik.uni-hannover.de
### last update: 19 January 2023
##############################
### parts of this script and the general logic for filtering XCMS results for isotopologue search
### are taken or adapted from Nett et al., AIChe Journal 2018, 64, 4319-4330
### see also https://github.com/Stanford-ChEMH-MCAC/d2o_metabolomics






#################################
### Settings and preparations ###
#################################

# load required libraries
library(xcms)
library(CAMERA)
library(tidyverse)
library(ecipex)
library(openxlsx)

# set working directory
setwd("D:/Seafile/LUH/projects/alliin/230114 allium metabolomics scripts for publication")


# LCMS file settings
file_list <- "sample_descriptions.csv" # csv file with names of LCMS files and their sample group (see more information below)
file_path <- "F:/OneDrive/LUH/projects/alliin/221117 schnittlauch orbitrap labelling all replicates/schnittlauch_high_res_full_ms/mzxml" # location where LCMS files are stored
file_extension <- "*.mzXML" # file type of input LCMS data


# program settings; these can be used to skip certain program parts for speeding up the analysis
run_XCMS <- TRUE
run_isotope_analysis <- TRUE
output_prefix <- "results_230119/allium_metabolomics_" # subfolder and file name prefix for output files


# chromatogram filter settings
min_RT <- 60 # minimum retention time in seconds
max_RT <- 600 # maximum retention time in seconds
RTwindow <- 5 # retention time tolerance in seconds when comparing two features 
min_mz <- 60 # minimum m/z value
max_mz <- 900 # maximum m/z value
min_intensity <- 1e3 # minimum peak intensity to include in analysis


# labelling search filter settings: these settings define which criteria mass feature pairs need to fulfil
# example for 13C search: for min_labelled_Cs = 1 and max_labelled_Cs = 3 the script would search for co-eluting mass features differing in mass by 1 * 1.0034, 2 * 1.0034 or 3 * 1.0034
# example for sulfur search: for min_labelled_S = 1 and max_labelled_S the script would search for co-eluting mass features that differ in mass by 1 * 1.9958 or 2 * 1.9958, and also fall within the isotope_ratio_tolerance threshold
min_labelled_Cs <- 1 # minimum number of carbons that have to be labelled
max_labelled_Cs <- 5 # maximum number of carbons that can be labelled
min_labelled_S <- 1 # minimum number of sulfur atoms
max_labelled_S <- 2 # maximum number of sulfur atoms
diff_ppm <- 2.5 # ppm tolerance when searching for specific mass differences
min_foldchange <- 10 # minimum fold change which the heavy feature must show in labelled vs. unlabelled samples
isotope_ratio_tolerance <- 0.02 # for sulfur search: absolute value by which the 34S/32S ratio can differ from the expected 4.3%/95.0% ( = 4.5%) ratio; example: a value of 0.02 leads to acceptance of a 34S/32S ratio of 0.045 +/- 0.02


# xcms settings
xcms_ppm <- 2.5 # ppm threshold for XCMS CentWave
xcms_mzdiff <- 0.001 # mzdiff parameter for XCMS CentWave
xcms_prefilter <- c(3,1000) # prefilter for XCMS CentWave; at least 3 peaks with at least 1000 intensity required
xcms_peakwidth <- c(5,20) # peakwidth for XCMS CentWave; typical peak width 5 - 20 s
xcms_snthresh <- 10 # snthresh parameter for XCMS CentWave; minimum required signal/noise ratio
xcms_noise <- 1000 # noise parameter for XCMS CentWave; minimum ion intensity
xcms_expand_rt <- 4 # retention time window in s for merging neighbouring peaks during refineChromPeaks
xcms_align_bin_size <- 0.6 # align bin size for RT adjustment using Obiwarp
xcms_group_bw <- 5 # peak grouping bandwidth
xcms_group_bin_size <- 0.001 # peak grouping bin size
xcms_group_min_fraction <- 0.4 # peak grouping minimum fraction of samples that need to contain a peak


# obtain theoretical values for identification of carbon and sulfur isotopes
# exact mass differences of 13C vs. 12C and 34S vs. 32S come from ecipex package
mass_diff_C <- ecipex('C')[[1]]$mass %>% diff # mass difference between 13C and 12C
mass_diffs_C <- outer(c(min_labelled_Cs:max_labelled_Cs), mass_diff_C) %>% as.numeric %>% sort # generate multiples of mass difference based on min_ and max_labelled_Cs

mass_diff_S <- ecipex('S')[[1]]$mass[1:2] %>% diff # mass difference between 34S and 32S
mass_diffs_S <- outer(c(min_labelled_S:max_labelled_S), mass_diff_S) %>% as.numeric %>% sort # generate multiples of mass difference based on min_ and max_labelled_S

isotope_ratio_S <- ecipex('S')[[1]]$abundance[2] / ecipex('S')[[1]]$abundance[1] # obtain the natural abundance ratio of 34S / 32S


# load csv file that contains names of LCMS files to be processed and their sample group (labelled vs. unlabelled)
# required format:
# sample_id,label
file_info <- read.csv(file_list)

# collect information for loading LCMS files
file_info <-
  file_info %>%
  arrange(sample_id) %>%
  mutate(full_name = # add complete absolute path for all samples
           list.files(path=file_path, 
                      pattern=file_extension,
                      full.names = T)
  )


# check that samples were correctly found on the file system and collect file names in my_files
my_files <-
  file_info %>%
  filter(str_detect(full_name, sample_id)) %>%
  pull(full_name)





###################################
### XCMS analysis of LCMS files ###
###################################
if (run_XCMS) {
  # load LCMS data from disk and restrict to desired retention time and m/z window
  raw_data <- readMSData(files = my_files, pdata = new("NAnnotatedDataFrame", file_info), mode = "onDisk")
  raw_data <- filterRt(raw_data, c(min_RT, max_RT))
  raw_data <- filterMz(raw_data, c(min_mz, max_mz))
  
  # perform peak search with XCMS CentWave
  cwp <- CentWaveParam(ppm = xcms_ppm,
                       mzdiff = xcms_mzdiff,
                       prefilter = xcms_prefilter,
                       peakwidth = xcms_peakwidth,
                       snthresh = xcms_snthresh,
                       noise = xcms_noise)
  xdata <- findChromPeaks(raw_data, param = cwp)
  
  # refine XCMS data
  xdata <- refineChromPeaks(xdata, MergeNeighboringPeaksParam(expandRt = xcms_expand_rt, ppm = xcms_ppm))
  
  xdata <- adjustRtime(xdata, param = ObiwarpParam(binSize = xcms_align_bin_size))
  
  pdp <- PeakDensityParam(sampleGroups = xdata@phenoData@data$label,
                          minFraction = xcms_group_min_fraction,
                          bw = xcms_group_bw,
                          binSize = xcms_group_bin_size)
  xdata <- groupChromPeaks(xdata, param = pdp)
  
  xdata <- fillChromPeaks(xdata) # fill missing values
  
  # convert xdata object to older xset object for use with CAMERA package
  # this requires that sampnames and sampclass are set manually
  xset <- as(xdata, "xcmsSet")
  my_classes <- # obtain sample classes from input data
    file_info %>%
    filter(str_detect(full_name, sample_id)) %>%
    pull(label) %>%
    as.character
  sampnames(xset) <- file_info$sample_id # set sample names based on input data sample ids
  sampclass(xset) <- my_classes # set sample classes (e.g. labelled or unlabelled)
  
  # save XCMS mass features to disk so that this part could be skipped when re-running the script
  save(xdata, file = "xdata.rdat")
  save(xset, file = "xset.rdat")
} else # load XCMS mass features from disk when run_XCMS is false
  load("xset.rdat")







##########################################
### Search for 13C and sulfur isotopes ###
##########################################
if (run_isotope_analysis) {
  # perform statistical analyses for all detected mass features and generate a feature table (diffreport from XCMS package), and annotate this feature table with CAMERA for potential adducts etc.
  feature_table <- annotateDiffreport(xset, class1="labelled", class2="unlabelled", missing = 0, ppm = xcms_ppm, mzabs = xcms_mzdiff)
  
  
  
  ##################
  ### 13C search ###
  ##################
  
  # identify which columns from the feature table contain the peak integrals?
  data_cols <- 
    feature_table %>%
    names %>%
    {. %in% (file_info %>% pull(sample_id))} %>%
    which
  
  # expand the feature table to include sample_id as an additional column
  long_table <-
    feature_table %>%
    gather(sample_id, intensity, data_cols) %>%
    left_join(file_info, by=c('sample_id')) %>%
    select(name, fold, tstat, pvalue, mzmed, rtmed, sample_id, label, intensity, isotopes)
  
  # remove features that are too small (i.e. min_intensity not reached) and fall out of the retention time window (min_RT / max_RT)
  filtered_table <-
    long_table %>%
    group_by(name, fold, tstat, pvalue, mzmed, rtmed, isotopes) %>%
    mutate(max_intensity = max(intensity)) %>% # add an additional column with the maximum intensity of a peak across all samples
    arrange(name) %>% 
    filter(max_intensity > min_intensity, rtmed > min_RT, rtmed < max_RT) %>% # filter out features that do not fulfill criteria
    arrange(-fold)
  
  # add columns with the mean intensity of every peak in the two sample groups (labelled vs. unlabelled)
  summary_table <-
    filtered_table %>%
    group_by(name, fold, tstat, pvalue, mzmed, rtmed, max_intensity, label, isotopes) %>%
    summarize(mean_intensity = mean(intensity)) %>%
    spread(label, mean_intensity) # this adds the mean intensity of all samples in the labelled group and the mean intensity of all samples in the unlabelled group
  
  
  # remove features that do not fall in the required mz window (min_mz / max_mz)
  summary_table_filtered <-
    summary_table %>%
    ungroup %>%
    mutate(foo = 1) %>%
    filter(mzmed <= max_mz, mzmed >= min_mz) 
  
  # create a complete pairwise comparison of all features against all other features (n x n) and filter out all features that don't show a mass difference matching to 13C vs. 12C
  feature_pairs <- 
    summary_table_filtered %>%
    left_join(., ., by='foo') %>% # form all feature pairs
    filter(mzmed.x < mzmed.y,  # ".x" is light feature, ".y" is heavy
           abs(rtmed.x - rtmed.y) <= RTwindow) %>% # retention time difference must fall within the RT window
    filter(isotopes.x == "" | grepl("\\[M\\]\\+", isotopes.x)) %>% # filter out isotopologues and multiply charged ions detected by CAMERA to reduce overall redundancy
    rowwise %>%
    mutate(mz_diff = abs(mzmed.x - mzmed.y), # add a column with the mass difference of the feature pair
           ppm = min(abs(mz_diff - mass_diffs_C) / mzmed.x * 1e6), # calculate the ppm difference of the observed mass difference to the closest theoretical mass difference for 13C vs. 12C; multiple 13Cs are also considered
           labelled_carbons = round(mz_diff)) %>% # calculate the putative number of labelled carbon atoms based on the observed mass difference
    filter(ppm <= diff_ppm) %>% # the calculated ppm difference for 13C vs. 12C must be below the threshold defined by diff_ppm
    arrange(-fold.y)
  
  # polish filtered feature pairs for output
  labelling_hits <- feature_pairs %>%
    ungroup %>% # required for calculating a scorer value
    filter(fold.y > min_foldchange, mzmed.y >= min_mz, mzmed.y <= max_mz) %>% # remove feature pairs where the heavy feature does not show a sufficient fold change (min_foldchange), and make sure that feature falls within mz range
    mutate(scorer = sqrt(rank(fold.y) * rank(max_intensity.x)), # calculate a scorer value based on the fold change of the heavy feature (ideally infinity - present in labelled samples, completely absent in unlabelled samples) and the maximum ion intensity
           intensity_ratio = labelled.y / labelled.x, # calculate the ratio of peak intensities of the heavy feature vs. the light feature in labelled samples; this represents the incorporation of label
           rtmed_min.x = rtmed.x / 60, # convert retention times from s to min 
           rtmed_min.y = rtmed.y / 60) %>% # convert retention times from s to min
    arrange(-scorer) %>% # sort by scorer value
    select(name.x, name.y, # only select columns required for output later
           mzmed.x, mzmed.y,
           rtmed_min.x, rtmed_min.y,
           fold.x, fold.y,
           mz_diff, ppm, labelled_carbons,
           intensity_ratio,
           max_intensity.x,
           scorer,
           isotopes.x,
           isotopes.y)
  # this completes the identification of features that show incorporation of 13C in labelled samples vs. unlabelled samples
  
  
  
  
  
  
  #####################
  ### Sulfur search ###
  #####################
  # next: search for feature pairs that indicate the presence of a sulfur based on mass difference (34S vs. 32S) and intensity ratio (4.3%/95.0%)
  # search is performed across all samples and results are pooled eventually (because sulfur metabolites can occur in unlabelled as well as labelled samples)
  for (i in 1:length(file_info$sample_id)) { # go through all LCMS files for the sulfur search
    # continue from the filtered table containing all features generated above and select only features in the current LCMS file and within the mz range (min_mz / max_mz)
    sulphur_table <- filtered_table %>%
      filter(sample_id == file_info$sample_id[i]) %>%
      ungroup %>%
      mutate(foo = 1) %>%
      filter(mzmed <= max_mz, mzmed >= min_mz)
    
    # create a pairwise comparison of all features (n x n) and filter out all features that do not show the expected mass difference 34S vs. 32S or a matching intensity ratio of 34S/32S
    sulphur_pairs_temp <- 
      sulphur_table %>%
      left_join(., ., by='foo') %>% # form all feature pairs
      filter(mzmed.x < mzmed.y,  # ".x" is light feature, ".y" is heavy
             abs(rtmed.x - rtmed.y) <= RTwindow) %>% # retention time difference must fall within the RT window
      filter(isotopes.x == "" | grepl("\\[M\\]\\+", isotopes.x)) %>% # filter out isotopologues and multiply charged ions detected by CAMERA to reduce overall redundancy
      rowwise %>%
      mutate(mz_diff = abs(mzmed.x - mzmed.y), # add a column with the mass difference of the feature pair
             ppm = min(abs(mz_diff - mass_diffs_S) / mzmed.x * 1e6), # calculate the ppm difference of the observed mass difference to the closest theoretical mass difference for 34S vs. 32S; multiple sulfur atoms are also considered
             isotope_ratio = intensity.y / intensity.x, # calculate the ratio of peak intensities of heavy feature vs. light feature
             isotope_multiple = round(isotope_ratio / isotope_ratio_S), # are there potentially multiple sulfur atoms, and how many?
             isotope_ratio_normalised = isotope_ratio / isotope_multiple) %>% # normalise the peak intensity ratio by the putative number of sulfur atoms
      rowwise %>%
      filter(ppm <= diff_ppm, # the calculated ppm difference for 34S vs. 32S must be below the threshold defined by diff_ppm
             isotope_ratio_normalised >= (isotope_ratio_S - isotope_ratio_tolerance) ^ (round(mz_diff) / 2), # the normalised isotope ratio must fall within the tolerance range (isotope_ratio_tolerance) of the theoretical isotope ratio (4.3%/95.0%)
             isotope_ratio_normalised <= (isotope_ratio_S + isotope_ratio_tolerance) ^ (round(mz_diff) / 2), # the normalised isotope ratio must fall within the tolerance range (isotope_ratio_tolerance) of the theoretical isotope ratio (4.3%/95.0%)
             isotope_multiple <= max_labelled_S) %>% # the maximum number of sulfur atoms must not exceed max_labelled_S
      arrange(-intensity.x) # sort by intensity
    
    
    # pool identified sulfur hits from each sample into the joined table sulphur_pairs
    if (i == 1) sulphur_pairs <- sulphur_pairs_temp
    else {
      sulphur_pairs <- full_join(sulphur_pairs, sulphur_pairs_temp)
    }
  }
  
  # after sulfur features were searched in all LCMS samples, remove all redundant hits that were found in multiple samples
  sulphur_pairs <- sulphur_pairs %>%
    distinct(name.x, name.y, .keep_all = TRUE)
  
  # get the names of all features that indicated 13C labelling
  labelling_features <- labelling_hits %>%
    pull(name.x, name.y) %>%
    unique
  
  # get the names of all features that indicated the presence of sulfur
  sulphur_features <- sulphur_pairs %>%
    pull(name.x) %>%
    unique
  
  # polish filtered sulfur pairs for output; add a column that indicates if a feature was also found during 13C labelling search
  sulphur_hits <- sulphur_pairs %>%
    ungroup %>% # required for calculating a scorer value
    mutate(scorer = rank(max_intensity.x), # calculate a scorer value for the sulfur hits based on the maximum peak intensity
           is_labelled = (name.x %in% labelling_features), # add a column that shows if the same feature was also found during the 13C labelling search
           rtmed_min.x = rtmed.x / 60, # convert retention times from s to min 
           rtmed_min.y = rtmed.y / 60, ) %>% # convert retention times from s to min 
    arrange(-scorer) %>% # sort by scorer value
    select(name.x, name.y, # only select columns required for output later
           mzmed.x, mzmed.y,
           rtmed_min.x, rtmed_min.y,
           fold.x, fold.y,
           mz_diff, ppm, isotope_ratio, isotope_multiple,
           max_intensity.x,
           scorer,
           is_labelled,
           isotopes.x,
           isotopes.y)
  
  # add columns to the labelling hits table that indicate if also a sulfur was found for the light and/or heavy feature
  labelling_hits_sulphur <- labelling_hits %>%
    mutate(sulphur_in_x = (name.x %in% sulphur_features),
           sulphur_in_y = (name.y %in% sulphur_features))
  # this completes the search for sulfur-containing mass features
  
  
  
  
  
  
  
  
  ##############
  ### Output ###
  ##############
  
  # here we generate some formatted xlsx tables that contain all found 13C labelling hits as well as all sulfur hits
  # also we generate a short text file that documents the settings and a quick summary how many hits were found and if the major expected peaks were correctly detected
  
  # formatting settings for Excel
  style_mass <- createStyle(numFmt = "0.0000") # mz values should have four decimal places
  style_RT <- createStyle(numFmt = "0.00") # retention times should have two decimal places
  style_general <- createStyle(numFmt = "0.00") # general numbers should have two decimal places
  style_intensity <- createStyle(numFmt = "0.00E+00") # peak intensities should use scientific notation
  style_intensity_ratio <- createStyle(numFmt = "0.0%") # peak intensity ratios should have one decimal place
  style_highlight <- createStyle(fontColour = "#006100", bgFill = "#C6EFCE") # cells should be highlighted with a light green background and dark green font
  style_no_highlight <- createStyle(fontColour = "#000000", bgFill = "#FFFFFF") # default format for non-highlighted cells (black font, white background)
  
  # generate xlsx output for 13C labelling hits
  wb <- createWorkbook()
  addWorksheet(wb, sheetName = "labelling_hits")
  writeData(wb, 1, labelling_hits_sulphur)
  addStyle(wb, 1, style = style_mass, rows = 1:10000, cols = 3:4, gridExpand = TRUE) # set desired style for the corresponding columns
  addStyle(wb, 1, style = style_mass, rows = 1:10000, cols = 9:10, gridExpand = TRUE) # set desired style for the corresponding columns
  addStyle(wb, 1, style = style_RT, rows = 1:10000, cols = 5:6, gridExpand = TRUE) # set desired style for the corresponding columns
  addStyle(wb, 1, style = style_general, rows = 1:10000, cols = 7:8, gridExpand = TRUE) # set desired style for the corresponding columns
  addStyle(wb, 1, style = style_general, rows = 1:10000, cols = 14, gridExpand = TRUE) # set desired style for the corresponding columns
  addStyle(wb, 1, style = style_intensity, rows = 1:10000, cols = 13, gridExpand = TRUE) # set desired style for the corresponding columns
  addStyle(wb, 1, style = style_intensity_ratio, rows = 1:10000, cols = 12, gridExpand = TRUE) # set desired style for the corresponding columns
  conditionalFormatting(wb, 1, rows = 1:10000, cols = 17:18, rule = "==TRUE", style_highlight) # set conditional formatting for highlighting cells
  conditionalFormatting(wb, 1, rows = 1:10000, cols = 17:18, rule = "!=TRUE", style_no_highlight) # set conditional formatting for highlighting cells
  setColWidths(wb, 1, cols = 1:16, widths = "auto") # adjust column width to the content
  saveWorkbook(wb, paste(output_prefix, "labelling.hits.xlsx", sep = ""), overwrite = TRUE) # save to file
  
  # generate xlsx output for sulfur hits
  wb <- createWorkbook()
  addWorksheet(wb, sheetName = "sulphur_hits")
  writeData(wb, 1, sulphur_hits)
  addStyle(wb, 1, style = style_mass, rows = 1:10000, cols = 3:4, gridExpand = TRUE) # set desired style for the corresponding columns
  addStyle(wb, 1, style = style_mass, rows = 1:10000, cols = 9:10, gridExpand = TRUE) # set desired style for the corresponding columns
  addStyle(wb, 1, style = style_RT, rows = 1:10000, cols = 5:6, gridExpand = TRUE) # set desired style for the corresponding columns
  addStyle(wb, 1, style = style_general, rows = 1:10000, cols = 7:8, gridExpand = TRUE) # set desired style for the corresponding columns
  addStyle(wb, 1, style = style_general, rows = 1:10000, cols = 14, gridExpand = TRUE) # set desired style for the corresponding columns
  addStyle(wb, 1, style = style_intensity, rows = 1:10000, cols = 13, gridExpand = TRUE) # set desired style for the corresponding columns
  addStyle(wb, 1, style = style_intensity_ratio, rows = 1:10000, cols = 11, gridExpand = TRUE) # set desired style for the corresponding columns
  conditionalFormatting(wb, 1, rows = 1:10000, cols = 15, rule = "==TRUE", style_highlight) # set conditional formatting for highlighting cells
  conditionalFormatting(wb, 1, rows = 1:10000, cols = 15, rule = "!=TRUE", style_no_highlight) # set conditional formatting for highlighting cells
  setColWidths(wb, 1, cols = 1:15, widths = "auto") # adjust column width to the content
  saveWorkbook(wb, paste(output_prefix, "sulphur_hits.xlsx", sep = ""), overwrite = TRUE) # save to file
  
  
  # generate a quick summary report text file
  sink(file = paste(output_prefix, "output.txt", sep = ""), type = "output")
  
  # add the date
  cat(output_prefix, "\n")
  cat("Date:", date(), "\n\n")
  
  # add the settings
  cat("Settings:\n")
  cat("Input files:", file_path, file_extension, "\n")
  cat("XCMS ppm:", xcms_ppm, "\n")
  cat("XCMS mzdiff:", xcms_mzdiff, "\n")
  cat("XCMS prefilter:", xcms_prefilter, "\n")
  cat("XCMS peakwidth:", xcms_peakwidth, "\n")
  cat("XCMS snthresh:", xcms_snthresh, "\n")
  cat("XCMS noise:", xcms_noise, "\n")
  cat("XCMS expand rt:", xcms_expand_rt, "\n")
  cat("XCMS align bin size:", xcms_align_bin_size, "\n")
  cat("XCMS grouping bw:", xcms_group_bw, "\n")
  cat("XCMS grouping bin size:", xcms_group_bin_size, "\n")
  cat("XCMS grouping min fraction:", xcms_group_min_fraction, "\n")
  cat("Peak RT window:", RTwindow, "\n")
  cat("Min intensity:", min_intensity, "\n")
  cat("ppm difference:", diff_ppm, "\n")
  cat("minimal fold change:", min_foldchange, "\n")
  cat("Isotope tolerance:", isotope_ratio_tolerance, "\n")
  cat("RT min:", round(min_RT / 60, digits = 1), "\n")
  cat("RT max:", round(max_RT / 60, digits = 1), "\n")
  cat("mz min:", min_mz, "\n")
  cat("mz max:", max_mz, "\n")
  
  cat("\n\n")
  
  # how many 13C labelling hits and sulfur hits were detected?
  cat("Total labelling hits:", nrow(labelling_hits_sulphur), "\n")
  cat("Total sulphur hits:", nrow(sulphur_hits), "\n")
  
  cat("\n\n")
  
  # report if the major expected 13C labelling events were detected
  cat("Labelling hits:\n")
  cat ("178/181 pair ", if ((labelling_hits_sulphur %>% filter(
    mzmed.x > 178.05, mzmed.x < 178.06,
    mzmed.y > 181.05, mzmed.y < 181.07, 
    rtmed_min.x > 1.9, rtmed_min.x < 2.2) %>% count) > 0) "" else "NOT ", "found\n", sep = "")
  
  cat ("162/165 pair ", if ((labelling_hits_sulphur %>% filter(
    mzmed.x > 162.05, mzmed.x < 162.065,
    mzmed.y > 165.06, mzmed.y < 165.075, 
    rtmed_min.x > 3.5, rtmed_min.x < 3.9) %>% count) > 0) "" else "NOT ", "found\n", sep = "")
  
  cat ("291/294 pair ", if ((labelling_hits_sulphur %>% filter(
    mzmed.x > 291.0, mzmed.x < 291.2,
    mzmed.y > 294.0, mzmed.y < 294.2, 
    rtmed_min.x > 6.7, rtmed_min.x < 7.2) %>% count) > 0) "" else "NOT ", "found\n", sep = "")
  
  cat ("208/212 pair ", if ((labelling_hits_sulphur %>% filter(
    mzmed.x > 208.0, mzmed.x < 208.2,
    mzmed.y > 212.0, mzmed.y < 212.2, 
    rtmed_min.x > 2.8, rtmed_min.x < 3.2) %>% count) > 0) "" else "NOT ", "found\n", sep = "")
  
  cat ("394/398 pair ", if ((labelling_hits_sulphur %>% filter(
    mzmed.x > 394.0, mzmed.x < 394.2,
    mzmed.y > 398.0, mzmed.y < 398.3, 
    rtmed_min.x > 6.4, rtmed_min.x < 7.0) %>% count) > 0) "" else "NOT ", "found\n", sep = "")
  
  cat ("180/183 pair ", if ((labelling_hits_sulphur %>% filter(
    mzmed.x > 180.0, mzmed.x < 180.2,
    mzmed.y > 183.0, mzmed.y < 183.2, 
    rtmed_min.x > 2.0, rtmed_min.x < 2.4) %>% count) > 0) "" else "NOT ", "found\n", sep = "")
  
  cat("\n\n")
  
  
  # report if the major sulfur-containing metabolites were detected
  cat("Sulphur hits:\n")
  cat ("178 ", if ((sulphur_hits %>% filter(
    mzmed.x > 178.05, mzmed.x < 178.06,
    rtmed_min.x > 1.9, rtmed_min.x < 2.2) %>% count) > 0) "" else "NOT ", "found\n", sep = "")
  
  cat ("181 (labelled) ", if ((sulphur_hits %>% filter(
    mzmed.x > 181.05, mzmed.x < 181.07,
    rtmed_min.x > 1.9, rtmed_min.x < 2.2) %>% count) > 0) "" else "NOT ", "found\n", sep = "")
  
  cat ("162 ", if ((sulphur_hits %>% filter(
    mzmed.x > 162.05, mzmed.x < 162.065,
    rtmed_min.x > 3.5, rtmed_min.x < 3.9) %>% count) > 0) "" else "NOT ", "found\n", sep = "")
  
  cat ("165 (labelled) ", if ((sulphur_hits %>% filter(
    mzmed.x > 165.05, mzmed.x < 165.075,
    rtmed_min.x > 3.5, rtmed_min.x < 3.9) %>% count) > 0) "" else "NOT ", "found\n", sep = "")
  
  cat ("291 ", if ((sulphur_hits %>% filter(
    mzmed.x > 291.0, mzmed.x < 291.2,
    rtmed_min.x > 6.7, rtmed_min.x < 7.2) %>% count) > 0) "" else "NOT ", "found\n", sep = "")
  
  cat ("294 (labelled) ", if ((sulphur_hits %>% filter(
    mzmed.x > 294.0, mzmed.x < 294.2,
    rtmed_min.x > 6.7, rtmed_min.x < 7.2) %>% count) > 0) "" else "NOT ", "found\n", sep = "")
  
  cat ("208 ", if ((sulphur_hits %>% filter(
    mzmed.x > 208.0, mzmed.x < 208.2,
    rtmed_min.x > 2.8, rtmed_min.x < 3.2) %>% count) > 0) "" else "NOT ", "found\n", sep = "")
  
  cat ("212 (labelled) ", if ((sulphur_hits %>% filter(
    mzmed.x > 212.0, mzmed.x < 212.2,
    rtmed_min.x > 2.8, rtmed_min.x < 3.2) %>% count) > 0) "" else "NOT ", "found\n", sep = "")
  
  cat ("394 ", if ((sulphur_hits %>% filter(
    mzmed.x > 394.0, mzmed.x < 394.2,
    rtmed_min.x > 6.4, rtmed_min.x < 7.0) %>% count) > 0) "" else "NOT ", "found\n", sep = "")
  
  cat ("398 (labelled) ", if ((sulphur_hits %>% filter(
    mzmed.x > 398.0, mzmed.x < 398.3,
    rtmed_min.x > 6.4, rtmed_min.x < 7.0) %>% count) > 0) "" else "NOT ", "found\n", sep = "")
  
  cat ("180 ", if ((sulphur_hits %>% filter(
    mzmed.x > 180.0, mzmed.x < 180.2,
    rtmed_min.x > 2.0, rtmed_min.x < 2.4) %>% count) > 0) "" else "NOT ", "found\n", sep = "")
  
  cat ("183 (labelled) ", if ((sulphur_hits %>% filter(
    mzmed.x > 183.0, mzmed.x < 183.2,
    rtmed_min.x > 2.0, rtmed_min.x < 2.4) %>% count) > 0) "" else "NOT ", "found\n", sep = "")
  
  sink()
}
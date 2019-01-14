parentDir <- "C:/Users/tenzl/Documents/AssaySupport/data-raw/"

test_Ablist <- readxl::read_xlsx("C:/Users/tenzl/Documents/AssaySupport/data-raw/c99_erbb4_neurexin.xlsx")

test_spectra <- AssaySupport::load_spectra(Dir = parentDir)

test_spectra_proc <- AssaySupport::preprocess_spectra(spec = test_spectra,
                                                      smooth_meth = "SavitzkyGolay",
                                                      smooth_halfWindowSize  = 5,
                                                      norm_meth = "TIC",
                                                      filter_spectra = NA,
                                                      baseline_meth = "TopHat",
                                                      align = FALSE)

test_peaks <- MALDIquant::detectPeaks(test_spectra_proc, 
                                      method = "SuperSmoother", 
                                      SNR = 3)

test_peak_df <- AssaySupport::generate_peakdf(test_spectra_proc, 
                                              pick_meth = "SuperSmoother", 
                                              SNR = 3, 
                                              binpeaks = FALSE) 


usethis::use_data(test_Ablist, test_spectra, test_spectra_proc, test_peaks, test_peak_df)

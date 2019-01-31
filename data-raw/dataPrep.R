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

# _c99_sarah.sqs: old sequence 
C99_3xFLAG <- c(C99 = "DAEFRHDSGYEVHHQKLVFFAEDVGSNKGAIIGLMVGGVVIATVIVITLVMLKKKQYTSIHHGVVEVDAAVTPEERHLSKMQQNGYENPTYKFFEQMQNGSDYKDHDGDYKDHDIDYKDDDDKGTLEVLFQ")
# C99-3xFLAG_delT.sqs: new sequence 
C99_deltaT <- c(C99T = "DAEFRHDSGYEVHHQKLVFFAEDVGSNKGAIIGLMVGGVVIATVIVITLVMLKKKQYTSIHHGVVEVDAAVTPEERHLSKMQQNGYENPTYKFFEQMQNGSDYKDHDGDYKDHDIDYKDDDDKGLEVLFQ")

# Substrates for Maria
Neurexin   <- c(Neurexin = "DAEFRVIRESSSTTGMVVGIVAAAALCILILLYAMYKYRNRDEGSYQVDETRNYISNSAQSNGTLMKEKQASSKSGHKKQKNKDKEYYVGSDYKDHDGDYKDHDIDYKDDDDKGTLEVLFQ")
ErbB4      <- c(ErbB4 = "IYYPWTGHSTLPQHARTPLIAAGVIGGLFILVIVGLTFAVYVRRKSIKKKRALRRFLETELVEPLTPSGTAYPYDVPDYASLGGPDYKDHDGDYKDHDIDYKDDDDK")
EphA4      <- c(EphA4 = "TNTVPSRIIGDGANSTVLLVSVSGSVVLVVILIAAFVISRRRSKYSKAKQEADEEKHLNQGVRTDYKDHDGDYKDHDIDYKDDDDK")


C99_mutations <- data.frame(name = c("Tottori", "Flemish", "Iowa", "Arctic", "Dutch", "Iberian"),
                            code = c("D7N", "A21G", "D23N", "E22G", "E22Q", "I45F"),
                            pos = c(7, 21, 23, 22, 22, 45),
                            substitute = c("N", "G", "N", "G", "Q", "F"), stringsAsFactors = FALSE)


usethis::use_data(test_Ablist, test_spectra, test_spectra_proc, test_peaks, test_peak_df, 
                  C99_3xFLAG, C99_deltaT, C99_mutations,
                  Neurexin, ErbB4, EphA4,
                  overwrite = TRUE)

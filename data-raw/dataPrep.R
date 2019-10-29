parentDir <- "C:/Users/tenzl/Documents/AssaySupport/data-raw/"

test_Ablist <- readxl::read_xlsx("C:/Users/tenzl/Documents/AssaySupport/data-raw/c99_erbb4_neurexin.xlsx")

test_spectra <- AssaySupport::load_spectra(Dir = parentDir)

test_spectra_proc <- AssaySupport::preprocess_spectra(spec = test_spectra,
                                                      smooth_meth = "SavitzkyGolay",
                                                      smooth_halfWindowSize  = 5,
                                                      norm_meth = "TIC",
                                                      filter_spectra = NA,
                                                      baseline_meth = "TopHat",
                                                      align = "none")

test_peaks <- MALDIquant::detectPeaks(test_spectra_proc, 
                                      method = "SuperSmoother", 
                                      SNR = 3)

test_peak_df <- AssaySupport::generate_peakdf(test_spectra_proc, 
                                              pick_meth = "SuperSmoother", 
                                              SNR = 3, 
                                              binpeaks = FALSE) 

test_result_df <- generate_resultdf(test_peak_df, 
                                    test_Ablist, 
                                    tol = 700, 
                                    tolppm = TRUE)

# _c99_sarah.sqs: old sequence 
C99_3xFLAG <- c(C99 = "DAEFRHDSGYEVHHQKLVFFAEDVGSNKGAIIGLMVGGVVIATVIVITLVMLKKKQYTSIHHGVVEVDAAVTPEERHLSKMQQNGYENPTYKFFEQMQNGSDYKDHDGDYKDHDIDYKDDDDKGTLEVLFQ")
# C99-3xFLAG_delT.sqs: new sequence 
C99_deltaT <- c(C99T = "DAEFRHDSGYEVHHQKLVFFAEDVGSNKGAIIGLMVGGVVIATVIVITLVMLKKKQYTSIHHGVVEVDAAVTPEERHLSKMQQNGYENPTYKFFEQMQNGSDYKDHDGDYKDHDIDYKDDDDKGLEVLFQ")

# Substrates for Maria
Neurexin   <- c(Neurexin = "DAEFRVIRESSSTTGMVVGIVAAAALCILILLYAMYKYRNRDEGSYQVDETRNYISNSAQSNGTLMKEKQASSKSGHKKQKNKDKEYYVGSDYKDHDGDYKDHDIDYKDDDDKGTLEVLFQ")
ErbB4      <- c(ErbB4 = "IYYPWTGHSTLPQHARTPLIAAGVIGGLFILVIVGLTFAVYVRRKSIKKKRALRRFLETELVEPLTPSGTAYPYDVPDYASLGGPDYKDHDGDYKDHDIDYKDDDDK")
EphA4      <- c(EphA4 = "TNTVPSRIIGDGANSTVLLVSVSGSVVLVVILIAAFVISRRRSKYSKAKQEADEEKHLNQGVRTDYKDHDGDYKDHDIDYKDDDDK")
#N_Cadherin <- c(N_Cadherin = "")


C99_mutations <- data.frame(name = c("Tottori", "Flemish", "Iowa", "Arctic", "Dutch", "Iberian"),
                            code = c("D7N", "A21G", "D23N", "E22G", "E22Q", "I45F"),
                            pos = c(7, 21, 23, 22, 22, 45),
                            substitute = c("N", "G", "N", "G", "Q", "F"), stringsAsFactors = FALSE)

C99T_commonSpecies <- generate_assigndf(C99_deltaT, fragmentList = list(Ab = list(start = c(1, 11), 
                                                                                  end = c(28,37:49), 
                                                                                  charge = 1), 
                                                                        AICD = list(start = 49:50, 
                                                                                    end = nchar(C99_deltaT[[1]]), 
                                                                                    charge = 1:2),
                                                                        Substrate = list(start = 1, 
                                                                                         end = nchar(C99_deltaT[[1]]), 
                                                                                         charge = 1:2)))


usethis::use_data(test_Ablist, test_spectra, test_spectra_proc, test_peaks, test_peak_df, test_result_df,
                  C99_3xFLAG, C99_deltaT, C99_mutations, 
                  C99T_commonSpecies,
                  Neurexin, ErbB4, EphA4,
                  overwrite = TRUE)

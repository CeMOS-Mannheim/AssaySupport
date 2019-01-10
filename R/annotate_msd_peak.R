#' Annotate peaks in MSD files
#'
#' @param label Character, label of the peak
#' @param file  Character, filepath to edit
#' @param peaks  a \code{MALDIquant::MassPeaks} object with corresponding peaks to file or a vector in form c(mz, intensity)
#' @param peakIdx Integer, index of peak in \code{peaks} or NA if mz and intensity are given directly
#' @param theoMz Numeric, theoretical mass of labeled peak (defaults to 0)
#' @param charge Integer, charge of peak (defaults to 1)
#' @param peakBaseline Numeric, intensity of baseline at peak position (not really sure what it does... defaults to 0)
#'
#' @return
#' writes annoteted file directly to disk
#' 
#' @export

annotate_msd_peak <- function(label, file, peaks, peakIdx = NA, theoMz = 0, charge = 1, peakBaseline = 0) {
  annotation_constructor <- function(name, peakMZ, peakIntensity, peakBaseline, charge, calcMZ) {
    anno <- paste0("<annotation peakMZ=\"", format(peakMZ, nsmall = 6, scientific = FALSE ),
                   "\" peakIntensity=\"", format(peakIntensity, nsmall = 6, scientific = FALSE), 
                   "\" peakBaseline=\"", format(peakBaseline, nsmall = 6, scientific = FALSE),
                   "\" charge=\"", charge, "\" calcMZ=\"", 
                   format(calcMZ, nsmall = 6, scientific = FALSE), "\">", name,"</annotation>")
    return(anno)
  }
  
  xml_doc <- xml2::read_xml(file, options = "")
  
  # get index of annotations node, create it if not found
  anno_idx <- which(xml2::xml_name(xml2::xml_children(xml_doc)) == "annotations")
  if(length(anno_idx) == 0) {
    xml2::xml_add_sibling(xml2::xml_children(xml_doc)[[length(xml2::xml_name(xml2::xml_children(xml_doc)))]], 
                          .value = xml2::read_xml("<annotations></annotations>"), 
                          .where = "after")
    anno_idx <- which(xml2::xml_name(xml2::xml_children(xml_doc)) == "annotations")
  }
  
  if(MALDIquant::isMassPeaks(peaks)) {
  anno <- annotation_constructor(name = label, 
                                 peakMZ = MALDIquant::mass(peaks)[peakIdx], 
                                 peakIntensity = MALDIquant::intensity(peaks)[peakIdx], 
                                 peakBaseline = peakBaseline, 
                                 charge = charge, 
                                 calcMZ = theoMz)
  } else {
    anno <- annotation_constructor(name = label, 
                                   peakMZ = peaks[1], 
                                   peakIntensity = peaks[2], 
                                   peakBaseline = peakBaseline, 
                                   charge = charge, 
                                   calcMZ = theoMz)
  }
  
  xml2::xml_add_child(xml2::xml_children(xml_doc)[[anno_idx]], .value = xml2::read_xml(anno))
  
  xml2::write_xml(xml_doc, file =file, encoding = "utf-8", options = "")
}


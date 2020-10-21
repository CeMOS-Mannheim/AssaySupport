getPlatePosition <- function(spec) {
  patch <- unlist(lapply(spec, function(x) {
    return(metaData(x)$patch)
  }))
  y <- vapply(stringr::str_sub(patch, 1, 1), function(x) {
    y <- which(LETTERS %in% x)
    return(y)
  }, FUN.VALUE = 1)
  
  x <- as.numeric(stringr::str_sub(patch, 2))
  
  return(tibble(ID = names(spec), x = x, y = y))
}

getPlateCalibarion <- function(peaks, cal_name = "CAL|Cal", ref, top_n = NA, tol = 250*1e-6, method = "quadratic") {
  nm <- names(peaks)
  if(!is.na(cal_name)) {
    nm_idx <- grepl(x = nm, pattern = cal_name)
    cat("Found", sum(nm_idx), "matches for cal_name:", cal_name, "\n")
    peaks <- peaks[nm_idx]
  }
  if(!is.na(top_n)) {
    peaks <- lapply(peaks, function(x) {
      int <- intensity(x)
      mz <- mass(x)
      idx <- order(int, decreasing = TRUE)
      
      mz_new <- mz[idx[1:top_n]]
      mz_idx <- order(mz_new)
      int_new <- int[idx[1:top_n]] 
      
      return(createMassPeaks(mass = mz_new[mz_idx], intensity = int_new[mz_idx]))
    })
    names(peaks) <- nm
  }
  
  cal <- determineWarpingFunctions(l = peaks, 
                                   reference = ref, 
                                   tolerance = tol, 
                                   method = method, 
                                   allowNoMatches = TRUE)
  return(cal)
}

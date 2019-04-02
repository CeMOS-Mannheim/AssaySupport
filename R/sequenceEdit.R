## helper to add protons 
addProton <- function(element_list, charge = 1) {
  if(is.null(element_list$H)) {
    element_list$H <- 0
  }
  
  element_list$H <- element_list$H+charge
  return(element_list)
}


#' Get the mass of ions from sequence
#' 
#' @description 
#' Caution:
#' The returned mass differs at the second decimal position from the values returned by Bruker's SequenceEditor's Charge states function!
#' This is no big deal for linear mode and should also be ok for most reflector mode applications but can not be used to determine exact masses for FT-ICR.
#'
#' @param sequence  character or (named) character vector, amino acid sequence
#' @param charge    numeric, number of charges of the ion
#' @param mode      character, can be "mono" (monoisotopic mass) or "avg" (average mass)
#' @param IAA       logical, see \code{OrgMassSpecR::ConvertPeptide}
#' @param unlist    logical, unlist result to return (named) vector of masses.
#'
#' @return list or single (named) numeric value      
#' @export
#'
#' @examples
#' 
#' getMassOfSequece(C99_3xFLAG)
getMassOfSequece <- function(sequence, charge = 1, mode = "avg", IAA = FALSE, unlist = TRUE) {
  if(!typeof(sequence) == "list") {
    sequencelist <- list(sequence)
  } else {
    sequencelist <- sequence
  }
  
  res_list <- vector("list", length(sequencelist) * length(charge))
  
  counter <- 0
  for(i in 1:length(sequencelist)) {
    for(j in 1:length(charge)) {
      counter <- counter + 1
      element_list <- addProton(OrgMassSpecR::ConvertPeptide(sequencelist[[i]], IAA = IAA), 
                                charge = charge[j])
      
      switch (mode,
              mono = {
                res <- OrgMassSpecR::MolecularWeight(element_list)
              },
              avg = {
                isotopicPattern <- OrgMassSpecR::IsotopicDistribution(element_list, 
                                                                      charge = charge[j])
                res <- isotopicPattern %>% 
                  dplyr::mutate(part = mz * percent) %>%
                  dplyr::summarise(avg = sum(part)/sum(percent)) %>%
                  dplyr::pull(avg)
              })
      res_list[[counter]] <- res
      if(!charge[j] > 1) {
        names(res_list)[counter] <- names(sequence)[i]
      } else {
        names(res_list)[counter] <- paste0(names(sequence)[i], "[", charge[j],"H+]")
      }
      
    }
  }
  if(!typeof(sequence) == "list" | unlist) {
    return(unlist(res_list))
  } else {
    return(res_list)
  }
}


#' Generate list of fragments from sequence
#'
#' @param sequence character, amino acid sequence
#' @param start    numeric vector, starting postions
#' @param end      numeric vector, end positions
#' @param prefix   character, prefix to construct name from. 
#'                 Names will be in the from "prefixStart_end"     
#'
#' @return         named list of sequence fragments
#' @export
#'
#' @examples
#' Ab42 <- "DAEFRHDSGYEVHHQKLVFFAEDVGSNKGAIIGLMVGGVVIA"
#' cutSequence(Ab42, 1, c(38, 40, 42), prefix = "Ab")
#' 
cutSequence <- function(sequence, start = 1, end = nchar(sequence), prefix = "Frag") {
  if(length(start) > 1 | length(end) >  1) {
    res_list <- vector("list", length = length(start) * length(end))
    pos = 0
    for(i in 1:length(start)) {
      for(j in 1:length(end)) {
        pos = pos + 1
        res_list[[pos]] <- substr(sequence, start[i], end[j])
        names(res_list)[pos] <- paste0(prefix, start[i], "_", end[j])
      }
    }
    return(res_list)
  }
  res <- substr(sequence, start, end)
  names(res) <- paste0(prefix, start, "_", end)
  return(res)
}

#' Generate data.frame with sequence fragments as used by \code{assign_species}
#'
#' @param sequencelist      Named list, sequences to process. Need to have the same character length for the default value of fragmentList to work
#' @param fragmentList      Named list, giving the name of the fragments and their start and end positions. See default values for an example.
#' @param noSubstrateName   Character, name for "no Substrate".
#'
#' @return
#' data.frame with columns Sequence, Species and mz
#' 
#' @examples
#' generate_assigndf(C99_3xFLAG)
#' 
#' @export

generate_assigndf <- function(sequencelist, 
                              noSubstrateName = "None",
                              fragmentList = list(Ab = list(start = 1 , end = 37:49, charge = 1),
                                                  AICD = list(start = 49:50, end = nchar(sequencelist[[1]]), charge = 1:2),
                                                  Substrate = list(start = 1, end = nchar(sequencelist[[1]]), charge = 1:2))) {
  
  res_df <- data.frame(Substrate = noSubstrateName,
                       Species = "Substrate",
                       mz = 0, 
                       stringsAsFactors = FALSE)
  
  for(i in 1:length(sequencelist)) {
    seqName <- names(sequencelist)[i]
    if(is.null(seqName)) {
      seqName <- "Unkown"
      warning("No name found for sequence number ", i,"\n")
    }
    for(j in 1:length(fragmentList)) {
      prefix <- names(fragmentList)[j]
      res <- getMassOfSequece(sequence = cutSequence(sequence = sequencelist[[i]], 
                                                     start = fragmentList[[j]][[1]], 
                                                     end = fragmentList[[j]][[2]], 
                                                     prefix = prefix), 
                              charge = fragmentList[[j]][[3]], 
                              mode = "avg", 
                              IAA = FALSE, 
                              unlist = TRUE)
      res_df <- rbind(res_df, data.frame(Substrate = seqName, Species = names(res), mz = res, stringsAsFactors = FALSE))
    }
    
  }
  rownames(res_df) <- NULL
  return(res_df)
}


#' Change single amino acid in sequence
#'
#' @param sequence   Character, sequence to edit.
#' @param pos        Integer vector, position at which \code{sequence} should be edited.
#' @param substitute Character vector, amino acid (sequece) to change the amino acid at \code{pos} with.
#' @param names      Character vector, names of generated sequences. 
#'                   Set to \code{NULL} to generate systemic names (e.g. D7N for Tottori)
#'
#' @return           mutated amino acid sequences
#'
#' @examples
#' pointMutateSequence(sequence = C99_deltaT, 
#'                     pos = C99_mutations$pos, 
#'                     substitute = C99_mutations$substitute, 
#'                     names = C99_mutations$name)
#'                     
#' @export
pointMutateSequence <- function(sequence, pos, substitute, names = NULL) {
  if(!length(pos) == length(substitute)) {
    stop("pos and substitute must have same length!\n")
  }
  
  res_list <- vector("list", length(pos))
  counter <- 0
  for(j in 1:length(pos)) {
    
    counter <- counter + 1
    start <- substr(sequence, 1, pos[j] - 1)
    end <- substr(sequence, pos[j] + 1, nchar(sequence))
    res_list[[counter]] <- paste0(start, substitute[j], end)
    
    if(any(is.null(names))) {
      names(res_list)[counter] <- paste0(substr(sequence, pos[j], pos[j]), pos[j], substitute[j])
    } else {
      if(!length(names) == length(res_list)) {
        stop("There has to be either no (automatic naming) or a name for each mutated sequence! \n")
      }
      names(res_list)[counter] <- names[counter]
    }
    
  }
  if(length(pos) == 1 | length(substitute) == 1){
    return(unlist(res_list))
  } else {
    return(res_list)
  }
  
}

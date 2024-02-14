#this function create a vector of spacers based on the wanted number of guides that, within each guide, respect the minimum distance given by the user. 
#from the pool of sequences given by the prediction and ranked by efficiency it takes one at a time a spacer that respect the parameters mentioned above, until it takes all the spacers needed
#this version check also if the spacer is extendable (it would not be in case it's at the very limit of the transcript and it would have an NA in the final table as a mistake)
remove_redundancy_array_v3 <- function(sequence ,spacers_pool, min_dist_array, guide_length, n_fastas, n_guides, name, min_dist_array_global, extension_length, restriction_sequences) {
  impossible <- 0
  #sequences of restriction enzymes to remove (to put in the inputs and remove from hardcoding)
  # restriction_sequences <- c("CGTCTC")
  #check that first spacer taken doesn't contain enzyme sequence
  spacer_position <- str_locate(sequence, revComp(spacers_pool[1]))#get position of reverse complemented spacer
  extended_spacer <- str_sub(sequence, spacer_position[1]-extension_length, (spacer_position[2]))#extract the reverse complemented spacer but with the extension
  extended_spacer <- revComp(extended_spacer)#reverse complement to spacer again
  #loop until finding first spacer without the enzyme sequence
  while (str_detect(extended_spacer, paste(restriction_sequences, collapse = "|"))) {
    spacers_pool <- list.remove(spacers_pool, 1)#removes the first n_guides spacers from the spacers_pool in order to avoid taking them again later in the program
    spacer_position <- str_locate(sequence, revComp(spacers_pool[1]))#get position of reverse complemented spacer
    extended_spacer <- str_sub(sequence, spacer_position[1]-extension_length, (spacer_position[2]))#extract the reverse complemented spacer but with the extension
    extended_spacer <- revComp(extended_spacer)#reverse complement to spacer again
  }
  #loop until finding first extendable spacer in case first spacers cannot be extended because at the end of the sequence
  while (nchar(extended_spacer) < 30 && length(spacers_pool != 0)) {
    spacers_pool <- list.remove(spacers_pool, 1)#removes the first n_guides spacers from the spacers_pool in order to avoid taking them again later in the program
    spacer_position <- str_locate(sequence, revComp(spacers_pool[1]))#get position of reverse complemented spacer
    extended_spacer <- str_sub(sequence, spacer_position[1]-extension_length, (spacer_position[2]))#extract the reverse complemented spacer but with the extension
    extended_spacer <- revComp(extended_spacer)#reverse complement to spacer again
  }
  spacers <- head(spacers_pool, 1)#takes the first n_guides spacers away from the spacers pool
  min_dist <- head(as.integer(min_dist_array), 1)
  min_dist_array <<- list.remove(min_dist_array, 1)
  min_dist_global <- head(as.integer(min_dist_array_global), 1)
  min_dist_array_global <<- list.remove(min_dist_array_global, 1)
  #add to say that if percentage minimum global distance results in less than 30 bp you have to take 23 as minimum global distance
  if (min_dist_global < 23) {
    min_dist_global <- 23
  }
  if (min_dist < 23) {
    min_dist <- 23
  }
  spacers_pool <- list.remove(spacers_pool, 1)#removes the first n_guides spacers from the spacers_pool in order to avoid taking them again later in the program
  
  print(paste("min_dist:",min_dist, "min_dist_global:", min_dist_global))
  for (n in 1:(n_fastas-1)) {#loops through n_fasta-n_guides times to take the spacers that are needed
    exit <- 0
    counter <- 1
    while ((counter < (length(spacers_pool))) && (exit == 0)) {#loops through all the spacers pool
      parsing <- length(spacers)#keeps track of the parsing position
      failure <- 0#flag to see if the spacers doesn't meet the requirements
      count_parse <- 1
      while ((parsing > 0) && (failure == 0)) {
        pos1 <- str_locate(sequence, revComp(spacers_pool[counter]))#get start and end position of one spacer
        pos2 <- str_locate(sequence, revComp(spacers[parsing]))#get start and end position of subsequent spacer
        dist <- abs(pos2[1]-pos1[1]) #distance between end of one spacer and start of subsequent spacer
        if ((count_parse %% n_guides == 0) && (dist <= min_dist) && (dist <= min_dist_global)) {#checks for in guides distance
          failure <- 1
        }
        if ((count_parse %% n_guides != 0) && (dist <= min_dist_global)) {#checks for global distance
          failure <- 1
        }
        #enzyme sequence test
        spacer_position <- str_locate(sequence, revComp(spacers_pool[counter]))#get position of reverse complemented spacer
        extended_spacer <- str_sub(sequence, spacer_position[1]-extension_length, (spacer_position[2]))#extract the reverse complemented spacer but with the extension
        extended_spacer <- revComp(extended_spacer)#reverse complement to spacer again
        if (str_detect(extended_spacer, paste(restriction_sequences, collapse = "|"))) {
          failure <- 1
        }
        #test that the spacer is extendable
        if (nchar(extended_spacer) < 30) {
          failure <- 1
        }
        parsing <- parsing-1#it makes sure that we check the distance for all already set spacers
        count_parse <- count_parse + 1#is used to check if we need to use global distance or in guide distance
      }
      if (failure == 0) {
        spacers <- append(spacers, spacers_pool[counter])#put distant enough spacer in the spacers list
        spacers_pool <- list.remove(spacers_pool, counter)#remove spacer from spacers_pool to avoid taking it multiple times
        exit <- 1
      }
      counter <- counter + 1
      if (counter == (length(spacers_pool)-1)) {
        impossible <- 1 #flag that goes on if there are not enough guides with the specified distance
      }
    }
  }
  if (impossible == 1) {
    print(paste0(name, ": Impossible to find enough spacers with the suggested minimum distance, try to lower the stringency of the parameters"))
  }
  return(spacers)
}

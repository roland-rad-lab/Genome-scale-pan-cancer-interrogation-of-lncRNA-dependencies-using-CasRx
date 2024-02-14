#this function extend every spacer by the number of bases specified by the user thorugh the "extension" parameter
extend_spacers <- function(sequence ,spacers, extension) {
  extended_spacers <- lapply(spacers, function(x){
    spacer_position <- str_locate(sequence, revComp(x))#get position of reverse complemented spacer
    extended_spacer <- str_sub(sequence, spacer_position[1]-extension, (spacer_position[2]))#extract the reverse complemented spacer but with the extension
    extended_spacer <- revComp(extended_spacer)#reverse complement to spacer again
    return(extended_spacer)
  })
  return(extended_spacers)
}
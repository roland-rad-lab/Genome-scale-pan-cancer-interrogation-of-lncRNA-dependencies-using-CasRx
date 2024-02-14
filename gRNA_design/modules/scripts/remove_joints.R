#this function keep only spacers that are fully contained in the exons of the targeted transcript
remove_joints <- function(sequence ,file, ids_by_exon) {
  spacers <- vector()
  spacers_pool <- as.vector(file$GuideSeq)
  for (spacer in spacers_pool) {
    for (exon in ids_by_exon$sequence) {
      if (grepl(revComp(spacer), exon)) {
        spacers <- append(spacers, spacer)
      }
    }
  }
  return(spacers)
}
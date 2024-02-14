trim_sides <- function(table, left_trim, right_trim){
  sequence <- apply(table, 1, function(x){
    length <- nchar(x[2])
    loss_5 <- trunc(length*left_trim)
    loss_3 <- trunc(length*right_trim)
    x[2] <- str_sub(x[2], start = loss_5, end = loss_3)
  })
  table["sequence"] <- as.data.frame(sequence)
  table
}


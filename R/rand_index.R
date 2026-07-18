#' Rand Index between two partitions
#'
#' Computes the  Rand Index (RI) between two partitions of the same
#' set of objects.
#'
#' @param part1,part2 Two vectors of group labels (numeric, character, or
#'   factor), of the same length, giving the partition of each object under
#'   two different classifications.
#'
#' @return A single numeric value (RI), equal to 1 for identical partitions
#'   and close to 0 for partitions no more similar than chance.
#'
#' @export 
rand_index<-function(part1,part2){
  part1<-dissim_partition(as.matrix(part1))
  part2<-dissim_partition(as.matrix(part2))
  p=dim(part1)[1]
return((sum(part1==part2)-p)/(p*(p-1)))
}

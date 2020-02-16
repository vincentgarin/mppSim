###########
# het_det #
###########

# function to determine the heterozygous score of parents

het_det <- function(x){

  a1 <- sample(strsplit(x = x[1], '')[[1]], size = 1)
  a2 <- sample(strsplit(x = x[2], '')[[1]], size = 1)

  paste0(a1, a2)

}

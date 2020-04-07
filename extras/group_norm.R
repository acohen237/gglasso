# Input: a vector of length n, another vector of length m
# The first vector is assumed to be partitioned (consecutively)
# into groups. The second vector describes the starting index of
# each group within the first vector

group_norm <- function(v, vg) {
  s <- 0
  for (i in 1:(length(vg)-1)) {
    subvec <- v[vg[i]:(vg[i+1]-1)]
    s <- s + sqrt(sum(subvec^2))
  }
  last_subvec <- v[tail(vg, n = 1):length(v)]
  s <- s + sqrt(sum(last_subvec^2))
  return(s)
}



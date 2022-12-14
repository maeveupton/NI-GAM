#--------------------THIS SECTION IS DATA CLEANING DONE IN OTHER CODE------------------------
#------Match Closest function doesn't exist on version R 3.6.3----
#--Will be used when matching similar long & lat values from both datasets--
match.closest <- function(x, table, tolerance=Inf, nomatch=NA_integer_) {
  lIdx <- findInterval(x, table, rightmost.closed=FALSE, all.inside=TRUE)
  rIdx <- lIdx + 1L
  lIdx[lIdx == 0L] <- 1L
  lDiff <- abs(table[lIdx] - x)
  rDiff <- abs(table[rIdx] - x)
  d <- which(lDiff >= rDiff)
  lIdx[d] <- rIdx[d]
  if (any(is.finite(tolerance))) {
    if (any(tolerance < 0L)) {
      warning(sQuote("tolerance"), " < 0 is meaningless. Set to zero.")
      tolerance[tolerance < 0L] <- 0L
    }
    if (length(nomatch) != 1L) {
      stop("Length of ", sQuote("nomatch"), " has to be one.")
    }
    tolerance <- rep_len(tolerance, length(table))
    lDiff[d] <- rDiff[d]
    lIdx[lDiff > tolerance[lIdx]] <- nomatch
  }
  lIdx
}

# Checking the data used in the different functions of DarkDiv
#
# \code{dataPrep} Checks the validity of the data provided for dark diversity estimations. It basically checks that the names of species are similar in the reference (if provided) and testing datasets.
#
#
# @param x Study data with sites in rows and species in columns. Can be a matrix or a data.set. Names for species must be provided.
# @param r Optional dataset for reference, with sites in rows and species in columns. If it is not provided, x is used by default. If it is provided, the names of the species must coincide with the names in x.
# @param removeAbsent Logical indicating what to do with species with zero occurrences in the indication matrix (i.e. for which no indication values can be estimated). removeAbsent = TRUE indicates that these species should be removed from results (giving dark diversity and pool matrices whose dimensions might not coincide with x). removeAbsent = FALSE indicates that these columns will be kept in the results, but filled with NAs. Default to TRUE.
# @param wa Logical indicating whether abundance should be considered for estimatiosn of dark diversity. Currently methods are only developed for wa=F, and setting wa=T will result in ending the function with an error. Defaults to FALSE


dataPrep <- function(x, r = x, removeAbsent, wa = F){
  x <- xOrig <- as.matrix(x)
  r <- rOrig <- as.matrix(r)
  if (is.null(colnames(x)) | is.null(colnames(r))) {
    stop("x and/or r do not have colnames. Column names must be the name of species.")
  }
  if (!identical(x,r)){
    if (!any(colnames(x) %in% colnames(r))){
      stop("x and/or r do not have common colnames. At least some of the species in x should be present in r.")
    }
    if (!all(colnames(x) %in% colnames(r))){
      warning("x does not contain exactly the same species as r.
              \nOnly those species present in r have been kept in x")
      x <- x[, colnames(r)]
    }
  }
  if (any(colSums(r>0) == 0)){
    removeSp <- which(colSums(r>0) == 0)
    r <- r[, -removeSp]
    x <- x[, -removeSp]
    if(removeAbsent){
      warning(paste("r included", length(removeSp), "species with zero occurrences.
                  \nThey have been removed from both r and x"))
    }else{
      warning(paste("r included", length(removeSp), "species with zero occurrences.
                  \nThey have been filled with NA"))
    }
  }
  return(list(x = x, r = r, xOrig = xOrig, rOrig = rOrig))
}

# \code{coocPrep} Makes the co-occurrence matrix needed for all the dark diversity methods used here
#
#
# @param r Dataset for reference, with sites in rows and species in columns.

coocPrep <- function(r){
  r <- ifelse(r > 0, 1, 0)
  ## co-occurrences
  M <- crossprod(r, r)
  return(M)
}


# \code{Beals} Estimates the Beals index (Beals smoothing) of dark diversity. The Beals index gives the probability that a species \emph{i} occurs in a site, based on the identity of the other species that are locally present, and their patterns of co-occurrence. For each pair of species i and j, the conditional probability that i occurs given that j is present is estimated and stored in the so called "indication matrix". Afterwards, these conditional probabilities are averaged for each site considering the species that were locally present, so that the result is a matrix containing the averaged probability for all species in all sites. The estimation of the Beals index presented here is based on the \code{beals} function of the package \code{vegan}.
#
#
# @param M Matrix containing the number co-occurrences between pairs of species.
# @param r Dataset for reference, with sites in rows and species in columns.
# @param x Study data with sites in rows and species in columns.

Beals <- function(M, x, r){
  C <- diag(M)
  M <- sweep(M, 2, replace(C, C == 0, 1), "/")
  diag(M) <- 0
  M <- as.matrix(M)
  S <- rowSums(r)
  b <- r
  for (i in 1:nrow(r)){
    b[i, ] <- rowSums(sweep(M, 2, r[i, ], "*"))
  }
  SM <- rep(S, ncol(r))
  SM <- SM - r
  b <- b/replace(SM, SM == 0, 1)
  t <- b
  colnames(M) <- paste("I", colnames(M), sep=".")
  rownames(M) <- paste("T", rownames(M), sep=".")
  ## for study set
  if (!identical(x,r)) {
    S <- rowSums(x)
    b <- x
    for (i in 1:nrow(x)) {
      b[i, ] <- rowSums(sweep(M, 2, x[i, ], "*"))
    }
    SM <- rep(S, ncol(x))
    SM <- SM - x
    b <- b/replace(SM, SM == 0, 1)
  }
  results <- list(indication = M, AllProbs = b, Pool = replace(b, x>0, 1),
                  Dark = replace(b, x > 0, NA), t = t)
  return(results)
}

# \code{BealsThres} Applies thresholds to raw beals values to transform the probabilities for each species in each site into binary presence/absence indication.
#
#
# @param Beals list coming from the Beals function.
# @param limit can be "quantile", "min", "const", "outlier" to determine threshold for species inclusion
# @param const constant for limit (as quantile or as minimal)
# @param r Dataset for reference, with sites in rows and species in columns.
# @param x Study data with sites in rows and species in columns.

BealsThres <- function(Beals, limit = NULL, const = 0.01, r, x){

  if (limit == "quantile" | limit == "const"){
    if(is.null(const) | const <0 | const>1){
      stop(paste("For limit of type", limit, "'const' must be specified and ranging between 0 and 1"))
    }
  }
  ## thresholds for each specis
  q <- numeric(ncol(r))
  t <- Beals$t
  b <- Beals$Pool
  for (j in 1:ncol(r)){
    q[j] <- 1
    if (limit == "min" & sum(r[, j]) > 0){
      q[j] <- min(t[r[, j] > 0, j])
    }
    if (limit == "outlier" & sum(r[, j]) > 0){
      q[j] <- stats::quantile(t[r[, j] > 0, j], probs = 0.25) - 1.5 * stats::IQR(t[r[, j] > 0, j])
      if (q[j] < min(t[r[, j] > 0, j])){
        q[j] <- min(t[r[, j] > 0, j])  ## outlier removal cannot be less than min
      }
    }
    if (limit == "quantile" & sum(r[, j]) > 0){
      q[j] <- stats::quantile(t[r[, j] > 0, j], probs = const, na.rm = T)
    }
    if (limit == "const" & sum(r[, j]) > 0){
      q[j] <- const
    }
  }
  threshold <- matrix(q, ncol=1, dimnames = list(colnames(r), "threshold"))
  ## Create copy for dark diversity estimation
  DD <- x
  DD[,] <- 0
  ## determining dark diversity
  for (i in 1:nrow(x)){
    DD[i, b[i, ] >= q & x[i, ] == 0] <- 1
  }
  return(list(indication = Beals$indication, AllProbs = DD,  Pool = DD + x,
              Dark = replace(DD, x > 0, NA)))
}



####FAVORABILITY
Favorability <- function(Beals, x){
  ##P is Beals raw index for each site & species
  ##n0 is the number of absences of that species
  ##n1 is the number of presences
  P <- Beals$Pool
  n1 <- matrix(rep(colSums(x > 0), nrow(x)), nrow=nrow(x), byrow=T)
  n0 <- nrow(x) - n1
  DD <- (P / (1 - P))/((n1 / n0) + (P / (1 - P)))
  return(list(indication = Beals$indication, AllProbs = DD, Pool = replace(DD, x>0, 1),
              Dark = replace(DD, x > 0, NA)))
}


# \code{Hypergeometric} Estimates dark diversity probability based on teh hypergeometric distribution.
#
#
# @param M Matrix containing the number co-occurrences between pairs of species.
# @param r Dataset for reference, with sites in rows and species in columns.
# @param x Study data with sites in rows and species in columns.

Hypergeometric <- function(M, x, r){
  C <- diag(M)
  N <- nrow(r) #total number of plots
  S <- length(C) #total number of species
  M1 <- matrix(rep(C, S), nrow = S)
  M <- stats::as.dist(M) #co-occurrence matrix
  M2 <- stats::as.dist(t(M1)) #total occurences of each species (in columns)
  M1 <- stats::as.dist(M1) #total occurences of each species (in rows)
  Mhat <- M1 * M2 / N # Expected number of occurences between pairs of species
  variance <- (M1 * M2 / N) * (N - M1) / N * (N - M2) / (N - 1)
  sdev <- sqrt(variance)
  M <- (M - Mhat) / sdev #Standardized Effect Size
  M <- replace(M, sdev == 0, 0)
  M <- as.matrix(M)
  # zero species cannot have indication 1, replacing by 0
  M[as.matrix(M1 * M2) == 0] <- 0
  colnames(M) <- paste("I",colnames(M),sep=".")
  rownames(M) <- paste("T",rownames(M),sep=".")
  ###Dark diversity matrix: Probability of absent species being in the local pool
  #Transform the indication values into probabilities:
  M <- stats::pnorm(M, mean=0, sd=1)
  # Average probabilities:
  b <- (x %*% M / rowSums(x))
  dimnames(b) <- dimnames(x)
  results <- list(indication = M, AllProbs = b, Pool = replace(b, x>0, 1),
                  Dark = replace(b, x > 0, NA))
  return(results)
}










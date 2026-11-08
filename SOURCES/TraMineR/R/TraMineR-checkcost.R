checkcost <- function(sma, seqdata, with.missing, indel, tol = NULL) {

  msg("checking 'sm' (size and triangle inequality)")

  if (max(sma) <= 0) {
    msg.warn("max substitution cost less or equal to 0!")
    return()
  }

  chk.zero.sm <- TRUE
  chk.triangleineq <- TRUE

  if (is.null(tol)) {
    tol <- 1e-07 * max(sma)
  }

  alphabet <- attr(seqdata, "alphabet")
  ## Gaps in sequences
  if (with.missing) {
    alphabet <- c(alphabet, attr(seqdata, "nr"))
  }
  alphsize <- length(alphabet)

  if (length(dim(sma)) == 2) {
    dim(sma) <- c(dim(sma), 1)
  } else {
    if (ncol(seqdata) != dim(sma)[3]) {
      stop(" [!] size of substitution cost matrix must be ", alphsize, "x",
        alphsize, "x", ncol(seqdata), call. = FALSE)
    }
  }

  if (dim(sma)[3] == 1) { #indel not used for time-varying costs

      if (!missing(indel) && any(indel < 0)) {
        stop(" [!] negative indel value(s)!")
      }

      if (!missing(indel) && any(indel == 0)) {
        msg.warn("[!] At least one indel value is 0!")
      }
  }

  for (i in 1:dim(sma)[3]) {
    prefix <- if (dim(sma)[3] > 1) paste0("At pos ",i,",") else ""
    sm <- sma[, , i]
    if (nrow(sm) != alphsize | ncol(sm) != alphsize) {
      stop(" [!] size of substitution cost matrix must be ", alphsize, "x",
        alphsize, call. = FALSE)
    }

    if (any(sm < 0)) {
      stop(" [!] Negative substitution costs are not allowed", call. = FALSE)
    }

    distsm <- as.dist(sm)
    if (chk.zero.sm && any(distsm == 0)) {
      kzero <- which(distsm == 0)
      ijzero <- rowcol(kzero,distsm) #rowcol defined below
      lkz <- length(kzero)
      isare <- if (lkz==1) "is 0!" else "are 0!"
      msg.warn("[!]", prefix, lkz, "substitution cost value(s)", isare)
      msg.warn("     for example, cost is zero for the pair (",
                alphabet[ijzero[1,2]],",",alphabet[ijzero[1,1]],")")
    
      chk.zero.sm <- FALSE
    }
    rm(distsm)

    if (any(diag(sm) != 0)) {
      stop(" [!] Non-zero diagonal elements in substitution cost matrix (sm)!")
    }

    if (chk.triangleineq){
        triangleineq <- checktriangleineq(sm, warn = FALSE, indices = TRUE, tol = tol)

        ## triangleineq contains a vector of problematic indices.
        if (!is.logical(triangleineq)) {
          message(" [!!] ",prefix," at least, one substitution cost doesn't respect the triangle inequality.\n",
            " [!!] replacing ", alphabet[triangleineq[1]], " with ", alphabet[triangleineq[3]],
            " (cost=", format(sm[triangleineq[1], triangleineq[3]]), ") and then ",
            alphabet[triangleineq[3]], " with ", alphabet[triangleineq[2]], " (cost=",
            format(sm[triangleineq[3], triangleineq[2]]), ")\n [!] costs less than replacing directly ",
            alphabet[triangleineq[1]], " with ", alphabet[triangleineq[2]], " (cost=",
            format(sm[triangleineq[1], triangleineq[2]]), ")\n", " [!] total difference ([",
            alphabet[triangleineq[1]], "=>", alphabet[triangleineq[3]], "] + [",
            alphabet[triangleineq[3]], "=>", alphabet[triangleineq[2]], "] - [",
            alphabet[triangleineq[1]], "=>", alphabet[triangleineq[2]], "]): ",
            format(sm[triangleineq[1], triangleineq[3]] + sm[triangleineq[3],
              triangleineq[2]] - sm[triangleineq[1], triangleineq[2]]))

          chk.triangleineq <- FALSE # avoid multiple messages in the time-varying case
        }
    }

    if (!missing(indel) && dim(sma)[3]==1) {
      if (length(indel) > 1) {
        if (length(indel) != alphsize) {
          stop(" [!] You must specify either a single indel or one indel per state.")
        }
        sm2 <- cbind(sm, indel)
        sm2 <- rbind(sm2, c(indel, 0))
        triangleineq <- checktriangleineq(sm2, warn = FALSE, indices = TRUE, tol = tol)
        ## triangleineq contains a vector of problematic indices.
        if (!is.logical(triangleineq)) {
          msg.warn("at least one indel cost does not respect the triangle inequality.")
        }
      } else if (any(sm > 2 * indel)) {
        msg.warn("At least one substitution cost greater than twice the indel cost.\n",
          "     Such substitution costs will never be used.")
      }
    }

    ## Testing for symmetric matrix
    if (any((sm - t(sm)) > tol)) {
      warning("Non-symmetric substitution cost matrix!", call. = FALSE)
    }
  }
}

## 1D index to 2D index
rowcol <- function (k, dist_obj) {
  #if (!inherits(dist_obj, "dist")) stop("please provide a 'dist' object")
  n <- attr(dist_obj, "Size")
  valid <- (k >= 1) & (k <= n * (n - 1) / 2)
  k_valid <- k[valid]
  j <- rep.int(NA_real_, length(k))
  j[valid] <- floor(((2 * n + 1) - sqrt((2 * n - 1) ^ 2 - 8 * (k_valid - 1))) / 2)
  i <- j + k - (2 * n - j) * (j - 1) / 2
  cbind(i, j)
}

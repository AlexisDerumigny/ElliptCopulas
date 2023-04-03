#' Fast estimation of Kendall's tau matrix
#'
#' Estimate Kendall's tau matrix using averaging estimators. Under
#' the structural assumption that Kendall's tau matrix is block-structured
#' with constant values in each off-diagonal block, this function estimates
#' Kendall's tau matrix ``fast'', in the sense that each interblock
#' coefficient is estimated in time \eqn{N \cdot n \cdot log(n)},
#' where \code{N} is the amount of pairs that are averaged.
#'
#'
#' @param dataMatrix matrix of size \code{(n,d)} containing \code{n} observations
#' of a \code{d}-dimensional random vector.
#'
#' @param averaging type of averaging used for fast estimation.
#' Possible choices are \itemize{
#'   \item \code{no}: no averaging;
#'   \item \code{all}: averaging all Kendall's taus in each block.
#'   \code{N} is then the number of entries in the block, i.e. the
#'   products of both dimensions.
#'   \item \code{diag}: averaging along diagonal blocks elements.
#'   \code{N} is then the minimum of the block's dimensions.
#'   \item \code{row}: averaging Kendall's tau along the smallest block side.
#'   \code{N} is then the minimum of the block's dimensions.
#'   \item \code{random}: averaging Kendall's taus along a random sample of \code{N} entries
#'   of the given block. These entries are chosen uniformly without replacement.
#' }
#'
#' @param blockStructure list of vectors.
#' Each vector corresponds to one group of variables
#' and contains the indexes of the variables that belongs to this group.
#' \code{blockStructure} must be a partition of \code{1:d},
#' where \code{d} is the number of columns in \code{dataMatrix}.
#'
#' @param N number of entries to average (n the random case.
#' By default, \code{N} is then the minimum of the block's dimensions.
#'
#'
#' @return matrix with dimensions depending on \code{averaging}.
#' \itemize{
#'   \item If \code{averaging = no},
#'   the function returns a matrix of dimension \code{(n,n)}
#'   which estimates the Kendall's tau matrix.
#'
#'   \item Else, the function returns a matrix of dimension
#'   \code{(length(blockStructure) , length(blockStructure))}
#'   giving the estimates of the Kendall's tau for each block with ones on the diagonal.
#'
#' }
#'
#' @author Rutger van der Spek, Alexis Derumigny
#' @concept Kendall correlation coefficient
#'
#' @references van der Spek, R., & Derumigny, A. (2022).
#' Fast estimation of Kendall's Tau and conditional Kendall's Tau matrices under structural assumptions.
#' \href{https://arxiv.org/abs/2204.03285}{arxiv:2204.03285}.
#'
#' @examples
#' # Estimating off-diagonal block Kendall's taus
#' matrixCor = matrix(c(1  , 0.5, 0.3 ,0.3, 0.3,
#'                      0.5,   1, 0.3, 0.3, 0.3,
#'                      0.3, 0.3,   1, 0.5, 0.5,
#'                      0.3, 0.3, 0.5,   1, 0.5,
#'                      0.3, 0.3, 0.5, 0.5,   1), ncol = 5 , nrow = 5)
#' dataMatrix = mvtnorm::rmvnorm(n = 100, mean = rep(0, times = 5), sigma = matrixCor)
#' blockStructure = list(1:2, 3:5)
#' estKTMatrix = list()
#' estKTMatrix$all = KTMatrixEst(dataMatrix = dataMatrix,
#'                               blockStructure = blockStructure,
#'                               averaging = "all")
#' estKTMatrix$row = KTMatrixEst(dataMatrix = dataMatrix,
#'                               blockStructure = blockStructure,
#'                               averaging = "row")
#' estKTMatrix$diag = KTMatrixEst(dataMatrix = dataMatrix,
#'                                blockStructure = blockStructure,
#'                                averaging = "diag")
#' estKTMatrix$random = KTMatrixEst(dataMatrix = dataMatrix,
#'                                  blockStructure = blockStructure,
#'                                  averaging = "random", N = 2)
#' InterBlockCor = lapply(estKTMatrix, FUN = function(x) {sin(x[1,2] * pi / 2)})
#'
#' # Estimation of the correlation between variables of the first group
#' # and of the second group
#' print(unlist(InterBlockCor))
#' # to be compared with the true value: 0.3.
#'
#' @export
#'
KTMatrixEst <- function(dataMatrix, blockStructure = NULL, averaging = "no", N = NULL)
{
  d = ncol(dataMatrix)
  n = nrow(dataMatrix)

  if (averaging == "no"){
    estimate <- wdm::wdm(dataMatrix, method = "kendall")
    return(estimate)
  }

  # We now assume that one of the averaging method is used.

  if (is.null(blockStructure))
  {
    stop("To use averaging estimators, a block structure must be specified.")
  }
  if (length(blockStructure) == 1){
    stop("To use averaging estimators, there must be more than one block.")
  }
  if ( sum(1:d %in% unlist(blockStructure)) != d | length(unlist(blockStructure)) != d)
  {
    stop(paste0("The block structure is not a partition of 1:", d ))
  }

  totalGroups = length(blockStructure)
  estimate = matrix(data = 1 , nrow = totalGroups , ncol = totalGroups)

  switch(
    averaging,

    "diag" = {

      for (g1 in 2:totalGroups) {
        for (g2 in 1:(g1-1)) {
          diagSize = min(length(blockStructure[[g1]]) ,
                         length(blockStructure[[g2]]) )

          vectorBlockKT = rep(NA, times = diagSize)
          for (j in 1:diagSize)
          {
            vectorBlockKT[j] = wdm::wdm(
              x = dataMatrix[ , blockStructure[[g1]][j] ] ,
              y = dataMatrix[ , blockStructure[[g2]][j] ] ,
              method = "kendall")
          }
          blockKT = mean(vectorBlockKT)
          estimate[g1,g2] = blockKT
          estimate[g2,g1] = blockKT
        }
      }
    } ,

    "all" = {
      for (g1 in 2:totalGroups) {
        for (g2 in 1:(g1-1)) {
          matrixBlockKT = matrix(nrow = length(blockStructure[[g1]]) ,
                                 ncol = length(blockStructure[[g2]]) )
          for (j1 in 1:length(blockStructure[[g1]]) )
          {
            for (j2 in 1:length(blockStructure[[g2]]) )
            {
              matrixBlockKT[j1,j2] = wdm::wdm(
                x = dataMatrix[ , blockStructure[[g1]][j1] ] ,
                y = dataMatrix[ , blockStructure[[g2]][j2] ] ,
                method = "kendall")
            }
          }
          blockKT = mean(matrixBlockKT)
          estimate[g1,g2] = blockKT
          estimate[g2,g1] = blockKT
        }
      }
    } ,

    "row" = {

      for (g1 in 2:totalGroups) {
        for (g2 in 1:(g1-1)) {
          g1Size = length(blockStructure[[g1]])
          g2Size = length(blockStructure[[g2]])

          rowSize = min(g1Size, g2Size)

          gSmall = ifelse(g1Size <= g2Size, g1, g2)
          gLarge = ifelse(g1Size <= g2Size, g2, g1)

          vectorBlockKT = rep(NA, times = rowSize)

          for (j in 1:rowSize)
          {
            vectorBlockKT[j] = wdm::wdm(
              x = dataMatrix[ , blockStructure[[gSmall]][j] ] ,
              y = dataMatrix[ , blockStructure[[gLarge]][1] ] ,
              method = "kendall")
          }
          blockKT = mean(vectorBlockKT)
          estimate[g1,g2] = blockKT
          estimate[g2,g1] = blockKT
        }
      }
    } ,

    "random" = {
      for (g1 in 2:totalGroups) {
        for (g2 in 1:(g1-1)) {
          matrixBlockKT = matrix(nrow = length(blockStructure[[g1]]) ,
                                 ncol = length(blockStructure[[g2]]) )
          numberEntries = length(matrixBlockKT)
          if (is.null(N)){
            size = min(length(blockStructure[[g1]]) , length(blockStructure[[g2]]))
          } else {
            size = N
          }
          selectedEntries = sample.int(n = numberEntries, size = size, replace = FALSE)
          for (j in 1:size)
          {
            j1 = (selectedEntries[j] - 1) %% length(blockStructure[[g1]]) + 1
            j2 = (selectedEntries[j] - 1) %/% length(blockStructure[[g1]]) + 1
            matrixBlockKT[j1,j2] = wdm::wdm(
              x = dataMatrix[ , blockStructure[[g1]][j1] ] ,
              y = dataMatrix[ , blockStructure[[g2]][j2] ] ,
              method = "kendall")
          }
          blockKT = mean(matrixBlockKT, na.rm = TRUE)
          estimate[g1,g2] = blockKT
          estimate[g2,g1] = blockKT
        }
      }
    }

  )

  if (!is.null(names(blockStructure))){
    colnames(estimate) <- rownames(estimate) <- names(blockStructure)
  }

  return(estimate)
}

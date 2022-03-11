#' Estimate the (averaging) Kendall's tau matrix
#'
#' @param dataMatrix n x d matrix of d-dimensional multivariate random variable
#' of n timepoints.
#'
#'
#' @param typeMatrixKT name of the matrix estimator used. Possible choices are
#' "all" (no averaging),
#' "aveDiag" (averaging along diagonal block elements),
#' "aveAll" (averaging all KT's within blocks),
#' "aveRow" (averaging elements on the first row along the smallest block side)
#' @param blockStructure list of groups. A group is a vector with
#' variable numbers. \code{blockStructure} must be a partition of 1:d, where d is
#' the amount of columns in \code{dataMatrix}.
#'
#'
#' @return matrix with dimensions depending on \code{typeMatrixKT}.
#' If \code{typeMatrixKT} = "all", the function returns a matrix of dimension
#' n x n, giving the Kendall's tau matrix.
#' Here, n is the number of rows in \code{dataMatrix}.
#' If \code{typeMatrixKT} = "aveDiag" or "aveAll" the function returns a matrix
#' of dimension \code{length(blockStructure)} x \code{length(blockStructure)},
#' giving the block estimates of the Kendall's tau with ones on the diagonal.
#'
#' @author Rutger van der Spek, Alexis Derumigny
#'
#' @examples
#' # Estimating off-diagonal block Kendall's taus
#' matrixCov = matrix(c(1, 0.5, 0.3 ,0.3,
#'                      0.5, 1, 0.3, 0.3,
#'                      0.3, 0.3, 1, 0.5,
#'                      0.3, 0.3, 0.5, 1), ncol = 4 , nrow = 4)
#' dataMatrix = mvtnorm::rmvnorm(n = 100, mean = rep(0, times = 4), sigma = matrixCov)
#' blockStructure = list(1:2, 3:4)
#' KTMatrixEst(dataMatrix = dataMatrix, blockStructure = blockStructure,
#'             typeMatrixKT = "aveAll")
#'
#'
#' @export
#'
KTMatrixEst <- function(dataMatrix, blockStructure = NULL, typeMatrixKT = "all")
{
  d = length( dataMatrix[1,] )
  n = length( dataMatrix[,1] )


  if(typeMatrixKT == "all")
  {
    estimate = matrix(data = 1 , nrow = d , ncol = d)
    for(j1 in 2:d)
    {
      for(j2 in 1:(j1-1))
      {
        vectorX1 = dataMatrix[ , j1]
        vectorX2 = dataMatrix[ , j2]
        estimate[j1,j2] = pcaPP::cor.fk(vectorX1 , vectorX2)
        estimate[j2,j1] = estimate[j1,j2]
      }
    }
  }


  if(typeMatrixKT == "aveDiag")
  {
    if(is.null(blockStructure))
    {
      stop(paste0("block structure not specified, when typeMatrixKT = ", typeMatrixKT))
    }

    if( sum(1:d %in% unlist(blockStructure)) == d & length(blockStructure) == d)
    {
      stop(paste0("block structure not a partition of 1:", d ))
    }

    totalGroups = length(blockStructure)
    estimate = matrix(data = 1 , nrow = totalGroups , ncol = totalGroups)
    for (g1 in 1:totalGroups)
    {
      for (g2 in 1:totalGroups)
      {
        if(g1 != g2)
        {
          diagSize = min(length(blockStructure[[g1]]),
                         length(blockStructure[[g2]]) )
          vectorBlockKT = rep(NA, times = diagSize)
          for (j in 1:diagSize)
          {
            vectorX1 = dataMatrix[ , blockStructure[[g1]][j] ]
            vectorX2 = dataMatrix[ , blockStructure[[g2]][j] ]
            vectorBlockKT[j] = pcaPP::cor.fk(vectorX1 , vectorX2)
          }
          blockKT = mean(vectorBlockKT)
          estimate[g1,g2] = blockKT
          estimate[g2,g1] = estimate[g1,g2]
        }
      }
    }
  }


  if(typeMatrixKT == "aveAll")
  {
    if(is.null(blockStructure))
    {
      stop(paste0("block structure not specified, when typeMatrixKT = ", typeMatrixKT))
    }

    if( sum(1:d %in% unlist(blockStructure)) == d & length(blockStructure) == d)
    {
      stop(paste0("block structure not a partition of 1:", d ))
    }

    totalGroups = length(blockStructure)
    print(totalGroups)
    estimate = matrix(data = 1 , nrow = totalGroups , ncol = totalGroups)
    for (g1 in 1:totalGroups)
    {
      for (g2 in 1:totalGroups)
      {
        if(g1 != g2)
        {
          matrixBlockKT = matrix(NA, nrow = length(blockStructure[[g1]]),
                                      ncol = length(blockStructure[[g2]]))
          for (j1 in 1:length(blockStructure[[g1]]) )
          {
            for (j2 in 1:length(blockStructure[[g2]]) )
            {
              vectorX1 = dataMatrix[ , blockStructure[[g1]][j1] ]
              vectorX2 = dataMatrix[ , blockStructure[[g2]][j2] ]
              matrixBlockKT[j1,j2] = pcaPP::cor.fk(vectorX1 , vectorX2)
            }
          }
          blockKT = mean(matrixBlockKT)
          estimate[g1,g2] = blockKT
          estimate[g2,g1] = estimate[g1,g2]
        }
      }
    }
  }


  if(typeMatrixKT == "aveRow")
  {
    if(is.null(blockStructure))
    {
      stop(paste0("block structure not specified, when typeMatrixKT = ", typeMatrixKT))
    }

    if( sum(1:d %in% unlist(blockStructure)) == d & length(blockStructure) == d)
    {
      stop(paste0("block structure not a partition of 1:", d ))
    }

    totalGroups = length(blockStructure)
    estimate = matrix(data = 1 , nrow = totalGroups , ncol = totalGroups)
    for (g1 in 1:totalGroups)
    {
      for (g2 in 1:totalGroups)
      {
        if(g1 != g2)
        {
          g1Size = length(blockStructure[[g1]])
          g2Size = length(blockStructure[[g2]])
          rowSize = min(g1Size, g2Size)
          gSmall = (g1Size <= g2Size) * g1 + (g1Size > g2Size) * g2
          gLarge = (g1Size <= g2Size) * g2 + (g1Size > g2Size) * g1
          vectorBlockKT = rep(NA, times = rowSize)
          for (j in 1:diagSize)
          {
            vectorX1 = dataMatrix[ , blockStructure[[gSmall]][j] ]
            vectorX2 = dataMatrix[ , blockStructure[[gLarge]][1] ]
            vectorBlockKT[j] = pcaPP::cor.fk(vectorX1 , vectorX2)
          }
          blockKT = mean(vectorBlockKT)
          estimate[g1,g2] = blockKT
          estimate[g2,g1] = estimate[g1,g2]
        }
      }
    }
  }


  return(estimate)
}






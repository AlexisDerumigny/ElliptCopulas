
test_that("KTMatrixEst outputs are well-formated", {
  matrixCor = matrix(c(1  , 0.5, 0.3 ,0.3,
                       0.5,   1, 0.3, 0.3,
                       0.3, 0.3,   1, 0.5,
                       0.3, 0.3, 0.5,   1), ncol = 4 , nrow = 4)

  dataMatrix = mvtnorm::rmvnorm(n = 20, mean = rep(0, times = 4),
                                sigma = matrixCor)

  blockStructure = list(1:2, 3:4)

  estKTMatrix = KTMatrixEst(dataMatrix = dataMatrix,
                            blockStructure = blockStructure,
                            averaging = "all")

  expect_identical(dim(estKTMatrix), c(2L,2L))
  expect_identical(diag(estKTMatrix), c(1,1))
  expect_identical(t(estKTMatrix), estKTMatrix)

  estKTMatrix = KTMatrixEst(dataMatrix = dataMatrix,
                            blockStructure = blockStructure,
                            averaging = "diag")

  expect_identical(dim(estKTMatrix), c(2L,2L))
  expect_identical(diag(estKTMatrix), c(1,1))
  expect_identical(t(estKTMatrix), estKTMatrix)

  estKTMatrix = KTMatrixEst(dataMatrix = dataMatrix,
                            blockStructure = blockStructure,
                            averaging = "row")

  expect_identical(dim(estKTMatrix), c(2L,2L))
  expect_identical(diag(estKTMatrix), c(1,1))
  expect_identical(t(estKTMatrix), estKTMatrix)

  estKTMatrix = KTMatrixEst(dataMatrix = dataMatrix,
                            averaging = "no")

  expect_identical(dim(estKTMatrix), c(4L,4L))
  expect_identical(diag(estKTMatrix), c(1,1,1,1))
  expect_identical(t(estKTMatrix), estKTMatrix)
})


test_that("KTMatrixEst reject inputs as intended", {
  matrixCor = matrix(c(1  , 0.5, 0.3 ,0.3,
                       0.5,   1, 0.3, 0.3,
                       0.3, 0.3,   1, 0.5,
                       0.3, 0.3, 0.5,   1), ncol = 4 , nrow = 4)

  dataMatrix = mvtnorm::rmvnorm(n = 20, mean = rep(0, times = 4),
                                sigma = matrixCor)

  expect_error( KTMatrixEst(dataMatrix = dataMatrix, # no block structure is provided
                            averaging = "all") )

  expect_error( KTMatrixEst(dataMatrix = dataMatrix,
                            blockStructure = list(c(1,2)), # only one block
                            averaging = "all") )

  expect_error( KTMatrixEst(dataMatrix = dataMatrix,
                            blockStructure = list(c(1,2,3,4)), # only one block
                            averaging = "all") )

  expect_error( KTMatrixEst(dataMatrix = dataMatrix,
                            blockStructure = list(c(1,2), c(1,3)), # variable 4 missing
                            averaging = "all") )

  expect_error( KTMatrixEst(dataMatrix = dataMatrix,
                            blockStructure = list(c(1,2), c(3)), # variable 4 missing
                            averaging = "all") )

  expect_error( KTMatrixEst(dataMatrix = dataMatrix,
                            blockStructure = list(c(1,2,4), c(3,4)), # variable 4 present in 2 blocks
                            averaging = "all") )
})


test_that("KTMatrixEst outputs have names if given", {
  matrixCor = matrix(c(1  , 0.5, 0.3 ,0.3,
                       0.5,   1, 0.3, 0.3,
                       0.3, 0.3,   1, 0.5,
                       0.3, 0.3, 0.5,   1), ncol = 4 , nrow = 4)

  dataMatrix = mvtnorm::rmvnorm(n = 20, mean = rep(0, times = 4),
                                sigma = matrixCor)

  estKTMatrix = KTMatrixEst(dataMatrix = dataMatrix,
                            blockStructure = list(1:2, 3:4),
                            averaging = "all")
  expect_null(colnames(estKTMatrix))
  expect_null(rownames(estKTMatrix))

  estKTMatrix = KTMatrixEst(dataMatrix = dataMatrix,
                            blockStructure = list(a = 1:2, b = 3:4),
                            averaging = "all")
  expect_identical(colnames(estKTMatrix), c("a", "b"))
  expect_identical(rownames(estKTMatrix), c("a", "b"))

  estKTMatrix = KTMatrixEst(dataMatrix = dataMatrix,
                            blockStructure = list(1:2, b = 3:4),
                            averaging = "diag")
  expect_identical(colnames(estKTMatrix), c("", "b"))
  expect_identical(rownames(estKTMatrix), c("", "b"))

})


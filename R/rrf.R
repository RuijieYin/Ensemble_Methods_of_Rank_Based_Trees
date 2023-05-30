#' Fit a Random Rank Forest model
#'
#' @param train_x A data frame or a matrix of predictors.
#' @param train_y A vector of response. If omitted, rrf will run in unsupervised mode.
#' @param n.tree Number of trees to grow.
#' @param mtry Number of gene pairs randomly sampled as input at each split.
#'
#' @return A rrf model
#' @export
#'
#' @examples
#' library(EnsembleRankTrees)
#' # load the example data:
#' data(liver_data)
#' model <- rrf(liver_data$train_x_liver, liver_data$train_y_liver)

rrf <- function (train_x, train_y = NULL, n.tree = 500, mtry = NULL) {

  # function: check if Check if a dataset has been row standardized (sample-wise)
  is_row_standardized <- function(data) {
    row_means <- rowMeans(data)
    row_vars <- apply(data, 1, var)
    all(abs(row_means) < 1e-10) && all(abs(row_vars - 1) < 1e-10)
  }

  # function: convert.genepairs
  # convert the raw expression values into gene pairs
  convert.genepairs <- function(input_dat = train_x) {
    # note: columns should be genes, rows are samples
    # the last column are the class labels
    # get gene pair data:
    input_expr = input_dat
    # get the number of columns
    n_col = dim(input_expr)[2]
    converted_input_expr = NULL
    for (i in 1:(n_col - 1)) {
      for (j in (i + 1):n_col) {
        converted_input_expr = cbind(converted_input_expr, as.numeric(input_expr[, i] < input_expr[, j]))
      }
    }
    data.list <- list("data.converted" = converted_input_expr)
    return(data.list)
  }


  if (!is.data.frame(train_x) && !is.matrix(train_x)) {
    # sanity check 1: if the input training data is a data frame or a matrix:
    stop("Please convert the input dataset to a data.frame or matrix object")
  } else if (is_row_standardized(train_x)) {
    # sanity check 2: if the input training data has been row standardized:
    stop("The dataset has been row standardized. Please do not row standardize the data.")
  } else {
    # remove symbols and special characters in column names, and replace them with "_"
    gsub("[[:punct:]]", "_", colnames(train_x))
    # convert input_dat to gene pairs:
    dat = convert.genepairs(train_x)$data.converted
    # merge expr and labels:
    data_train = cbind.data.frame(dat, as.factor(train_y))
    colnames(data_train)[dim(data_train)[2]] <- "subtype"
    rrf.model = rfsrc(subtype ~ .,
                  data = data_train,
                  ntree = n.tree,
                  mtry = mtry
    )
  }
  model.list <- list("model" = rrf.model)
  return(rrf.model)
}
















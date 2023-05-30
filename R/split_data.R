#' Split a dataset into train and validation data, to be used in
#' extract_rules and select_rules
#'
#'
#' @param train_x A data frame or a matrix of predictors.
#' @param train_y A vector of response.
#' @param seed_num The seed number to be used for random sampling process.
#' Setting a fixed seed number ensures reproducibility of the results.
#' @param valid_perc The percentage of the original training data that
#' can be used a validation data. By default, 30% samples in the original
#' training data is used as validation data.
#'
#' @return A training dataset and a validation dataset (the last column are the
#' class labels)
#' @export
#'
#' @examples
#' library(EnsembleRankTrees)
#' # load the example data:
#' data(liver_data)
#' # get training data and validation data
#' train_data <- split.train(train_x_liver, train_y_liver)$train.data
#' validation_data <- split.train(train_x_liver, train_y_liver)$validation.data

split.train <- function(train_x, train_y, seed_num = 1, valid_perc = 0.3) {

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

  # convert train_x to gene pairs:
  dat = convert.genepairs(train_x)$data.converted
  # merge train_x and class labels:
  data_input = cbind.data.frame(dat, as.factor(train_y))
  # rename the last column as phenotype:
  colnames(data_input)[dim(data_input)[2]] = "phenotype"

  # extract the phenotypes
  phenotypes = unique(train_y)

  # x%: training data (x% training; 1-x% validation)
  data_class_1 <- data_input[data_input$phenotype == phenotypes[1],]
  set.seed(seed_num)
  data_class_1.training <- data_class_1[sample(nrow(data_class_1), nrow(data_class_1)-round(nrow(data_class_1)*valid_perc)), ]
  data_class_1.validation <- data_class_1[!rownames(data_class_1) %chin% rownames(data_class_1.training),]

  data_class_2 <- data_input[data_input$phenotype == phenotypes[2],]
  set.seed(seed_num)
  data_class_2.training <- data_class_2[sample(nrow(data_class_2), nrow(data_class_2)-round(nrow(data_class_2)*valid_perc)), ]
  data_class_2.validation <- data_class_2[!rownames(data_class_2) %chin% rownames(data_class_2.training),]

  training_data <- rbind(data_class_1.training, data_class_2.training)
  validation_data <- rbind(data_class_1.validation, data_class_2.validation)

  # class type conversion
  # re-factorize
  names(training_data) <- make.names(names(training_data))
  training_data$phenotype <- factor(training_data$phenotype)
  names(validation_data) <- make.names(names(validation_data))
  validation_data$phenotype <- factor(validation_data$phenotype)




  data.split <- list("train.data" = training_data,
                     "validation.data" = validation_data)

  return(data.split)
}


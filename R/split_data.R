#' Split a dataset into train and validation data, to be used in \n
#' extract_rules and select_rules
#'
#'
#' @param train_x A data frame or a matrix of predictors.
#' @param train_y A vector of response.
#' @param seed_num The seed number to be used for random sampling process. \n
#' Setting a fixed seed number ensures reproducibility of the results.
#' @param valid_perc The percentage of the original training data that \n
#' can be used a validation data. By default, 30% samples in the original \n
#' training data is used as validation data.
#'
#' @return A training dataset and a validation dataset (the last column are the \n
#' class labels)
#' @export
#'
#' @examples
#' library(EnsembleRankTrees)
#' # load the example data:
#' data(train_x_liver, train_y_liver)
#' # get training data and validation data
#' train_data <- split.train(train_x_liver, train_y_liver)$train.data
#' validation_data <- split.train(train_x_liver, train_y_liver)$validation.data

split.train <- function(train_x, train_y, seed_num = 1, valid_perc = 0.3) {
  # extract the phenotypes
  phenotypes = unique(train_y)

  # merge:
  data_input = cbind.data.frame(train_x, as.factor(train_y))

  # x%: training data (1x% training; 1-x% validation)
  data_class_1 <- data_input[data_input$phenotype == phenotypes[1],]
  set.seed(seed_num)
  data_class_1.train.validation <- data_class_1[sample(nrow(data_class_1), round(nrow(data_class_1)*valid_perc)), ]
  set.seed(seed_num)
  data_class_1.training <- data_class_1[sample(nrow(data_class_1.train.validation), nrow(data_class_1)-round(nrow(data_class_1)*valid_perc)), ]
  data_class_1.validation <- data_class_1[!rownames(data_class_1.train.validation) %chin% rownames(data_class_1.training),]

  data_class_2 <- data_input[data$phenotype == phenotypes[2],]
  set.seed(seed_num)
  data_class_2.train.validation <- data_class_2[sample(nrow(data_class_2), round(nrow(data_class_2)*valid_perc)), ]
  set.seed(seed_num)
  data_class_2.training <- data_class_2[sample(nrow(data_class_2.train.validation), nrow(data_class_2)-round(nrow(data_class_2)*valid_perc)), ]
  data_class_2.validation <- data_class_2[!rownames(data_class_2.train.validation) %chin% rownames(data_class_2.training),]

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


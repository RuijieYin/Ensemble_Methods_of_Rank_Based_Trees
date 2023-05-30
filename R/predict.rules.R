#' Prediction function for rule-based methods
#'
#' @param data Test dataset
#' @param decision.rules The selected decision rules from select.rules
#'
#' @return
#' @export
#'
#' @examples
#' predict.rules(liver_data$test_x_liver, selection)$Prediction.results


predict.rules <- function(data, decision.rules) {
  for (i in 1:nrow(decision.rules)) {
    if (i > nrow(decision.rules)) {
      break
    } else {
      pred_class <- matrix(NA, nrow = dim(data)[1], ncol = 1)

      for (k in 1:dim(data)[1]) {
        # define x
        x <- data[k,]
        inter.table <- data.frame("label" = as.character(),
                                  "Performance Score" = as.numeric(),
                                  "Indicator" = as.numeric())

        for (w in 1:nrow(decision.rules))  {
          if (eval(parse(text = decision.rules[w,]$Rules)) == T) {
            inter.table[w,] <- list("label" = decision.rules[w,][,2],
                                    "Performance Score" = decision.rules[w,][,4],
                                    "Indicator" = 1)
            # indicator = 1: the sample's label is determined by 'class label'
          } else if (eval(parse(text = decision.rules[w,]$Rules)) == F) {
            inter.table[w,] <- list("label" = decision.rules[w,][,3],
                                    "Performance Score" = decision.rules[w,][,4],
                                    "Indicator" = 0)
            # indicator = 0: the sample's label is determined by 'else class'
          }
        }

        # determine the class label of the kth sample in validation data:
        label <- table(as.character(inter.table[,1]))
        if (length(label) != 1) {
          pred_class[k,1] <- as.character(inter.table[,1][which.max(inter.table[,3])])
        } else{
          pred_class[k,1] <- names(label[label == max(label)])
        }
      }
    }
  }

  prediction.results <- data.frame("sample no." = 1:nrow(data),
                                   "predicted class labels" = pred_class[,1])

  # # evaluation
  # # see performance: Cohen's Kappa
  # # 1st row: reference labels
  # # 2nd row: predicted labels
  # prediction <- matrix(prediction.results[,2], nrow = 1)
  #
  # compare_labels <- rbind(matrix(data[,dim(data)[2]], nrow = 1),
  #                         prediction)
  #
  # remove <- which(prediction[1,] == "UNS")
  #
  # if (length(remove) == 0) {
  #   compare_labels <- compare_labels
  # } else if (length(remove) != 0) {
  #   compare_labels <- compare_labels[,-remove]
  # }
  #
  # # Concordance: kappa statistic
  # # will need a matrix where the diagonal elements of the matrix are
  # # the agreeing elements; the discordant observations are on the off-diagonal.
  # # A confusion matrix:
  # con_table <- confusionMatrix(data = as.factor(compare_labels[2,]),
  #                              reference = as.factor(compare_labels[1,]))

  prediction.labels.kappa <- list("Prediction.results" = prediction.results)
                                  # "Results.table" = con_table)
  return(prediction.labels.kappa)
}

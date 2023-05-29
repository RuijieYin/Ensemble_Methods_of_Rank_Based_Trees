#' Select a set of rules that can be used to explain the model
#'
#' @param validation.data A separate validation data aimed to help select rules
#' from all decision rules obtained. Can be obtained from function split_data.
#' @param rule.table The object obtained from function extract_rule.
#'
#' @return A list of selected rules.
#' @return A plot of increase in accuracy on the validation dataset when adding rules to the ensemble.
#' @export
#'
#' @examples
#' selected_rules <- select.rule(validation_data, rules_table)$Selected.Rules
#' # show the line plot:
#' select.rule(validation_data, rules_table)$Plot.Num.Accuracy


source("R/treeRule_Obtain_Judge.R")

select.rules <- function(validation.data = validation_data, rule.table = rules_table) {

  # in validation data: rows are samples and columns are vars

  # remove rules with performance score = 0
  rule.table <- rule.table[rule.table$`Performance Score` != 0,]

  # the probability of predicting correctly by random guessing is 0.5*0.5+0.5*0.5 = 0.5
  rules.order.temp <- rule.table[rule.table$`Performance Score` > 0.5,]

  # order the rules by their performance scores:
  rules.order <- rules.order.temp[order(rules.order.temp$`Performance Score`, decreasing = T),]

  # initialize a data frame to store # of rules v.s. accuracy
  accuracy.table <- data.frame("Num.of.Rules" = as.numeric(),
                               "Accuracy" = as.numeric())

  # initialize accuracy
  accuracy <- 0

  # initialize rules will be retained:
  rules.retained.new <- data.frame("Rules" = character(),
                                   "Class Label" = character(),
                                   "Else Class" = character(),
                                   "Performance Score" = numeric(),
                                   "Performance Score (Class)" = numeric(),
                                   "Performance Score (Else Class)" = numeric())

  # initialize rules retained in the previous selection step:
  rules.retained.previous <- data.frame("Rules" = character(),
                                        "Class Label" = character(),
                                        "Else Class" = character(),
                                        "Performance Score" = numeric(),
                                        "Performance Score (Class)" = numeric(),
                                        "Performance Score (Else Class)" = numeric())

  for (i in 1:nrow(rules.order)) {

    if (i > nrow(rules.order)) {
      break
    } else {
      rules.add <- rules.order[i,]
      rules.retained.new <- rbind(rules.retained.previous, rules.add)

      pred_class <- matrix(NA, nrow = dim(validation.data)[1], ncol = 1)

      for (k in 1:dim(validation.data)[1]) {
        # define x
        x <- validation.data[k,]
        inter.table <- data.frame("label" = as.character(),
                                  "Performance Score" = as.numeric(),
                                  "Indicator" = as.numeric())

        for (w in 1:nrow(rules.retained.new))  {
          if (eval(parse(text = rules.retained.new[w,]$Rules)) == T) {
            inter.table[w,] <- list("label" = rules.retained.new[w,][,2],
                                    "Performance Score" = rules.retained.new[w,][,4],
                                    "Indicator" = 1)
            # indicator = 1: the sample's label is determined by 'class label'
          } else if (eval(parse(text = rules.retained.new[w,]$Rules)) == F) {
            inter.table[w,] <- list("label" = rules.retained.new[w,][,3],
                                    "Performance Score" = rules.retained.new[w,][,4],
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


      comp_table <- cbind(pred_class, as.matrix(validation.data$subtype, ncol = 1))

      accuracy.temp <- sum(comp_table[,1] == comp_table[,2])/nrow(comp_table)

      if (accuracy.temp <= accuracy) {
        accuracy <- accuracy
        rules.retained.previous <- rules.retained.previous
        accuracy.table[i,] <- list("Num.of.Rules" = i,
                                   "Accuracy" = accuracy)
      } else if (accuracy.temp > accuracy) {
        accuracy <- accuracy.temp
        accuracy.table[i,] <- list("Num.of.Rules" = i,
                                   "Accuracy" = accuracy)
        rules.retained.previous <- rules.retained.new
      }

    }
  }

  df <- data.frame(x = accuracy.table[,1], y = accuracy.table[,2])
  plot.list <- ggplot(data = df, aes(x = x, y = y, group=1)) +
    geom_line() +
    geom_point() +
    xlab("Decision rules") +
    ylab("Accuracy") +
    theme_classic()


  # # output the final selected rules:
  final.decision.rules <- rules.retained.previous

  # final number of rules:
  final_num <- nrow(final.decision.rules)

  optimal.rules <- list("Plot.Num.Accuracy" = plot.list,
                        "Selected.Rules" = final.decision.rules,
                        "Final.Num" = final_num)
  return(optimal.rules)

}

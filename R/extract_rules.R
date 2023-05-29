#' Extract rules containing comparing gene pairs from Random Rank Forest
#'
#' @param train_x A data frame or a matrix of predictors.
#' @param train_y A vector of response.
#' @param n.trees Number of trees to grow.
#' @param mtry Number of gene pairs randomly sampled as input at each split.
#' @param node_depth_max Maximum depth to which a tree should be grown. The default behavior is that this parameter is ignored.
#'
#' @return A list of rules extracted from the fitted Random Rank Forest model
#' @export
#'
#' @examples
#' library(EnsembleRankTrees)
#' # load the example data:
#' data(train_x_liver, train_y_liver)
#' rules_table <- extract.rules(train_x_liver, train_y_liver)

extract.rules <- function (train_x, train_y, n.trees = 500, mtry = 10, node_depth_max = 3) {
  # grow n.tree of trees
  set.seed(1)
  data_train = cbind.data.frame(train_x, as.factor(train_y))
  tree.collection <- rfsrc(subtype~., data = data_train,
                           ntree = n.trees, mtry = mtry, nodedepth = node_depth_max,
                           bootstrap = "by.root", samptype = 'swr',
                           membership = T,
                           # grow class balanced trees
                           # case.wt = randomForestSRC:::make.wt(data$subtype),
                           sampsize = randomForestSRC:::make.size(data$subtype))

  # summary(tree.collection)

  # calculate the total number of rules obtained from all trees
  total_num_rules <- length(getTreeRule(tree.collection)$tree.id)


  # initialize a data frame to store (rules[[j]], class, else_class, as.numeric(perform.score))
  final_result <- data.frame("Rules" = character(),
                             "Class Label" = character(),
                             "Else Class" = character(),
                             "Performance Score" = numeric(),
                             "Performance Score.1" = numeric(),
                             "Performance Score.2" = numeric())

  # extract tree rules from each tree:
  rules <- vector("list", total_num_rules)

  for (j in 1:total_num_rules) {

    # extract rule j from tree i:
    rules[[j]] <- parseRule(getTreeRule(tree.collection)$treeRule[[j]])
    # the ith tree where rule j comes from:
    tree.i <- getTreeRule(tree.collection)$tree.id[j]

    # index of inbag sample to grow tree i:
    index.inbag <- which(tree.collection$inbag[,tree.i] != 0)
    # find inbag samples suffice the jth rule:
    x <- data # should tell what x is
    all.sample.fit.rule <- which(eval(parse(text = rules[[j]])))
    inbag.sample.fit.rule <- Reduce(intersect, list(index.inbag, all.sample.fit.rule))

    # determine the class label for this rule
    # should consider the tied class:
    class_label_train <- table(data[inbag.sample.fit.rule, dim(data)[2]])
    class <- names(class_label_train[class_label_train == max(class_label_train)])

    # see if this rule has an "else": the majority class that does not suffice jth rule in the inbag sample,
    # use majority vote to determine "else":
    inbag.sample.not.fit.rule <- index.inbag[-inbag.sample.fit.rule]
    # find the class labels of samples that does not suffice jth rule;
    # first find and remove the majority class in 'class'label' section:
    majority_class_index <- which(as.character(data[inbag.sample.not.fit.rule, dim(data)[2]],class) %in% class)
    else_label_train <- table(as.character(data[inbag.sample.not.fit.rule, dim(data)[2]])[-majority_class_index])
    else_class_temp <- names(else_label_train[else_label_train == max(else_label_train)])
    # should consider tied votes:
    if (length(else_class_temp) != 1) {
      else_class <- NA
    } else {
      else_class <- else_class_temp
    }

    # calculate the performance score for rule j in tree i using oob data from tree i:
    # get class balanced OOB sample:
    out.of.bag <- which(tree.collection$inbag[,tree.i] == 0)
    oob.sample <- data[out.of.bag,]
    oob.sample.balanced <- oob.sample[complete.cases(oob.sample), ]
    #oob.sample.balanced <- downSample(x = data[,-dim(data)[2]],
    # y = factor(data[,dim(data)[2]]))
    colnames(oob.sample.balanced)[dim(oob.sample.balanced)[2]] <- "subtype"

    # 2 scenarios: the rule has or does not have an 'else'

    # 1st scenario: the rule j does not have an 'else':
    if (is.na(else_class) == T) {
      # store results in the iteration:
      class_label <- matrix(NA, nrow = dim(oob.sample.balanced)[1],
                            ncol = 2)
      # 1st element: see if sample suffice the above criteria, 1/0
      # 2nd element: see if sample is indeed the predicted class, 1/0

      for (k in 1:dim(oob.sample.balanced)[1]) {
        if (is.na(eval(parse(text = rules[[j]]))[out.of.bag[k]]) == F &&
            eval(parse(text = rules[[j]]))[out.of.bag[k]] == T &&
            oob.sample.balanced[k,]$subtype == class) {
          class_label[k,] <- c(1,1)
        } else if (is.na(eval(parse(text = rules[[j]]))[out.of.bag[k]]) == F &&
                   eval(parse(text = rules[[j]]))[out.of.bag[k]] == T &&
                   oob.sample.balanced[k,]$subtype != class) {
          class_label[k,] <- c(1,0)
        } else {
          class_label[k,] <- c(0,0)
        }
      }

      # see notes for details on how to calculate performance score
      perform.score.1 <- (sum(class_label[,1])/dim(class_label)[1])*
        (sum(class_label[class_label[,1]==1,2])/sum(class_label[,1]))

      perform.score.2 <- 0

      perform.score <- perform.score.1 + perform.score.2
    }
    # 2nd scenario: the rule j has an 'else':
    else if (is.na(else_class) == F) {
      # part 1:
      # store results in the iteration:
      class_label <- matrix(NA, nrow = dim(oob.sample.balanced)[1],
                            ncol = 2)
      # 1st element: see if sample suffice the above criteria, 1/0
      # 2nd element: see if sample is indeed the predicted class, 1/0

      for (k in 1:dim(oob.sample.balanced)[1]) {
        if (is.na(eval(parse(text = rules[[j]]))[out.of.bag[k]]) == F &&
            eval(parse(text = rules[[j]]))[out.of.bag[k]] == T &&
            oob.sample.balanced[k,]$subtype == class) {
          class_label[k,] <- c(1,1)
        } else if (is.na(eval(parse(text = rules[[j]]))[out.of.bag[k]]) == F &&
                   eval(parse(text = rules[[j]]))[out.of.bag[k]] == T &&
                   oob.sample.balanced[k,]$subtype != class) {
          class_label[k,] <- c(1,0)
        } else {
          class_label[k,] <- c(0,0)
        }
      }

      # see notes for details on how to calculate performance score
      perform.score.1 <- (sum(class_label[,1])/dim(class_label)[1])*
        (sum(class_label[class_label[,1]==1,2])/sum(class_label[,1]))

      # part 2:
      # store results in the iteration:
      class_label <- matrix(NA, nrow = dim(oob.sample.balanced)[1],
                            ncol = 2)
      # 1st element: see if sample suffice the above criteria, 1/0
      # 2nd element: see if sample is indeed the predicted class, 1/0

      for (k in 1:dim(oob.sample.balanced)[1]) {
        if (is.na(eval(parse(text = rules[[j]]))[out.of.bag[k]]) == F &&
            eval(parse(text = rules[[j]]))[out.of.bag[k]] == F &&
            oob.sample.balanced[k,]$subtype == else_class) {
          class_label[k,] <- c(1,1)
        } else if (is.na(eval(parse(text = rules[[j]]))[out.of.bag[k]]) == F &&
                   eval(parse(text = rules[[j]]))[out.of.bag[k]] == F &&
                   oob.sample.balanced[k,]$subtype != else_class) {
          class_label[k,] <- c(1,0)
        } else {
          class_label[k,] <- c(0,0)
        }
      }

      # see notes for details on how to calculate performance score
      perform.score.2 <- (sum(class_label[,1])/dim(class_label)[1])*
        (sum(class_label[class_label[,1]==1,2])/sum(class_label[,1]))

      perform.score <- perform.score.1 + perform.score.2
    }

    final_result[j,] <- list("Rules" = rules[[j]],
                             "Class Label" = class,
                             "Else Class" = else_class,
                             "Performance Score" = perform.score,
                             "Performance Score.1" = perform.score.1,
                             "Performance Score.2" = perform.score.2)

  }
  # remove rules with NaN performance scores
  final_result_temp_1 <- final_result[!is.na(final_result[,4]),]

  # remove rules with performance scores = 0
  final_result_temp_2 <- final_result_temp_1[final_result_temp_1[,4] != 0,]

  # remove duplicated rules
  final_result_all <- final_result_temp_2[!duplicated(final_result_temp_2),]

  # remove duplicated rules by column (attempt 2)
  final_result_all <- final_result_all %>% distinct(final_result_all[,1], .keep_all = TRUE)

  colnames(final_result_all) <- c("Rules", "Class Label", "Else Class", "Performance Score",
                                  "Performance Score (Class)", "Performance Score (Else Class)")
  rules.all <- list("Rules" = final_result_all)
  # function ends here
}

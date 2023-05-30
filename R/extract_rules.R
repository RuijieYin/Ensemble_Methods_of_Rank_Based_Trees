#' Extract rules containing comparing gene pairs from Random Rank Forest
#'
#' @param train_x A data frame or a matrix of predictors.
#' @param train_y A vector of response.
#' @param n.trees Number of trees to grow.
#' @param mtry Number of gene pairs randomly sampled as input at each split.
#' @param node_depth_max Maximum depth to which a tree should be grown.
#' The default behavior is that this parameter is ignored.
#'
#' @return A list of rules extracted from the fitted Random Rank Forest model
#' @export
#'
#' @examples
#' library(EnsembleRankTrees)
#' # load the example data:
#' data(train_x_liver, train_y_liver)
#' # get training data and validation data
#' train_data <- split.train(train_x_liver, train_y_liver)$train.data
#' validation_data <- split.train(train_x_liver, train_y_liver)$validation.data
#' extract rules:
#' rules_table <- extract.rules(train_x_liver, train_y_liver)

# source("R/treeRule_Obtain_Judge.R")

library(data.table)
library(data.tree)
library(DiagrammeR)
library(dplyr)
library(ggplot2)
library(parallel)
library(randomForestSRC)


extract.rules <- function (train_x, train_y, n.trees = 500, mtry = 10, node_depth_max = 3) {
  ## select rules randomly from trees, rule index must be sequential
  rfolds <- function (max.rules.tree, lfc) {

    ntree <- length(lfc)
    max.rules <- max.rules.tree * ntree
    tree <- rep(1:ntree, lfc)
    lfc <- unlist(sapply(lfc, function(lc){1:lc}))
    idx <- sample(1:length(lfc), size = min(max.rules, length(lfc)))
    tree <- sort(tree[idx])
    lfc <- unlist(tapply(tree, tree, function(z) {1:length(z)}))
    cbind(tree, lfc)

  }

  getTreeRule.short <- function(object, tree.id = b){

    tolerance = sqrt(.Machine$double.eps)

    ## pull xvar.names
    xvar.names <- object$xvar.names
    xvar.factor <- object$xvar.factor

    ## pull the arrays
    native.array <- object$native.array
    native.f.array <- object$native.f.array

    ## may add processing needed for factors
    f.ctr <- 0
    factor.flag <- FALSE


    ## define the display tree
    display.tree <- native.array[native.array$treeID == tree.id,, drop = FALSE]

    converted.tree <- display.tree
    vars.id <- data.frame(var = c("<leaf>", xvar.names), parmID = 0:length(xvar.names), stringsAsFactors = FALSE)
    converted.tree$var <- vars.id$var[match(display.tree$parmID, vars.id$parmID)]

    ## special symbol to be used for encoding the counter for variables (see note below)
    special <- "999_999"

    var.count <- 1:nrow(converted.tree)
    lapply(unique(converted.tree$var), function(vv) {
      pt <- converted.tree$var == vv
      var.count[which(pt)] <<- 1:sum(pt)
    })

    converted.tree$var_count <- var.count
    converted.tree$var_conc <- paste0(converted.tree$var, special, converted.tree$var_count)

    ## preliminary
    from_node <- ""
    network <- data.frame()
    num.children <- data.frame(converted.tree, children = 0)
    num.children <- num.children[num.children$var != "<leaf>",, drop = FALSE]
    num.children <- num.children[!duplicated(num.children$var_conc),, drop = FALSE]
    num_children <- as.list(rep(0, nrow(num.children)))
    names(num_children) <- num.children$var_conc


    ## loop (using lapply)
    lapply(1:nrow(converted.tree), function(i) {
      rowi <- converted.tree[i, ]
      xs <- converted.tree$contPT[converted.tree$var_conc == from_node]
      if (i == 1){
        from_node <<- rowi$var_conc
      }
      else{
        ## develop the split encoding
        if (num_children[[from_node]] == 0) {#left split
          side <- "<="
          contPT.pretty <- round(as.numeric(xs, 3))
          split_ineq_pretty <- paste0(side, contPT.pretty)
        }
        else {#right split
          side <- ">"
          split_ineq_pretty <- ""
        }

        if (is.numeric(xs)) {
          xs <- xs + tolerance
        }
        split_ineq <- paste0(side, xs)

        ## update the network
        to_node <- rowi$var_conc
        new_node <- list(from = from_node, to = to_node, split = split_ineq, split.pretty = split_ineq_pretty)
        network <<- data.frame(rbind(network, new_node, stringsAsFactors = FALSE))
        num_children[[from_node]] <<- num_children[[from_node]] + 1
        if(rowi$var != "<leaf>")
          from_node <<- to_node
        else{
          if(i != nrow(converted.tree)){
            while(num_children[[from_node]] == 2){
              from_node <<- network$from[network$to == from_node]
            }
          }
        }
      }
    })



    data.tree.network <- data.tree::FromDataFrameNetwork(network, "split")

    ctr <- 0
    treerule <- varUsed <- varNames <- list()

    lapply(data.tree.network$leaves, function(node) {


      ## pull relevant information
      path_list <- node$path
      var_list <- sapply(path_list, function(x){strsplit(x, special)[[1]][1]})
      var_list[length(var_list)] <- ""
      node_iter <- data.tree.network

      ## make boolean string operator - save the list of variable names
      varnames <- NULL
      call <- lapply(2:(length(path_list)), function(i) {
        node_iter <<- node_iter[[path_list[[i]]]]
        str <- node_iter$split ###XXXXXXXXXXXXXXXXXX change to node_iter$split.pretty for pretty digits. Since we use quantiles (range 0~100), I think keep node_iter$split.pretty as integers may OK
        varnames <<- c(varnames, var_list[i-1])
        ## numeric boolean operator
        if (!any(grepl("\\{", str))) {
          str <- paste0("", paste0(var_list[i-1], str))
        }
        ## complementary pair boolean operator
        else {
          str <- gsub("\\{", "", str)
          str <- gsub("\\}", "", str)
          str <- strsplit(str, ",")[[1]]
          str <- paste("==", str, sep = "")
          str <- paste0("(",paste(paste0("", var_list[i-1], str), collapse = "|"),")")
        }
        str
      })
      names(varnames) <- NULL

      ## update the counter and save the results
      ctr <<- ctr + 1
      treerule[[ctr]] <<- call
      varUsed[[ctr]] <<- sort(unique(varnames))
      varNames[[ctr]] <<- varnames

    })

    list(treeRule = treerule, varUsed = varUsed, varNames = varNames)

  }


  getTreeRule <- function(object){
    ntree <- object$ntree

    lfc <- object$leaf.count[1:ntree]
    treeRuleSeq <- rfolds(ntree, lfc)


    xvar.names <- object$forest$xvar.names
    xvar.factor <- object$forest$xvar.factor
    p <- length(xvar.names)

    ## extract the data and process it
    ## missing data not allowed
    ## convert the data to numeric mode, apply the na.action protocol
    xvar <- object$forest$xvar
    xvar <- randomForestSRC:::finalizeData(xvar.names, xvar, miss.flag = FALSE)

    arr <- object$forest$nativeArray
    arrf <- object$forest$nativeFactorArray[[1]]
    pt.arr <- is.element(paste0(arr$treeID, ":", arr$nodeID),
                         paste0(treeRuleSeq[, 1], ":", treeRuleSeq[, 2]))

    ptf.arr <- arr$mwcpSZ != 0
    arr <- arr[pt.arr,, drop = FALSE]

    ntreeSeq <- sort(unique(arr$treeID))

    ## now reduce the object to the minimal information
    object <- list(xvar.names = xvar.names,
                   xvar.factor = xvar.factor,
                   native.array = arr)

    treeRuleO <- mclapply(ntreeSeq, function(b) {
      getTreeRule.short(object, tree.id = b)
    })

    treeRule <- unlist(lapply(treeRuleO, function(oo) {oo$treeRule}), recursive = FALSE)
    if (!is.null(treeRule)) {## convert list of lists to a list
      treeRule <- lapply(treeRule, function(oo) {unlist(oo)})
    }

    varUsed <- unlist(lapply(treeRuleO, function(oo) {oo$varUsed}), recursive = FALSE)
    varNames <- unlist(lapply(treeRuleO, function(oo) {oo$varNames}), recursive = FALSE)

    tree.id <- unlist(lapply(1:length(treeRuleO), function(j) {rep(j, length(treeRuleO[[j]]$treeRule))}))

    rm(treeRuleO)

    list(treeRule = treeRule,
         varUsed = varUsed,
         varNames = varNames,
         tree.id = tree.id,
         treeSeq = sort(unique(treeRuleSeq[, 1])))

  }
  parseRule <- function(rule) {

    anyC <- grepl("\\(", rule)

    if (sum(anyC) == 0) {
      paste0("x$", rule, collapse=" & ")
    }
    else {
      unlist(sapply(1:length(rule), function(j) {
        rj <- rule[j]
        if (anyC[j]) {
          rj <- sub("\\(", "", rj)
          rj <- sub("\\)", "", rj)
          rj <- strsplit(rj, "\\|")[[1]]
          unlist(lapply(rj, function(rr) {
            paste0("x$", rr)
          }))
        }
        else {
          paste0("x$", rj)
        }
      }))
    }

  }


  # grow n.tree of trees
  set.seed(1)
  data = cbind.data.frame(train_x, as.factor(train_y))
  tree.collection <- rfsrc(subtype~., data = data,
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

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


# source("R/treeRule_Obtain_Judge.R")

library(data.table)
library(data.tree)
library(DiagrammeR)
library(dplyr)
library(ggplot2)
library(parallel)
library(randomForestSRC)

select.rules <- function(validation.data = validation_data, rule.table = rules_table) {
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

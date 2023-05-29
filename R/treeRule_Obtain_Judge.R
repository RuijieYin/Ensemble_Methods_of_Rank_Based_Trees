library(randomForestSRC)
library(parallel)
library(data.tree)
library(DiagrammeR)
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


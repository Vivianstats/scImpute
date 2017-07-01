imputation_model1 <-
function (count, point, parslist, drop_thre = 0.5, method = 2, 
    ncores) 
{
    count = as.matrix(count)
    count_imp = count
    valid_genes = which( (rowSums(count) > point * ncol(count)) & 
                           complete.cases(parslist) )
    count = count[valid_genes, , drop = FALSE]
    parslist = parslist[valid_genes, , drop = FALSE]
    I = nrow(count)
    J = ncol(count)
    droprate = t(sapply(1:I, function(i) {
        wt = calculate_weight(count[i, ], parslist[i, ])
        return(wt[, 1])
    }))
    impute = function(I, cellid, count, parslist, droprate, method) {
        yobs = count[, cellid]
        yimpute = rep(0, I)
        geneid_drop = which(droprate[, cellid] >= drop_thre)
        if (length(geneid_drop) == 0 | length(geneid_drop) == I) {
            yimpute = yobs
        }
        else {
            geneid_obs = setdiff(1:I, geneid_drop)
            xx = count[geneid_obs, -cellid]
            yy = count[geneid_obs, cellid]
            weight = 1 - parslist[geneid_obs, "rate"]
            set.seed(cellid)
            model_cv = cv.glmnet(xx, yy, alpha = 1, weights = weight, 
                nfolds = 5)
            vars_nonzero = predict(model_cv, type = "nonzero", 
                s = "lambda.1se")
            if (class(vars_nonzero) == "list") {
                return(yobs)
            }
            else {
                vars_nonzero = vars_nonzero[, 1]
                xselect = xx[, vars_nonzero, drop = FALSE]
                model_ols = lm(yy ~ xselect, weights = weight)
                ximpute = count[geneid_drop, -cellid]
                ximpute = ximpute[, vars_nonzero, drop = FALSE]
                xnew = cbind(rep(1, length(geneid_drop)), ximpute)
                ynew = xnew %*% matrix(coef(model_ols), ncol = 1)
                if (method == 2) {
                  check = apply(ximpute, 1, function(x) {
                    sum(x == log10(1.01))
                  })
                  ynew[check == ncol(ximpute)] = point
                }
                yimpute[geneid_drop] = ynew
                yimpute[geneid_obs] = yobs[geneid_obs]
                return(yimpute)
            }
        }
    }
    res = mclapply(1:J, function(cellid) {
        if (cellid %% 50 == 0) 
            print(cellid)
        y = try(impute(I, cellid, count, parslist, droprate, 
            method), silent = TRUE)
        if (class(y) == "try-error") {
            print(y)
            y = count[, cellid]
        }
        return(y)
    }, mc.cores = ncores)
    res = Reduce(cbind, res)
    count_imp[valid_genes, ] = res
    count_imp[count_imp < point] = point
    return(count_imp)
}




imputation_model1_bytype <-
  function (count, labels, point, parslist, drop_thre = 0.5, method = 2, 
            ncores) 
  {
    impute = function(I, cellid, cellIDs, count, parslist, droprate, method) {
      yobs = count[, cellid]
      yimpute = rep(0, I)
      geneid_drop = which(droprate[, cellid] >= drop_thre)
      if (length(geneid_drop) == 0 | length(geneid_drop) == I) {
        yimpute = yobs
      }
      else {
        geneid_obs = setdiff(1:I, geneid_drop)
        xx = count[geneid_obs, setdiff(cellIDs, cellid)]
        yy = count[geneid_obs, cellid]
        weight = 1 - parslist[geneid_obs, "rate"]
        set.seed(cellid)
        model_cv = cv.glmnet(xx, yy, alpha = 1, weights = weight, 
                             nfolds = 5)
        vars_nonzero = predict(model_cv, type = "nonzero", 
                               s = "lambda.1se")
        if (class(vars_nonzero) == "list") {
          return(yobs)
        }
        else {
          vars_nonzero = vars_nonzero[, 1]
          xselect = xx[, vars_nonzero, drop = FALSE]
          model_ols = lm(yy ~ xselect, weights = weight)
          ximpute = count[geneid_drop, setdiff(cellIDs, cellid)]
          ximpute = ximpute[, vars_nonzero, drop = FALSE]
          xnew = cbind(rep(1, length(geneid_drop)), ximpute)
          ynew = xnew %*% matrix(coef(model_ols), ncol = 1)
          if (method == 2) {
            check = apply(ximpute, 1, function(x) {
              sum(x == log10(1.01))
            })
            ynew[check == ncol(ximpute)] = point
          }
          yimpute[geneid_drop] = ynew
          yimpute[geneid_obs] = yobs[geneid_obs]
          return(yimpute)
        }
      }
    }
    count = as.matrix(count)
    count_imp = count
    valid_genes = which( (rowSums(count) > point * ncol(count)) & 
                           complete.cases(parslist) )
    count = count[valid_genes, , drop = FALSE]
    parslist = parslist[valid_genes, , drop = FALSE]
    I = nrow(count)
    J = ncol(count)
    droprate = t(sapply(1:I, function(i) {
      wt = calculate_weight(count[i, ], parslist[i, ])
      return(wt[, 1])
    }))
    type_list = lapply(unique(labels), function(x) {
      return(which(labels == x))
    })
    for (kk in 1:length(type_list)) {
      cellIDs = type_list[[kk]]
      res = mclapply(cellIDs, function(cellid) {
        if (cellid %% 50 == 0) print(cellid)
        y = try(impute(I, cellid, cellIDs, count, parslist, droprate, 
                       method), silent = TRUE)
        if (class(y) == "try-error") {
          print(y)
          y = count[, cellid]
        }
        return(y)
      }, mc.cores = ncores)
      res = Reduce(cbind, res)
      count_imp[valid_genes, cellIDs] = res
    }
    count_imp[count_imp < point] = point
    return(count_imp)
  }

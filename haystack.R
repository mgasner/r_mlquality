library(veritable.r)
library(impute)
library(digest)

options(veritable.server = "192.168.0.185", veritable.username = "veritable")

source("~/r_utilities/memoize.R")

# 3 clusters, one at -1,-1, one at 0,0, and one at 1,1
# distractor dimensions are 0 mean, stdev .1

signal <- function (n) {
  mix <- runif(n)
  mix
}

signal.column <- function (d, mix, std, means = c(-0.25, 0, 0.25), offset = 0) {
  n <- length(mix)
  mat <- matrix(nrow = n, ncol = d)
  for (j in 1:d) {
    for (i in seq_along(mix)) {
      p = mix[i]
      frac = 1.0 / length(means)
      
      for (k in seq_along(means)) {
        if (p < frac) {
          mat[i,j] <- offset + rnorm(1, means[k], std)
          break
        } else {
          p <- p - frac
        }
      }
    }
  }  
  df <- as.data.frame(mat)
  df
}

distractor.column <- function (d, n, std, means = c(-10, -5, 0, -5, 10)) {
  mix <- runif(n)
  signal.column(d, mix, std, means)
  #mix <- rep(0.5, n)
  #signal.column(d, mix, std, c(0))
}

haystack <- function (d.sig, d.noise, n, std.sig, std.noise) {
  # assemble the signal data frame
  components <- signal(n)
  sig.df <- signal.column(d.sig, components, std.sig)
  
  # assemble the noise data frame
  if (d.noise > 0) {
    noise.df <- distractor.column(d.noise, n, std.noise)
    # paste them together and return
    df <- cbind(sig.df, noise.df)
  } else {
    df <- sig.df
  }
  
  df
}

omit <- function(df, p, cols) {
	nrows <- nrow(df)
  df <- as.matrix(df)
	for (i in 1:nrows) {
	  for (j in cols) {
	    if (runif(1) < (p / 100)) {
	      df[i, j] <- NA
	    }
	  }
	}
	df <- as.data.frame(df)
	df
}

score <- function(full.df, missing.df, imputed.df, cols, type = "abs") {
  err = 0
  count = 0
  
  N = dim(full.df)[1]
  D = dim(full.df)[2]
  
  #vars = vector(mode = "numeric", length = length(cols))
  #for (i in 1:length(cols)) {
  #  vars[i] = var(full.df[, cols[i]])
  #}

  #print(vars)  
  for (i in 1:N) {
    for (j in cols) {
      if (is.na(missing.df[i, j])) {
        if (type == "abs") {
          delta <- abs(full.df[i, j] - imputed.df[i, j])
          #print(paste(i,j,delta))
          #err <- err + (delta / vars[j])
          err <- err + delta
          count <- count + 1
        } else {
          stop ("Error in score, didn't recognize type:", type)
        }
      }
    }
  }
  err / count
}

gen.haystack.df <- function(seed, d.sig, d.max, N, omit.percent, std.sig, std.noise) {
  set.seed(seed)
  full.df <- haystack(d.sig, d.max - d.sig, N, std.sig, std.noise)
  missing.df <- omit(full.df, omit.percent, 1:d.sig)
  
  list(full = full.df, missing = missing.df)
}

mem.gen.haystack.df <- memoize(gen.haystack.df)

gen.imputation.val <- function(d.max, d.sig, d.noise, N, method, train_seed, inf_seed, omit.percent, std.sig, std.noise) {
  print(paste("GETTING DATA: ", d.max,  d.sig, d.noise, N, method, train_seed, inf_seed, omit.percent, std.sig, std.noise))

   res = mem.gen.haystack.df(train_seed, d.sig, d.max, N, omit.percent, std.sig, std.noise)
   full.df = res[["full"]][, 1:(d.sig+d.noise)]
   missing.df = res[["missing"]][, 1:(d.sig+d.noise)]
      
   set.seed(inf_seed)
      
   if (method == "hastie") {
    
    res <- impute.knn(t(missing.df), k = 1)
    imputed.df <- as.data.frame(t(res[["data"]]))
    .Random.seed <<- res[["rng.state"]]
    
   } else if (method == "avg") {
     colavgs <- c()
     for (col in names(full.df)) {
      colavgs <- c(colavgs, mean(full.df[,col]))
     }

     mat <- as.matrix(missing.df)

     for (i in 1:dim(full.df)[1]) {
       for (j in 1:dim(full.df)[2]) {
         if (is.na(mat[i, j])) {
          mat[i,j] <- colavgs[j]
        }
       }
     }
     imputed.df <- as.data.frame(mat)
     
   } else if (method == "random") {
     mat <- as.matrix(missing.df)
     full.mat <- as.matrix(full.df)
     for (i in 1:dim(full.df)[1]) {
       for (j in 1:dim(full.df)[2]) {
         if (is.na(mat[i, j])) {
           mat[i,j] <- full.mat[sample(1:dim(full.df)[1], 1), j]
         }
       }
     }
     imputed.df <- as.data.frame(mat)
   } else if (method == "veritable") {
    mat <- as.matrix(missing.df)
    ds <- as.veritable.dataset(missing.df)
    datatypes(ds) <- "continuous"
    h <- start.analysis(ds, samples = 10, iterations = 500, max_time = 3600)
    res <- wait.for.analysis(h, ds)
    plot(feature.zmatrix(res))
    ph <- list()
    for (i in 1:dim(missing.df)[1]) {
      if (length(which(is.na(missing.df[i,]))) != 0) {
        ph[[i]] <- start.predictions(h, ds, 1000, fixed.features = missing.df[i, names(missing.df)[-which(is.na(missing.df[i,]))]], predicted.features = names(missing.df)[which(is.na(missing.df[i,]))])
      } else {
        ph[[i]] <- "NONE"
      }
    }
    pr <- list()
    for (i in 1:dim(missing.df)[1]) {
      print(paste("GETTING PREDICTIONS FOR ROW", i))
      if (class(ph[[i]])[1] == "veritable.predictions.handle") {
        pr[[i]] <- wait.for.predictions(ph[[i]], ds)
      } else {
        pr[[i]] <- missing.df[i,]
      }
    }
    
    for (i in 1:dim(missing.df)[1]) {
      print(paste("REPLACING MISSING VALUES FOR ROW", i))
      if (class(pr[[i]])[1] == "veritable.predictions") {
        for (j in which(is.na(missing.df[i,]))) {
          mat[i,j] <- median(pr[[i]]@data[,j])
        }
      }
    }
    imputed.df <- as.data.frame(mat)
   } else {
    stop("unknown method: ", method)
   }
   
   #hist(abs(as.matrix(full.df - imputed.df)[,1:d.sig]))
   s = score(full.df, missing.df, imputed.df, 1:d.sig)   
   s
}

mem.gen.imputation.val <- memoize(gen.imputation.val)

if(exists(".Random.seed")) rm(.Random.seed)
set.seed(0)

#methods <- list(
#  list("hastie"
#)
methods <- c("hastie")
#methods <- c("hastie", "avg", "random")
#methods <- c("hastie", "avg", "random", "veritable")
colors <- list(hastie = "#FF0000FF", avg = "#00FF00FF", random = "#0000FFFF", veritable = "#8000FFFF")
d.sig <- 10
d.noise <- seq(from = 0, to = 500, by = 10)
#d.noise = 5
#d.noise <- c(0, 1, 5, 15)
d.max <- max(d.noise) + d.sig
std.sig <- 0.05
std.noise <- 0.05
#std.noise <- 100
N <- 100
#omit.percent <- seq(1, 20, by = 1)
omit.percent <- 5
#omit.percent <- c(10, 50, 70)

res = mem.gen.haystack.df(0, d.sig, d.max, N, 10, std.sig, std.noise)
full.df = res[["full"]][,]
missing.df = res[["missing"]][,]

pdf(file = paste("hastie_converge", std.noise, ".pdf", sep = ""), onefile = TRUE)

#lapply(omit.percent, function (pct) {
#results = lapply(methods, function(method) {
#            list(method, lapply(d.noise, function(d.noise) {
#              mem.gen.imputation.val(d.max, d.sig, d.noise, N, method, 0, 0, pct, std.sig, std.noise)
#            }))
#          })
#
#plot(d.noise, rep(NA, length(d.noise)), ylim=c(0, 1), main = paste(pct, "% omitted", sep = ""))
#lapply(results, function(res) {
#    # needs to be labeled with res[[1]]
#    print(paste(res[[1]], res[[2]]))
#    lines(d.noise, res[[2]], col = colors[[res[[1]]]])
#})
#})

results <- list()
for (i in 1:100) {
  result <- c()
  for (j in d.noise) {
    result <- c(result, mem.gen.imputation.val(d.max, d.sig, j, N, methods[[1]], 0, 0, omit.percent, std.sig, std.noise))
  }
  results[[i]] <- result
}
lapply(results, function (r) plot(d.noise, r, type ="l"))
dev.off()
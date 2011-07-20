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
  print(paste("GETTING DATA: ", d.max,  d.sig, d.noise, N, method[["name"]], train_seed, inf_seed, omit.percent, std.sig, std.noise))

   res = mem.gen.haystack.df(train_seed, d.sig, d.max, N, omit.percent, std.sig, std.noise)
   full.df = res[["full"]][, 1:(d.sig+d.noise)]
   missing.df = res[["missing"]][, 1:(d.sig+d.noise)]
      
   set.seed(inf_seed)
      
   if (method[["name"]] == "hastie") {
    
    res <- impute.knn(t(missing.df), k = method[["k"]])
    imputed.df <- as.data.frame(t(res[["data"]]))
    .Random.seed <<- res[["rng.state"]]
    
   } else if (method[["name"]] == "avg") {
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
     
   } else if (method[["name"]] == "random") {
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
   } else if (method[["name"]] == "veritable") {
    mat <- as.matrix(missing.df)
    ds <- as.veritable.dataset(missing.df)
    datatypes(ds) <- "continuous"
    h <- start.analysis(ds, samples = method[["samples"]], iterations = method[["iterations"]], max_time = method[["max_time"]])
    res <- wait.for.analysis(h, ds)
    plot(feature.zmatrix(res))
    ph <- list()
    for (i in 1:dim(missing.df)[1]) {
      if (length(which(is.na(missing.df[i,]))) != 0) {
        ph[[i]] <- start.predictions(h, ds, method[["predictions"]], fixed.features = missing.df[i, names(missing.df)[-which(is.na(missing.df[i,]))]], predicted.features = names(missing.df)[which(is.na(missing.df[i,]))])
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

methods <- list(
  list(name = "hastie", k = 1, color = "#FF0000FF"),
  list(name = "avg", color = "#00FF00FF"),
  list(name = "random", color = "#0000FFFF"),
  list(name = "veritable",
       samples = 10,
       iterations = 500,
       max_time = 3600,
       predictions = 1000,
       color = "#8000FFFF")
)

SEED <- 0

experiment <- list(
  name = "5pctomit",
  methods = methods,
  d.sig = 10,
  d.noise = 0:10,
  d.max = 20,
  std.sig = 0.05,
  std.noise = 0.05,
  N = 100,
  omit.percent = 5
)

run.experiment <- function (experiment) {
  pdf(file = paste(experiment[["name"]], ".pdf", sep = ""), onefile = TRUE)
  
  lapply(experiment[["omit.percent"]], function (pct) {
    results <- lapply(experiment[["methods"]], function (method) {
      lapply(experiment[["d.noise"]], function (d.noise) {
        mem.gen.imputation.val(experiment[["d.max"]], experiment[["d.sig"]], d.noise, experiment[["N"]], method, 0, 0, pct, experiment[["std.sig"]], experiment[["std.noise"]])
      })
    })
  })
  
  plot(d.noise, rep(NA, length(d.noise)), ylim=c(0, 1), main = paste(pct, "% omitted", sep = ""))
  lapply(results, function(res) {
    # needs to be labeled with res[[1]]
    print(paste(res[[1]][["name"]], res[[2]]))
    lines(d.noise, res[[2]], col = res[[1]][["color"]])
  })
  dev.off()
}
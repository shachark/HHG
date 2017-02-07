hhg.example.datagen = function(n, example) {
  if (example == '') {
  } else if (example == '4indclouds') {
    .datagen4indclouds(n)
  } else if (example == '2Parabolas') {
    .datagen2Parabolas(n)
  } else if (example == 'W') {
    .datagenW(n)
  } else if (example == 'Parabola') {
    .datagenParabola(n)
  } else if (example == 'Diamond') {
    .datagenDiamond(n)
  } else if (example == 'Circle') {
    .datagenCircle(n)
  } else if (example == 'TwoClassUniv') {
    .datagenTwoClassUniv(n)
  } else if (example == 'FourClassUniv') {
    .datagenFourClassUniv(n)
  } else if (example == 'TwoClassMultiv') {
    .datagenTwoClassMultiv(n)
  } else {
    stop('Unexpected example specified. Please consult the documentation.')
  }
}

.datagen4indclouds = function(n) {
  dx = rnorm(n) / 3
  dy = rnorm(n) / 3
  cx = sample(c(-1, 1), size = n, replace = T)
  cy = sample(c(-1, 1), size = n, replace = T)
  u = cx + dx
  v = cy + dy 
  return (rbind(u, v))
}

.datagen2Parabolas = function(n) {
  x = seq(-1, 1, length = n)
  y = (x ^ 2 + runif(n) / 2) * (sample(c(-1, 1), size = n, replace = T))
  return (rbind(x, y))
}

.datagenW = function(n) {
  x = seq(-1, 1, length = n)
  u = x + runif(n)/3
  v =  4*( ( x^2 - 1/2 )^2 + runif(n)/500 )
  return (rbind(u,v))
}

.datagenParabola = function(n) {
  x = seq(-1, 1, length = n)
  y = (x ^ 2 + runif(n)) / 2
  return (rbind(x,y))
}

.datagenDiamond = function(n) {
  x = runif(n, min = -1, max = 1)
  y = runif(n, min = -1, max = 1)

  theta = -pi / 4
  rr = rbind(c(cos(theta), -sin(theta)),
             c(sin(theta),  cos(theta)))
  tmp = cbind(x, y) %*% rr
  u = tmp[,1]
  v =  tmp[,2]
  return (rbind(u, v))
}

.datagenCircle = function(n) {
  x = seq(-1, 1, length = n)
  u = sin(x * pi) + rnorm(n) / 8
  v = cos(x * pi) + rnorm(n) / 8
  return (rbind(u, v))
}

.datagenTwoClassUniv = function(n) {
  y = as.double(runif(n) < 0.5)
  x = y * rnorm(n, mean = -0.2) + (1 - y) * rnorm(n, mean = 0.2)
  return (list(x = x, y = y))
}

.datagenFourClassUniv = function(n) {
  y = as.double(sample(x = 0:3, size = n, replace = T))
  x = (y == 1) * rnorm(n, mean = -0.4) + 
      (y == 2) * rnorm(n, mean = -0.2) +
      (y == 3) * rnorm(n, mean =  0.2) +
      (y == 4) * rnorm(n, mean =  0.4)
  return (list(x = x, y = y))
}

.datagenTwoClassMultiv = function(n) {
  m = 10
  x = matrix(as.double((runif(n * m) < 0.4) + (runif(n * m) < 0.4)), ncol = m)
  y = as.double(xor(rowSums(x[, 1:5] > 0) > 2, rowSums(x[, 6:10] > 0) > 2))
  return (list(x = x, y = y))
}

xdp.gof.competitors = function(x, null.cdf, ...) 
{
  if (!is.numeric(x) && !is.ordered(x)) {
    stop('x is expected to be a numeric or ordered vector')
  }
  
  n = length(x)
  test_type = .UV_GOF_EXISTING
  
  # y is used to store values of the given null CDF at points between every pair of consecutive 
  # observations (and I pad with an extra 0 at the beginning).
  if (!is.numeric(null.cdf)) {
    # If the user only supplies a function as the CDF, then I compute it at the midpoints.
    xs = sort(x)
    y = do.call(null.cdf, list(xs[1:(n - 1)] + diff(xs) / 2, ...))
  }
  y = as.matrix(as.double(c(0, y)), nrow = length(y), ncol = 1)
  
  # Dx and Dy are not used
  Dx = 0
  Dy = 0
  
  w_sum = as.double(0)
  w_max = as.double(2)
  
  # Can make these parameter at some point
  nr.perm = 0
  total.nr.tests = 1
  is.sequential = T
  alpha.hyp = NULL
  alpha0 = NULL
  beta0 = NULL
  eps = NULL
  nr.threads = 1
  tables.wanted = F
  perm.stats.wanted = F  
  
  extra_params = 0
  is_sequential = as.integer(is.sequential)
  
  wald = .configure.wald.sequential(total.nr.tests, is.sequential, alpha.hyp, alpha0, beta0, eps)
  alpha_hyp = as.double(wald$alpha.hyp)
  alpha0 = as.double(wald$alpha0)
  beta0 = as.double(wald$beta0)
  eps = as.double(wald$eps)
  
  nr_perm = as.integer(nr.perm)
  nr_threads = as.integer(nr.threads)
  tables_wanted = as.integer(tables.wanted)
  perm_stats_wanted = as.integer(perm.stats.wanted)
  
  res = .Call('HHG_R_C', test_type, Dx, Dy, y, w_sum, w_max, extra_params, is_sequential, alpha_hyp, alpha0, beta0, eps, nr_perm, nr_threads, tables_wanted, perm_stats_wanted)
  ret = .organize.results(res, n = n, nr.perm, tables.wanted, perm.stats.wanted, grid.len = 0, extra.stats.wanted = F)
  
  names(ret)[names(ret) == 'sum.chisq'] = 'cvm.chisq'
  names(ret)[names(ret) == 'max.chisq'] = 'ks.chisq'
  names(ret)[names(ret) == 'sum.lr'   ] = 'cvm.lr'
  names(ret)[names(ret) == 'max.lr'   ] = 'ks.lr'
  
  return (ret)
}

xdp.ks.competitors = function(x, y, nr.perm = 0, total.nr.tests = 1,
  is.sequential = T, alpha.hyp = NULL, alpha0 = NULL, beta0 = NULL, eps = NULL, 
  nr.threads = 1) 
{
  # Can make these parameters at some point
  w.max = 0
  w.sum = 2
  tables.wanted = F
  perm.stats.wanted = F  
  
  test_type = .UV_KS_CVM_KS
  is_sequential = as.integer(is.sequential)
  
  wald = .configure.wald.sequential(is.sequential, total.nr.tests, alpha.hyp, alpha0, beta0, eps)
  alpha_hyp = as.double(wald$alpha.hyp)
  alpha0 = as.double(wald$alpha0)
  beta0 = as.double(wald$beta0)
  eps = as.double(wald$eps)
  
  nr_perm = as.integer(nr.perm)
  nr_threads = as.integer(nr.threads)
  tables_wanted = as.integer(tables.wanted)
  perm_stats_wanted = as.integer(perm.stats.wanted)
  
  # Dx is used to store ranks of x (a permutation of 1:n)
  Dx = as.matrix(as.double(rank(x, ties.method = 'random')), nrow = length(x), ncol = 1)

  # y is passed as numbers in 0:(K - 1)
  if (is.factor(y)) {
    y = as.numeric(levels(y))[y]
  }
  y = as.matrix(as.double(y), nrow = length(y), ncol = 1)

  # Dy is not used
  Dy = 0
  
  w_sum = as.double(w.sum)
  w_max = as.double(w.max)
  
  extra_params = as.double(0)
  
  res = .Call('HHG_R_C', test_type, Dx, Dy, y, w_sum, w_max, extra_params, is_sequential, alpha_hyp, alpha0, beta0, eps, nr_perm, nr_threads, tables_wanted, perm_stats_wanted)
  ret = .organize.results(res, n = nrow(Dx), nr.perm, tables.wanted, perm.stats.wanted, grid.len = 0, extra.stats.wanted = F)
  
  names(ret)[names(ret) == 'sum.chisq'] = 'cvm.chisq'
  names(ret)[names(ret) == 'max.chisq'] = 'ks.chisq'
  names(ret)[names(ret) == 'sum.lr'   ] = 'cvm.lr'
  names(ret)[names(ret) == 'max.lr'   ] = 'ks.lr'
  
  return (ret)
}

hhg.2.sample.competitors = function(Dx, y, nr.perm = 10000, 
  is.sequential = F, seq.total.nr.tests = 1, seq.alpha.hyp = NULL, seq.alpha0 = NULL, 
  seq.beta0 = NULL, seq.eps = NULL, nr.threads = 0, perm.stats.wanted = F)
{
  # Argument checking is negligent at this point...
  if (!is.vector(y) && !is.factor(y)) {
    stop('y is expected to be a numeric or factor vector with values in {0, 1}')
  }
  if (!is.double(Dx) || !is.matrix(Dx) || nrow(Dx) != ncol(Dx) || nrow(Dx) != length(y)) {
    stop('Dx is expected to be a square matrix of doubles, and must have the same number of rows/cols as the vector y')
  }
  if (nr.perm < 0) {
    stop('nr.perm should not be negative')
  }
  if (is.numeric(y) && !all(y %in% c(0, 1))) {
    stop('y is expected to be a numeric or factor vector with values in {0, 1}')
  }
  if (is.factor(y)) {
    y = as.numeric(levels(y))[y]
    if (!all(y %in% c(0, 1))) {
      stop('y is expected to be a numeric or factor vector with values in {0, 1}')
    }
  }
  
  test_type = .MV_TS_EXISTING
  
  dummy.Dy = 0
  y = as.matrix(as.double(y), nrow = length(y), ncol = 1)
  extra_params = as.double(0)
  is_sequential = as.integer(is.sequential)
  
  wald = .configure.wald.sequential(is.sequential, seq.total.nr.tests, seq.alpha.hyp, seq.alpha0, seq.beta0, seq.eps)
  alpha_hyp = as.double(wald$alpha.hyp)
  alpha0 = as.double(wald$alpha0)
  beta0 = as.double(wald$beta0)
  eps = as.double(wald$eps)
  
  nr_perm = as.integer(nr.perm)
  nr_threads = as.integer(nr.threads)
  tables_wanted = as.integer(F)
  perm_stats_wanted = as.integer(perm.stats.wanted)
  
  res = .Call('HHG_R_C', test_type, Dx, dummy.Dy, y, 0, 2, extra_params, is_sequential, alpha_hyp, alpha0, beta0, eps, nr_perm, nr_threads, tables_wanted, perm_stats_wanted)
  ret = .organize.results(res, n = nrow(Dx), nr.perm, tables_wanted, perm.stats.wanted, grid.len = 0)
  
  ret$edist = ret$sum.chisq
  ret$ht = ret$sum.lr
  ret$perm.pval.edist = ret$perm.pval.hhg.sc
  ret$perm.pval.ht = ret$perm.pval.hhg.sl

  ret$sum.chisq = NULL
  ret$sum.lr = NULL
  ret$max.chisq = NULL
  ret$max.lr = NULL
  ret$perm.pval.hhg.sc = NULL
  ret$perm.pval.hhg.sl = NULL
  ret$perm.pval.hhg.mc = NULL
  ret$perm.pval.hhg.ml = NULL
  
  if (perm.stats.wanted) {
    ret$extras.perm.stats = ret$extras.perm.stats[, 1:2]
    names(ret$extras.perm.stats) = c('edist', 'ht')
  }
  
  return (ret)
}

dynamic.slicing.ks.test = function(x, y, variant = 'ds', lambda = 1) {
  w.max = 0
  w.sum = 2
  tables.wanted = F
  perm.stats.wanted = F  
  nr.threads = 1
  
  if (variant == 'ds') {
    test_type = .UV_KS_DS
  } else if (variant == 'mds') {
    test_type = .UV_KS_MDS
  } else {
    stop('Unexpected variant specified')
  }
  
  nr.perm = 0
  is.sequential = F
  total.nr.tests = 1
  alpha.hyp = NULL
  alpha0 = NULL
  beta0 = NULL
  eps = NULL
    
  wald = .configure.wald.sequential(is.sequential, total.nr.tests, alpha.hyp, alpha0, beta0, eps)
  alpha_hyp = as.double(wald$alpha.hyp)
  alpha0 = as.double(wald$alpha0)
  beta0 = as.double(wald$beta0)
  eps = as.double(wald$eps)

  is_sequential = as.integer(is.sequential)  
  nr_perm = as.integer(nr.perm)
  nr_threads = as.integer(nr.threads)
  tables_wanted = as.integer(tables.wanted)
  perm_stats_wanted = as.integer(perm.stats.wanted)
  
  # y is passed as numbers in 0:(K - 1), sorted according to the order of x
  if (is.factor(y)) {
    y = as.numeric(levels(y))[y]
  }
  y = as.matrix(as.double(y), nrow = length(y), ncol = 1)
  y = as.matrix(y[order(x)])
  
  # Dx and Dy are not used
  Dx = 0
  Dy = 0
  
  w_sum = as.double(w.sum)
  w_max = as.double(w.max)
  
  extra_params = as.double(lambda)

  res = .Call('HHG_R_C', test_type, Dx, Dy, y, w_sum, w_max, extra_params, is_sequential, alpha_hyp, alpha0, beta0, eps, nr_perm, nr_threads, tables_wanted, perm_stats_wanted)
  ret = .organize.results(res, n = nrow(Dx), nr.perm, tables.wanted, perm.stats.wanted, grid.len = 0, extra.stats.wanted = F)
    
  return (ret)
}

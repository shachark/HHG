# Main R file for the Heller-Heller-Gorfine test package

# The extended test of independence
hhg.test.ext = function(Dx, Dy, 
  score = 'ADP', score.control = list(K = 2, w.sum = 0, w.max = 2),
  nr.perm = 10000, is.sequential = F, seq.total.nr.tests = 1, seq.alpha.hyp = NULL, 
  seq.alpha0 = NULL, seq.beta0 = NULL, seq.eps = NULL,
  nr.threads = 0, tables.wanted = F, perm.stats.wanted = F)
{
  # Argument checking is negligent at this point...
  if (!is.double(Dx) || !is.double(Dy) || !is.matrix(Dx) || !is.matrix(Dy) || 
        nrow(Dx) != ncol(Dx) || nrow(Dx) != nrow(Dy) || nrow(Dy) != ncol(Dy)) {
    stop('Dx and Dy are expected to be square matrices of doubles, and must have the same number of rows/cols')
  }
  if (nr.perm < 0) {
    stop('nr.perm should not be negative')
  }
  
  test_type = .MV_IND_HHG_EXTENDED

  # Default params
  w.sum = 0
  w.max = 2
  correct.mi.bias = F
  
  # y is not used
  dummy.y = matrix(0, nrow(Dy), 1)
  
  if (score == 'ADP') {
    w.sum = as.double(score.control$w.sum)
    w.max = as.double(score.control$w.max)
    if (score.control$K == 2) {
      extra_params = as.double(c(.UV_IND_ADP2, 2, correct.mi.bias))
    } else {
      extra_params = as.double(c(.UV_IND_ADP, score.control$K, correct.mi.bias))
    }
  } else if (score == 'DDP') {
    w.sum = as.double(score.control$w.sum)
    w.max = as.double(score.control$w.max)
    if (score.control$K == 2) {
      extra_params = as.double(c(.UV_IND_DDP2, 2, correct.mi.bias))
    } else if (score.control$K == 3) {
      extra_params = as.double(c(.UV_IND_DDP3, 3, correct.mi.bias))
    } else {
      extra_params = as.double(c(.UV_IND_DDP, score.control$K, correct.mi.bias))
    }
  } else {
    stop('Unexepeced score encountered')
  }
  
  is_sequential = as.integer(is.sequential)
  wald = .configure.wald.sequential(is.sequential, seq.total.nr.tests, seq.alpha.hyp, seq.alpha0, seq.beta0, seq.eps)
  alpha_hyp = as.double(wald$alpha.hyp)
  alpha0 = as.double(wald$alpha0)
  beta0 = as.double(wald$beta0)
  eps = as.double(wald$eps)
  
  nr_perm = as.integer(nr.perm)
  nr_threads = as.integer(nr.threads)
  tables_wanted = as.integer(tables.wanted)
  perm_stats_wanted = as.integer(perm.stats.wanted)
  
  res = .Call('HHG_R_C', test_type, Dx, Dy, dummy.y, w.sum, w.max, extra_params, is_sequential, alpha_hyp, alpha0, beta0, eps, nr_perm, nr_threads, tables_wanted, perm_stats_wanted)
  ret = .organize.results(res, n = nrow(Dx), nr.perm, tables.wanted, perm.stats.wanted, grid.len = 0, is.extended.test = T)
  return (ret)
}

# The general test of independence (with or without handling of ties)
hhg.test = function(Dx, Dy, ties = T, w.sum = 0, w.max = 2, nr.perm = 10000, 
  is.sequential = F, seq.total.nr.tests = 1, seq.alpha.hyp = NULL, seq.alpha0 = NULL, 
  seq.beta0 = NULL, seq.eps = NULL, nr.threads = 0, tables.wanted = F, perm.stats.wanted = F)
{
  # Argument checking is negligent at this point...
  if (!is.double(Dx) || !is.double(Dy) || !is.matrix(Dx) || !is.matrix(Dy) || 
        nrow(Dx) != ncol(Dx) || nrow(Dx) != nrow(Dy) || nrow(Dy) != ncol(Dy)) {
    stop('Dx and Dy are expected to be square matrices of doubles, and must have the same number of rows/cols')
  }
  if (w.sum < 0 || w.max < 0) {
    stop('w.sum and w.max should be greater or equal to zero')
  }
  if (nr.perm < 0) {
    stop('nr.perm should not be negative')
  }
  
  if (ties) {
    test_type = .MV_IND_HHG
  } else {
    test_type = .MV_IND_HHG_NO_TIES
  }

  dummy.y = matrix(0, nrow(Dy), 1)
  extra_params = as.double(0)
  is_sequential = as.integer(is.sequential)
  
  wald = .configure.wald.sequential(is.sequential, seq.total.nr.tests, seq.alpha.hyp, seq.alpha0, seq.beta0, seq.eps)
  alpha_hyp = as.double(wald$alpha.hyp)
  alpha0 = as.double(wald$alpha0)
  beta0 = as.double(wald$beta0)
  eps = as.double(wald$eps)
  
  nr_perm = as.integer(nr.perm)
  nr_threads = as.integer(nr.threads)
  tables_wanted = as.integer(tables.wanted)
  perm_stats_wanted = as.integer(perm.stats.wanted)
  
  res = .Call('HHG_R_C', test_type, Dx, Dy, dummy.y, w.sum, w.max, extra_params, is_sequential, alpha_hyp, alpha0, beta0, eps, nr_perm, nr_threads, tables_wanted, perm_stats_wanted)
  ret = .organize.results(res, n = nrow(Dx), nr.perm, tables.wanted, perm.stats.wanted, grid.len = 0)
  return (ret)
}

# The extended K-sample test (with y_i in 0:(K-1))
hhg.test.k.sample.ext = function(Dx, y, 
  score = 'XDP', score.control = list(K = 2, w.sum = 0, w.max = 2),                                 
  nr.perm = 10000, is.sequential = F, seq.total.nr.tests = 1, seq.alpha.hyp = NULL, 
  seq.alpha0 = NULL, seq.beta0 = NULL, seq.eps = NULL, 
  nr.threads = 0, tables.wanted = F, perm.stats.wanted = F)
{
  # Argument checking is negligent at this point...
  if (!is.vector(y)) {
    stop('y is expected to be a vector')
  }
  if (!is.double(Dx) || !is.matrix(Dx) || 
        nrow(Dx) != ncol(Dx) || nrow(Dx) != length(y)) {
    stop('Dx is expected to be a square matrix of doubles, and must have the same number of rows/cols as the vector y')
  }
  if (nr.perm < 0) {
    stop('nr.perm should not be negative')
  }
  if (!is.numeric(y) && !is.factor(y) && !is.logical(y)) {
    stop('y is expected to be a numeric, logical, or factor vector with values in {0, 1, ... K-1}')
  }
  if (is.factor(y)) {
    y = as.numeric(levels(y))[y]
  } else if (is.logical(y)) {
    y = as.numeric(y)
  }
  if (any(y != round(y))) {
    stop('y is expected to be a numeric, logical, or factor vector with values in {0, 1, ... K-1}')
  }
  if (max(y) > 10) {
    warning('the K-sample test is appropriate for small values of K, consider using the general test, implemented by hhg.test')
  }
  
  test_type = .MV_KS_HHG_EXTENDED

  # Default params
  w.sum = 0
  w.max = 2
  correct.mi.bias = F
  
  if (score == 'KW') {
    extra_params = as.double(.UV_KS_KW)
  } else if (score == 'AD') {
    extra_params = as.double(.UV_KS_AD)
  } else if (score == 'CvM/KS') {
    extra_params = as.double(.UV_KS_CVM_KS)
  } else if (score == 'dCov') {
    extra_params = as.double(.UV_KS_DCOV)
  } else if (score == 'DS') {
    extra_params = as.double(c(.UV_KS_DS, score.control$lambda))
  } else if (score == 'MDS') {
    extra_params = as.double(.UV_KS_MDS)
  } else if (score == 'XDP') {
    w.sum = as.double(score.control$w.sum)
    w.max = as.double(score.control$w.max)
    if (score.control$K == 2) {
      extra_params = as.double(c(.UV_KS_XDP2, 2, correct.mi.bias))
    } else if (score.control$K == 3) {
      extra_params = as.double(c(.UV_KS_XDP3, 3, correct.mi.bias))
    } else {
      extra_params = as.double(c(.UV_KS_XDP, score.control$K, correct.mi.bias))
    }
  } else {
    stop('Unexepeced score encountered')
  }

  # y is passed as numbers in 0:(K - 1)
  if (is.factor(y)) {
    y = as.numeric(levels(y))[y]
  }
  y = as.matrix(as.double(y), nrow = length(y), ncol = 1)
  
  # Dy is not used
  Dy = 0
  
  wald = .configure.wald.sequential(is.sequential, seq.total.nr.tests, seq.alpha.hyp, seq.alpha0, seq.beta0, seq.eps)

  nr_perm = as.integer(nr.perm)
  is_sequential = as.integer(is.sequential)
  alpha_hyp = as.double(wald$alpha.hyp)
  alpha0 = as.double(wald$alpha0)
  beta0 = as.double(wald$beta0)
  eps = as.double(wald$eps)
  
  nr_threads = as.integer(nr.threads)
  tables_wanted = as.integer(tables.wanted)
  perm_stats_wanted = as.integer(perm.stats.wanted)
  
  res = .Call('HHG_R_C', test_type, Dx, Dy, y, w.sum, w.max, extra_params, is_sequential, alpha_hyp, alpha0, beta0, eps, nr_perm, nr_threads, tables_wanted, perm_stats_wanted)
  ret = .organize.results(res, n = nrow(Dx), nr.perm, tables.wanted, perm.stats.wanted, grid.len = 0, is.extended.test = T)
  return (ret)
}

# The K-sample test (with y_i in 0:(K-1))
hhg.test.k.sample = function(Dx, y, w.sum = 0, w.max = 2, nr.perm = 10000, 
  is.sequential = F, seq.total.nr.tests = 1, seq.alpha.hyp = NULL, seq.alpha0 = NULL, 
  seq.beta0 = NULL, seq.eps = NULL, nr.threads = 0, tables.wanted = F, perm.stats.wanted = F)
{
  # Argument checking is negligent at this point...
  if (!is.vector(y)) {
    stop('y is expected to be a vector')
  }
  if (!is.double(Dx) || !is.matrix(Dx) || 
        nrow(Dx) != ncol(Dx) || nrow(Dx) != length(y)) {
    stop('Dx is expected to be a square matrix of doubles, and must have the same number of rows/cols as the vector y')
  }
  if (w.sum < 0 || w.max < 0) {
    stop('w.sum and w.max should be greater or equal to zero')
  }
  if (nr.perm < 0) {
    stop('nr.perm should not be negative')
  }
  if (!is.numeric(y) && !is.factor(y) && !is.logical(y)) {
    stop('y is expected to be a numeric, logical, or factor vector with values in {0, 1, ... K-1}')
  }
  if (is.factor(y)) {
    y = as.numeric(levels(y))[y]
  } else if (is.logical(y)) {
    y = as.numeric(y)
  }
  if (any(y != round(y))) {
    stop('y is expected to be a numeric, logical, or factor vector with values in {0, 1, ... K-1}')
  }
  if (max(y) > 10) {
    warning('the K-sample test is appropriate for small values of K, consider using the general test, implemented by hhg.test')
  }
  
  test_type = .MV_KS_HHG
  
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
  tables_wanted = as.integer(tables.wanted)
  perm_stats_wanted = as.integer(perm.stats.wanted)
    
  res = .Call('HHG_R_C', test_type, Dx, dummy.Dy, y, w.sum, w.max, extra_params, is_sequential, alpha_hyp, alpha0, beta0, eps, nr_perm, nr_threads, tables_wanted, perm_stats_wanted)
  ret = .organize.results(res, n = nrow(Dx), nr.perm, tables.wanted, perm.stats.wanted, grid.len = 0)
  return (ret)
}

# The 2-sample test (with y_i in {0, 1})
hhg.test.2.sample = function(Dx, y, w.sum = 0, w.max = 2, nr.perm = 10000, 
  is.sequential = F, seq.total.nr.tests = 1, seq.alpha.hyp = NULL, seq.alpha0 = NULL, 
  seq.beta0 = NULL, seq.eps = NULL, nr.threads = 0, tables.wanted = F, perm.stats.wanted = F)
{
  # Argument checking is negligent at this point...
  if (!is.numeric(y) && !is.factor(y) && !is.logical(y)) {
    stop('y is expected to be a numeric, logical, or factor vector with values in {0, 1, ... K-1}')
  }
  if (!is.double(Dx) || !is.matrix(Dx) || nrow(Dx) != ncol(Dx) || nrow(Dx) != length(y)) {
    stop('Dx is expected to be a square matrix of doubles, and must have the same number of rows/cols as the vector y')
  }
  if (w.sum < 0 || w.max < 0) {
    stop('w.sum and w.max should be greater or equal to zero')
  }
  if (nr.perm < 0) {
    stop('nr.perm should not be negative')
  }
  if (is.factor(y)) {
    y = as.numeric(levels(y))[y]
  } else if (is.logical(y)) {
    y = as.numeric(y)
  }
  if (!all(y %in% c(0, 1))) {
    stop('y is expected to be a numeric, logical, or factor vector with values in {0, 1}')
  }
  
  test_type = .MV_TS_HHG
  
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
  tables_wanted = as.integer(tables.wanted)
  perm_stats_wanted = as.integer(perm.stats.wanted)
    
  res = .Call('HHG_R_C', test_type, Dx, dummy.Dy, y, w.sum, w.max, extra_params, is_sequential, alpha_hyp, alpha0, beta0, eps, nr_perm, nr_threads, tables_wanted, perm_stats_wanted)
  ret = .organize.results(res, n = nrow(Dx), nr.perm, tables.wanted, perm.stats.wanted, grid.len = 0)
  return (ret)
}

# The univariate distribution-free test
xdp.test = function(x, y, variant = 'DDP', K = 3, correct.mi.bias = F, w.sum = 0, w.max = 2) {
  if (variant == 'DDP') {
    if (K != as.integer(K) || K < 2) {
      stop('K must be an integer greater than 1')
    } else if (K == 2) {
      x.variant = 'spr.obs'
    } else if (K == 3) {
      x.variant = 'ppr.33.obs'
    } else if (K == 4) {
      x.variant = 'tpr.obs'
    } else {
      x.variant = 'ddp.obs'
    }
  } else if (variant == 'ADP') {
    if (K != as.integer(K) || K < 2) {
      stop('K must be an integer greater than 1')
    } else if (K == 2) {
      x.variant = 'spr.all'
    } else if (K == 3) {
      # One could use 'ppr.33.all' that has the same complexity, but in practice it is much slower
      x.variant = 'ddp.all'
    } else if (K == 4) {
      # tpr.all is too time consuming
      x.variant = 'ddp.all'
    } else {
      x.variant = 'ddp.all'
    }
  }

  ret = .hhg.test.udfree(x = x, y = y, variant = x.variant, K = K, correct.mi.bias = correct.mi.bias)
  
  if (((variant == 'DDP') && (K > 4)) || ((variant == 'ADP') && (K > 2))) {
    ret$max.chisq = NA
    ret$max.lr    = NA

    if (!is.null(ret$perm.pval.hhg.mc)) {
      ret$perm.pval.hhg.mc = NA
      ret$perm.pval.hhg.ml = NA
    }
  }
    
  return (ret)
}

.hhg.test.udfree = function(x, y, variant = 'ppr.33.obs', w.sum = 0, w.max = 2,
  nr.perm = 0, K = 3, correct.mi.bias = F, total.nr.tests = 1, 
  is.sequential = F, alpha.hyp = NULL, alpha0 = NULL, beta0 = NULL, eps = NULL, 
  nr.threads = 1, tables.wanted = F, perm.stats.wanted = F)
{
  if (!is.vector(y)) {
    stop('y is expected to be a vector')
  }
  if (!is.vector(y) || length(x) != length(y)) {
    stop('x is expected to be a vector, and must have the same length as the vector y')
  }
  if (nr.perm < 0) {
    stop('nr.perm should not be negative')
  }
  if (!is.numeric(y) && !is.ordered(y)) {
    stop('y is expected to be a numeric or ordered vector')
  }
  if (!is.numeric(x) && !is.ordered(x)) {
    stop('x is expected to be a numeric or ordered vector')
  }
  if (K < 2 || K > length(x)) {
    stop('K should be strictly between 2 and length(x)')
  }

  if (variant == 'spr.obs') {
    test_type = .UV_IND_DDP2
  } else if (variant == 'spr.all') {
    test_type = .UV_IND_ADP2
  } else if (variant == 'ppr.22.obs') {
    test_type = .UV_IND_DDP3_C
  } else if (variant == 'ppr.22.all') {
    test_type = .UV_IND_ADP3_C
  } else if (variant == 'ppr.33.obs') {
    test_type = .UV_IND_DDP3
  } else if (variant == 'ppr.33.all') {
    test_type = .UV_IND_ADP3
  } else if (variant == 'tpr.obs') {
    test_type = .UV_IND_DDP4
  } else if (variant == 'tpr.all') {
    test_type = .UV_IND_ADP4
  } else if (variant == 'ddp.obs') {
    test_type = .UV_IND_DDP
  } else if (variant == 'ddp.all') {
    test_type = .UV_IND_ADP
  } else {
  	stop('Unexpected variant specified.')
  }
  
  # Dx is used to store x
  Dx = as.matrix(as.double(rank(x, ties.method = 'random')), nrow = length(x), ncol = 1)
  y  = as.matrix(as.double(rank(y, ties.method = 'random')), nrow = length(y), ncol = 1)

  # Dy is not used
  Dy = 0
  
  # For historical reasons, the high-k DDP and ADP variants work on 1-based ranks, while the small-k
  # variants work on 0-based ranks.
  if (!(variant == 'ddp.obs') && !(variant == 'ddp.all')) {
    Dx = Dx - 1
    y = y - 1
  }
  
  w_sum = as.double(w.sum)
  w_max = as.double(w.max)

  extra_params = as.double(c(K, correct.mi.bias))
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
  
  res = .Call('HHG_R_C', test_type, Dx, Dy, y, w_sum, w_max, extra_params, is_sequential, alpha_hyp, alpha0, beta0, eps, nr_perm, nr_threads, tables_wanted, perm_stats_wanted)
  ret = .organize.results(res, n = nrow(Dx), nr.perm, tables.wanted, perm.stats.wanted, grid.len = 0)
  return (ret)
}

# The k-sample version of the XDP test
xdp.test.k.sample = function(x, y, ddp.K = 3, w.sum = 0, w.max = 2) 
{
  # The interface may need more work..
  
  if (!is.numeric(x) && !is.ordered(x)) {
    stop('x is expected to be a numeric or ordered vector')
  }
  if (!is.numeric(y) && !is.factor(y) && !is.logical(y)) {
    stop('y is expected to be a numeric, logical, or factor vector with values in {0, 1, ... K-1}')
  }
  if (is.factor(y)) {
    y = as.numeric(levels(y))[y]
  } else if (is.logical(y)) {
    y = as.numeric(y)
  }
  if (any(y != round(y))) {
    stop('y is expected to be a numeric, logical, or factor vector with values in {0, 1, ... K-1}')
  }
  if (max(y) > 10) {
    warning('the K-sample test is appropriate for small values of K, consider using the general test, implemented by xdp.test')
  }
  if (!is.vector(y) || length(x) != length(y)) {
    stop('x is expected to be a vector, and must have the same length as the vector y')
  }
  if (ddp.K != as.integer(ddp.K) || ddp.K < 2 || ddp.K > length(x)) {
    stop('K must be an integer between 2 and length(x)')
  } 
  
  if (ddp.K == 2) {
    test_type = .UV_KS_XDP2
  } else if (ddp.K == 3) {
    test_type = .UV_KS_XDP3
  } else {
    test_type = .UV_KS_XDP
  }

  # Dx is used to store ranks of x (a permutation of 1:n)
  Dx = as.matrix(as.double(rank(x, ties.method = 'random')), nrow = length(x), ncol = 1)

  # y is passed as numbers in 0:(K - 1)
  y = as.matrix(as.double(y), nrow = length(y), ncol = 1)

  # Dy is not used
  Dy = 0
  
  w_sum = as.double(w.sum)
  w_max = as.double(w.max)

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
  
  extra_params = as.double(c(ddp.K))
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
  
  res = .Call('HHG_R_C', test_type, Dx, Dy, y, w_sum, w_max, extra_params, is_sequential, alpha_hyp, alpha0, beta0, eps, nr_perm, nr_threads, tables_wanted, perm_stats_wanted)
  ret = .organize.results(res, n = nrow(Dx), nr.perm, tables.wanted, perm.stats.wanted, grid.len = 0)
  
  if (ddp.K > 3) {
    ret$max.chisq = NA
    ret$max.lr    = NA
    
    if (!is.null(ret$perm.pval.hhg.mc)) {
      ret$perm.pval.hhg.mc = NA
      ret$perm.pval.hhg.ml = NA
    }
  }
  
  return (ret)
}

# The 1-sample (goodness-of-fit) version of the XDP test
xdp.test.gof = function(x, null.cdf, K = 3, w.sum = 0, w.max = 2, ...) 
{
  if (!is.numeric(x) && !is.ordered(x)) {
    stop('x is expected to be a numeric or ordered vector')
  }
  if (K != as.integer(K) || K < 2 || K > length(x)) {
    stop('K must be an integer between 2 and length(x)')
  } 
  
  if (K == 2) {
    test_type = .UV_GOF_XDP2
  } else if (K == 3) {
    test_type = .UV_GOF_XDP3
  } else {
    test_type = .UV_GOF_XDP
  }
  
  n = length(x)
  
  # y is used to store values of the given null CDF at points between every pair of consecutive 
  # observations (and I pad next with an extra 0 at the beginning).
  if (!is.numeric(null.cdf)) {
    # If the user only supplies a function as the CDF, then I compute it at the midpoints.
    xs = sort(x)
    y = do.call(null.cdf, list(xs[1:(n - 1)] + diff(xs) / 2, ...))
  }
  
  if (any(y > 1 || y < 0)) {
    stop('CDF values at sample points lie outside [0, 1]')
  }

  if (any(y == 1 || y == 0)) {
    warning('The given sample seems to lie outside the support of the (assumed continuous) null distribution. Rejecting.')
    ret = list()
    ret$sum.chisq = Inf
    ret$sum.lr    = Inf
    ret$max.chisq = Inf
    ret$max.lr    = Inf
    return (ret)
  }

  y = as.matrix(as.double(c(0, y)), nrow = length(y), ncol = 1)
  
  # Dx and Dy are not used
  Dx = 0
  Dy = 0
  
  w_sum = as.double(w.sum)
  w_max = as.double(w.max)
  
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
  
  extra_params = as.double(K)
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
  
  res = .Call('HHG_R_C', test_type, Dx, Dy, y, w_sum, w_max, extra_params, is_sequential, alpha_hyp, alpha0, beta0, eps, nr_perm, nr_threads, tables_wanted, perm_stats_wanted)
  ret = .organize.results(res, n = n, nr.perm, tables.wanted, perm.stats.wanted, grid.len = 0)
  
  if (K > 3) {
    ret$max.chisq = NA
    ret$max.lr    = NA
    
    if (!is.null(ret$perm.pval.hhg.mc)) {
      ret$perm.pval.hhg.mc = NA
      ret$perm.pval.hhg.ml = NA
    }
  }
  
  return (ret)
}

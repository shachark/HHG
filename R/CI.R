# HHG conditional independence tests

# NOTE: unlike the unconditional DDP/ADP tests, these are not distribution-free
# because they depend on the conditional marginals (and for multivariate z, also
# on its marginal distribution)
xdp.ci.test = function(x, y, Dz = NULL, z = NULL, variant = 'ddp.all', K = 5, correct.mi.bias = F,
  kern = list(type = 'NN', h = ceiling(length(x) / 10)), 
  lsb.kern = list(type = 'NN', h = ceiling(length(x) / 10)), 
  nr.perm = 10000, total.nr.tests = 1, is.sequential = T, alpha.hyp = NULL, 
  alpha0 = NULL, beta0 = NULL, eps = NULL, nr.threads = 0, tables.wanted = F, 
  perm.stats.wanted = F)
{
  # Argument checking is negligent at this point...
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
  if (is.null(Dz) && is.null(z)) {
    stop('Either Dz or z must be provided')
  }
  if (!is.null(z) && (!is.vector(z) || !is.double(z) || length(z) != length(x) || is.unsorted(z))) {
    stop('z has to be a sorted vector of doubles whose length is length(x)')
  }
  if (!is.null(Dz) && (!is.double(Dz) || !is.matrix(Dz) || any(dim(Dz) != c(length(x), length(x))))) {
    stop('Dz has to be a distance matrix with the same dimensions as x')
  }
  if (K < 2 || K > length(x)) {
    stop('K should be strictly between 2 and length(x)')
  }
  if (!is.list(lsb.kern) || !is.list(kern)) {
    stop('kern and lsb.kern must be lists specifying the kernels to use')
  }
  if (variant != 'ddp.all') {
    stop('For the moment, I have only implemented the variant \'ddp.all\'')
  }
  if (perm.stats.wanted && nr.perm == 0) {
    stop('perm.stats.wanted must be false with nr.perm == 0')
  }
  if (lsb.kern$type != 'NN') {
    stop('Only the nearest neighbor kernel is supported for the locally smoothed bootstrap at this time')
  } else {
    if (lsb.kern$h < 0 || lsb.kern$h > length(x) - 1) {
      stop('lsb.kern$h should be greater than zero and lass than length(x) - 1')
    }
  }
  
  # Maybe add these later to implementation and as parameters here
  w.sum = 0
  w.max = 0
  
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
    
  # Dx and Dy are used to store x and y, respectively
  # FIXME computing the ranks is redundant now
  Dx = as.matrix(as.double(rank(x, ties.method = 'random')), nrow = length(x), ncol = 1)
  Dy = as.matrix(as.double(rank(y, ties.method = 'random')), nrow = length(y), ncol = 1)
  
  if (is.null(Dz)) {
    # Univariate z
    stop('For the moment, I did not implement the univariate z optimized test; pass Dz instead.')
  } else {
    # Multivariate z
    if (kern$type == 'NN') {
      test_type = .CI_UDF_ADP_MVZ_NN
      if (kern$h > nrow(Dx) - 1) {
        stop('Kernel neighborhood size must be < number of samples - 1')
      }
      extra_params = as.double(c(K, correct.mi.bias, kern$h, lsb.kern$h))
      res = .Call('HHG_R_C', test_type, Dx, Dy, Dz, w.sum, w.max, extra_params, is_sequential, alpha_hyp, alpha0, beta0, eps, nr_perm, nr_threads, tables_wanted, perm_stats_wanted)
    } else if (kern$type == 'Gaussian') {
      stop('For the moment, I did not implement the Gaussian kernel; use the NN kernel instead.')
    } else {
      stop('Unexpected kernel type given')
    }
  }

  n = nrow(Dx)
  if (kern$type != 'grid NN') {
    h.vec.len = 0
  }
  ret = .organize.results(res, n, nr.perm, tables.wanted, perm.stats.wanted, h.vec.len)
  return (ret)
}

hhg.ci.test = function(Dx, Dy, Dz = NULL, z = NULL,
  kern = list(type = 'NN', h = ceiling(nrow(Dx) / 10)), 
  lsb.kern = list(type = 'NN', h = ceiling(nrow(Dx) / 10)), 
  ties = T, w.sum = 0, w.max = 2, nr.perm = 10000, total.nr.tests = 1, 
  is.sequential = T, alpha.hyp = NULL, alpha0 = NULL, beta0 = NULL, eps = NULL,
  nr.threads = 0, tables.wanted = F, perm.stats.wanted = F)
{
  # Argument checking is negligent at this point...
  if (!is.double(Dx) || !is.double(Dy) || !is.matrix(Dx) || !is.matrix(Dy) || 
        nrow(Dx) != ncol(Dx) || nrow(Dx) != nrow(Dy) || nrow(Dy) != ncol(Dy)) {
    stop('Dx and Dy are expected to be square matrices of doubles, and must have the same number of rows/cols')
  }
  if (is.null(Dz) && is.null(z)) {
    stop('Either Dz or z must be provided')
  }
  if (!is.null(z) && (!is.vector(z) || !is.double(z) || length(z) != nrow(Dx) || is.unsorted(z))) {
    stop('z has to be a sorted vector of doubles whose length is nrow(Dx)')
  }
  if (!is.null(Dz) && (!is.double(Dz) || !is.matrix(Dz) || any(dim(Dz) != dim(Dx)))) {
    stop('Dz has to be a distance matrix with the same dimensions as Dx')
  }
  if (w.sum < 0 || w.max < 0) {
    stop('w.sum and w.max should be greater or equal to zero')
  }
  if (nr.perm < 0) {
    stop('nr.perm should not be negative')
  }
  if (!is.list(lsb.kern) || !is.list(kern)) {
    stop('kern and lsb.kern must be lists specifying the kernels to use')
  }
  if (perm.stats.wanted && nr.perm == 0) {
    stop('perm.stats.wanted must be false with nr.perm == 0')
  }
  if (lsb.kern$type != 'NN') {
    stop('Only the nearest neighbor kernel is supported for the locally smoothed bootstrap at this time')
  }
  
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
  
  if (is.null(Dz)) {
    # Univariate z
    if (kern$type == 'NN') {
      test_type = .CI_UVZ_NN
      if (kern$h > nrow(Dx) - 1) {
        stop('Kernel neighborhood size must be < number of samples - 1')
      }
      extra_params = as.double(c(kern$h, lsb.kern$h))
      res = .Call('HHG_R_C', test_type, Dx, Dy, as.matrix(z), w.sum, w.max, extra_params, is_sequential, alpha_hyp, alpha0, beta0, eps, nr_perm, nr_threads, tables_wanted, perm_stats_wanted)
    } else if (kern$type == 'Gaussian') {
      test_type = .CI_UVZ_GAUSSIAN
      extra_params = as.double(c(kern$sig, lsb.kern$h))
      res = .Call('HHG_R_C', test_type, Dx, Dy, as.matrix(z), w.sum, w.max, extra_params, is_sequential, alpha_hyp, alpha0, beta0, eps, nr_perm, nr_threads, tables_wanted, perm_stats_wanted)
    } else {
      stop('Unexpected kernel type given')
    }
  } else {
    # Multivariate z
    if (kern$type == 'NN') {
      test_type = .CI_MVZ_NN
      if (kern$h > nrow(Dx) - 1) {
        stop('Kernel neighborhood size must be < number of samples - 1')
      }
      extra_params = as.double(c(kern$h, lsb.kern$h))
      res = .Call('HHG_R_C', test_type, Dx, Dy, Dz, w.sum, w.max, extra_params, is_sequential, alpha_hyp, alpha0, beta0, eps, nr_perm, nr_threads, tables_wanted, perm_stats_wanted)
    } else if (kern$type == 'Gaussian') {
      test_type = .CI_MVZ_GAUSSIAN
      extra_params = as.double(c(kern$sig, lsb.kern$h))
      res = .Call('HHG_R_C', test_type, Dx, Dy, Dz, w.sum, w.max, extra_params, is_sequential, alpha_hyp, alpha0, beta0, eps, nr_perm, nr_threads, tables_wanted, perm_stats_wanted)
    } else if (kern$type == 'grid NN') {
      test_type = .CI_MVZ_NN_GRID_BW
      h.vec.len = length(kern$h.vec)
      extra_params = as.double(c(h.vec.len, lsb.kern$h, kern$h.vec))
      res = .Call('HHG_R_C', test_type, Dx, Dy, Dz, w.sum, w.max, extra_params, is_sequential, alpha_hyp, alpha0, beta0, eps, nr_perm, nr_threads, tables_wanted, perm_stats_wanted)
    } else {
      stop('Unexpected kernel type given')
    }
  }

  n = nrow(Dx)
  if (kern$type != 'grid NN') {
    h.vec.len = 0
  }
  ret = .organize.results(res, n, nr.perm, tables.wanted, perm.stats.wanted, h.vec.len)
  return (ret)
}

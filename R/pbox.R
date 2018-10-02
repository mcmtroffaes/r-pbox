#' An S4 class to represent a discrete cumulative distribution function.
#'
#' This class represents a probability mass function over n points on
#' the real line, where the probability masses are uniformly
#' distributed.
#'
#' @slot points A sorted numeric vector of length one or more
#' @export
setClass(
  "DiscreteCdf",
  representation(points = "numeric"),
  prototype(points = 0),
  validity = function(object) {
    if (length(object@points) < 1L) {
      "must contain at least one point"
    } else if (is.unsorted(object@points)) {
      "points must be sorted"
    } else {
      TRUE
    }
  }
)

#' Create a new DiscreteCdf from a given set of sorted points.
#'
#' @param points A sorted numeric vector of length one or more
#' @export
discrete_cdf <- function(points) new("DiscreteCdf", points = points)

#' Number of points.
#' @export
setMethod("length", "DiscreteCdf", function(x) length(x@points))

#' Distribution mean.
#' @export
setMethod("mean", "DiscreteCdf", function(x) mean(x@points))

#' Distribution median.
#' @export
setMethod("median", "DiscreteCdf", function(x) median(x@points))

#' Create cumulative distribution function.
eval_discrete_cdf <- function(cdf) {
  function(x) {
    # (note: implementation is simple but slow)
    n <- length(cdf)
    for (i in 1:n) {
      if (x < cdf@points[i]) return((i - 1) / n)
    }
    return(1)
  }
}

# a p-box is represented by two cdf sequences
# u@points[1] ... u@points[n]
# l@points[1] ... l@points[n]

#' An S4 class to represent a discrete p-box.
#'
#' This class represents the set of all cumulative distribution
#' functions contained between two discrete cumulative distribution
#' functions.
#'
#' @slot u Upper discrete cumulative distribution function
#' @slot l Lower discrete cumulative distribution function
#'
#' @export
setClass(
  "DiscretePBox",
  representation(u = "DiscreteCdf", l = "DiscreteCdf"),
  validity = function(object) {
    if (length(object@u) != length(object@l)) {
      "ucdf and lcdf must have same length"
    } else if (any(object@u@points > object@l@points)) {
      "ucdf must be left of lcdf"
    } else {
      TRUE
    }
  }
)

#' Create a new DiscretePBox from two sets of sorted points.
#'
#' @param u A sorted numeric vector of length one or more,
#'     representing the upper cdf
#' @param l A sorted numeric vector of length one or more,
#'     representing the lower cdf
#'
#' @export
discrete_pbox <- function(upoints, lpoints)
  new("DiscretePBox", u = discrete_cdf(upoints), l = discrete_cdf(lpoints))

#' Number of points.
#'
#' @export
setMethod("length", "DiscretePBox", function(x) length(x@u@points))

#' Create a discrete pbox bounding a given normal distribution.
#'
#' @param mean Mean
#' @param sd Standard deviation
#' @param bot Quantile for left cut-off point (should be less than 1/steps)
#' @param top Quantile for right cut-off point (should be more than 1-1/steps)
#' @param steps Number of discretization steps
#'
#' @export
discrete_pbox_norm <- function(mean = 0, sd = 1, bot = 0.001, top = 0.999, steps = 200) {
  stopifnot(bot >= 0, top <= 1)
  ps <- c(
    min(bot, 1 / steps),
    seq(0, 1, length.out = steps + 1)[2:steps],
    max(top, 1 - 1 / steps)
  )
  points <- sapply(ps, function(x) qnorm(x, mean = mean, sd = sd))
  discrete_pbox(upoints = points[-(steps + 1)], lpoints = points[-1])
}

#' Create a discrete pbox bounding a given uniform distribution.
#'
#' @param min Left side of the range
#' @param max Right side of the range
#' @param steps Number of discretization steps
#'
#' @export
discrete_pbox_unif <- function(min = 0, max = 1, steps = 200) {
  ps <- seq(0, 1, length.out = steps + 1)
  points <- sapply(ps, function(x) qunif(x, min = min, max = max))
  discrete_pbox(upoints = points[-(steps + 1)], lpoints = points[-1])
}

eval_discrete_pbox_u <- function(pbox) eval_discrete_cdf(pbox@u)

eval_discrete_pbox_l <- function(pbox) eval_discrete_cdf(pbox@l)

setMethod("mean", "DiscretePBox", function(x) c(mean(x@u), mean(x@l)))

setMethod("median", "DiscretePBox", function(x) c(median(x@u), median(x@l)))

# Williamson & Downs, p. 126-127
discrete_pbox_central_moment <- function(pbox, k) {
  mu <- mean(pbox@u@points)
  ml <- mean(pbox@l@points)
  points1 <- pbox@l@points - mu
  points2 <- pbox@u@points - mu
  points3 <- pbox@l@points - ml
  points4 <- pbox@u@points - ml
  mml <- pmin(points1, points2, points3, points4)**k
  mmu <- pmax(points1, points2, points3, points4)**k
  c(mean(pmin(mml, mmu)), mean(pmax(mml, mmu)))
}

setMethod(
  "sd", "DiscretePBox",
  function(x) sqrt(discrete_pbox_central_moment(x, 2))
)

# unary negation
setMethod(
  "-", signature(e1 = "DiscretePBox", e2 = "missing"),
  function(e1) discrete_pbox(-rev(e1@l@points), -rev(e1@u@points))
)

setMethod(
  "+", signature(e1 = "DiscretePBox", e2 = "numeric"),
  function(e1, e2) discrete_pbox(e1@u@points + e2, e1@l@points + e2)
)

setMethod(
  "+", signature(e1 = "numeric", e2 = "DiscretePBox"),
  function(e1, e2) e2 + e1
)

setMethod(
  "-", signature(e1 = "DiscretePBox", e2 = "numeric"),
  function(e1, e2) e1 + (-e2)
)

setMethod(
  "-", signature(e1 = "numeric", e2 = "DiscretePBox"),
  function(e1, e2) e1 + (-e2)
)

setMethod(
  ">=", signature(e1 = "DiscretePBox", e2 = "numeric"),
  function(e1, e2) e1@u@points[1] >= e2
)

setMethod(
  ">", signature(e1 = "DiscretePBox", e2 = "numeric"),
  function(e1, e2) e1@u@points[1] > e2
)

setMethod(
  "<=", signature(e1 = "DiscretePBox", e2 = "numeric"),
  function(e1, e2) e1@l@points[length(e1)] <= e2
)

setMethod(
  "<", signature(e1 = "DiscretePBox", e2 = "numeric"),
  function(e1, e2) e1@l@points[length(e1)] < e2
)

setMethod(
  "/", signature(e1 = "numeric", e2 = "DiscretePBox"),
  function(e1, e2) {
    stopifnot(e1 >= 0, e2 > 0)
    discrete_pbox(rev(e1 / e2@l@points), rev(e1 / e2@u@points))
  }
)

setMethod(
  "/", signature(e1 = "DiscretePBox", e2 = "numeric"),
  function(e1, e2) {
    stopifnot(e2 > 0)
    discrete_pbox(e1@u@points / e2, e1@l@points / e2)
  }
)

# apply function on every combination of elements of xs1 and xs2, and sort
sortedfunc <- function(func, xs1, xs2) {
  ys1 <- rep(xs1, length(xs2))
  ys2 <- rep(xs2, rep(length(xs1), length(xs2)))
  sort(func(ys1, ys2))
}

# implementation of Williamson & Downs, Figure 14, page 127
discrete_pbox_convolution <- function(func) function(pbox1, pbox2) {
    n <- length(pbox1)
    stopifnot(n == length(pbox2))
    ixs <- (0:(n - 1)) * n
    discrete_pbox(
      upoints = sortedfunc(func, pbox1@u@points, pbox2@u@points)[ixs + 1],
      lpoints = sortedfunc(func, pbox1@l@points, pbox2@l@points)[ixs + n]
    )
  }

setGeneric("%iadd%", function(e1, e2) standardGeneric("%iadd%"))
setGeneric("%imul%", function(e1, e2) standardGeneric("%imul%"))
setGeneric("%isub%", function(e1, e2) standardGeneric("%isub%"))
setGeneric("%idiv%", function(e1, e2) standardGeneric("%idiv%"))

setMethod("%iadd%", "DiscretePBox",
          function(e1, e2) discrete_pbox_convolution(`+`)(e1, e2))

setMethod("%imul%", "DiscretePBox", function(e1, e2) {
  stopifnot(e1 >= 0, e2 >= 0)
  discrete_pbox_convolution(`*`)(e1, e2)
})

setMethod("%isub%", "DiscretePBox", function(e1, e2) e1 + (-e2))

setMethod("%idiv%", "DiscretePBox", function(e1, e2) e1 * (1 / e2))

# implementation of Williamson & Downs, page 120-121 (cases i and ii)
discrete_pbox_frechet <- function(func) function(pbox1, pbox2) {
    n <- length(pbox1)
    stopifnot(n == length(pbox2))
    discrete_pbox(
      upoints = sapply(1:n, function(i)
        max(func(pbox1@u@points[1:i], pbox2@u@points[i:1]))),
      lpoints = sapply(1:n, function(i)
        min(func(pbox1@l@points[i:n], pbox2@l@points[n:i])))
    )
  }

setGeneric("%fadd%", function(e1, e2) standardGeneric("%fadd%"))
setGeneric("%fmul%", function(e1, e2) standardGeneric("%fmul%"))
setGeneric("%fsub%", function(e1, e2) standardGeneric("%fsub%"))
setGeneric("%fdiv%", function(e1, e2) standardGeneric("%fdiv%"))

setMethod("%fadd%", "DiscretePBox", function(e1, e2) discrete_pbox_frechet(`+`)(e1, e2))

setMethod("%fmul%", "DiscretePBox", function(e1, e2) {
  stopifnot(e1 >= 0, e2 >= 0)
  discrete_pbox_frechet(`*`)(e1, e2)
})

setMethod("%fsub%", "DiscretePBox", function(e1, e2) e1 %fadd% (-e2))

setMethod("%fdiv%", "DiscretePBox", function(e1, e2) e1 %fmul% (1 / e2))

setMethod("plot", "DiscretePBox", function(x, xlab = "x", ylab = "cumulative probability", ...) {
  n <- length(x)
  cs <- rep(2, n)
  uxs <- c(rep(x@u@points, cs), x@l@points[length(x)])
  uys <- rep((0:n) / n, c(1, cs))
  lxs <- c(x@u@points[1], rep(x@l@points, cs))
  lys <- rep((0:n) / n, c(cs, 1))
  plot(uxs, uys, col = 2, type = "l", xlab = xlab, ylab = ylab, ...)
  lines(lxs, lys, col = 1, ...)
  legend("topleft",
    legend = c("upper cdf", "lower cdf"),
    col = c(2, 1),
    lty = 1
  )
})

test1 <- function() {
  pb <- discrete_pbox_norm(steps = 10)
  plot(pb)
  xs <- seq(-5, 5, 0.01)
  lines(xs, sapply(xs, pnorm), type = "l", col = 3)
}

test2 <- function() {
  pb <- discrete_pbox_unif(min = -2, max = 3, steps = 10)
  plot(pb)
  xs <- seq(-4, 4, 0.01)
  lines(xs, sapply(xs, function(x) punif(x, min = -2, max = 3)), type = "l", col = 3)
}

# W&D figure 19
test3 <- function() {
  pb1 <- discrete_pbox_unif(min = 1, max = 2, steps = 40)
  pb2 <- discrete_pbox_unif(min = 1, max = 2, steps = 40)
  pb3 <- pb1 %iadd% pb2
  plot(pb3)
}

# W&D figure 21
test4 <- function() {
  pb1 <- discrete_pbox_unif(min = 0, max = 1, steps = 40)
  pb2 <- discrete_pbox_unif(min = 0, max = 1, steps = 40)
  pb3 <- pb1 %fadd% pb2
  plot(pb3)
}

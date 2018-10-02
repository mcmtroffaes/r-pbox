# a cdf sequence is a sequence of numbers x1 ... xn so that
# cdf(x)=0   for x<x1
# cdf(x)=1/n for x1<=x<x2
# cdf(x)=2/n for x2<=x<x3
# ...
# cdf(x)=1   for xn<=x

check_discrete_cdf = function(object) {
    if (length(object@points) < 1L) "must contain at least one point"
    else if (is.unsorted(object@points)) "points must be sorted"
    else TRUE
}

setClass("DiscreteCdf", representation(points="numeric"),
         list(points=0),
         validity=check_discrete_cdf)

setMethod("show", "DiscreteCdf", function(object) cat("points for cdf =", object@points, "\n"))

discrete_cdf = function(points) new("DiscreteCdf", points=points)

setMethod("length", "DiscreteCdf", function(x) length(x@points))

setMethod("mean", "DiscreteCdf", function(x) mean(x@points))

setMethod("median", "DiscreteCdf", function(x) median(x@points))

# calculate the cdf as a function of x according the above
# discrete representation (note: implementation is simple but slow)
eval_discrete_cdf = function(cdf) function(x) {
    n = length(cdf)
    for (i in 1:n) { if (x < cdf@points[i]) return ((i - 1) / n) }
    return(1)
}

# a p-box is represented by two cdf sequences
# u@points[1] ... u@points[n]
# l@points[1] ... l@points[n]

check_discrete_pbox = function(object) {
    if (length(object@u) != length(object@l)) "ucdf and lcdf must have same length"
    else if (any(object@u@points > object@l@points)) "ucdf must be left of lcdf"
    else TRUE
}

setClass("DiscretePBox", representation(u="DiscreteCdf", l="DiscreteCdf"),
         validity=check_discrete_pbox)

discrete_pbox = function(upoints, lpoints)
    new("DiscretePBox", u=discrete_cdf(upoints), l=discrete_cdf(lpoints))

setMethod("length", "DiscretePBox", function(x) length(x@u@points))

setMethod("show", "DiscretePBox", function(object) {
    cat("points for upper cdf =", object@u@points, "\n");
    cat("points for lower cdf =", object@l@points, "\n")
})

# bot = bottom quantile (cut off at 1/n)
# top = top quantile (cut off at 1-1/n)
discrete_pbox_norm = function(mean=0, sd=1, bot=0.001, top=0.999, n=200) {
    stopifnot(bot > 0, top < 1)
    ps = c(min(bot, 1/n), seq(0, 1, length.out=n+1)[2:n], max(top, 1-1/n))
    points = sapply(ps, function(x) qnorm(x, mean=mean, sd=sd))
    discrete_pbox(upoints=points[-(n+1)], lpoints=points[-1])
}

discrete_pbox_unif = function(min=0, max=1, n=200) {
    ps = seq(0, 1, length.out=n+1)
    points = sapply(ps, function(x) qunif(x, min=min, max=max))
    discrete_pbox(upoints=points[-(n+1)], lpoints=points[-1])
}

eval_discrete_pbox_u = function(pbox) eval_discrete_cdf(pbox@u)

eval_discrete_pbox_l = function(pbox) eval_discrete_cdf(pbox@l)

setMethod("mean", "DiscretePBox", function(x) c(mean(x@u), mean(x@l)))

setMethod("median", "DiscretePBox", function(x) c(median(x@u), median(x@l)))

# Williamson & Downs, p. 126-127
discrete_pbox_central_moment = function(pbox, k) {
    mu = mean(pbox@u@points)
    ml = mean(pbox@l@points)
    points1 = pbox@l@points - mu
    points2 = pbox@u@points - mu
    points3 = pbox@l@points - ml
    points4 = pbox@u@points - ml
    mml = pmin(points1, points2, points3, points4) ** k
    mmu = pmax(points1, points2, points3, points4) ** k
    c(mean(pmin(mml, mmu)), mean(pmax(mml, mmu)))
}

setMethod("sd", "DiscretePBox",
          function(x) sqrt(discrete_pbox_central_moment(x, 2)))

# unary negation
setMethod("-", signature(e1="DiscretePBox", e2="missing"),
          function(e1) discrete_pbox(-rev(e1@l@points), -rev(e1@u@points)))

setMethod("+", signature(e1="DiscretePBox", e2="numeric"),
          function(e1, e2) discrete_pbox(e1@u@points+e2, e1@l@points+e2))

setMethod("+", signature(e1="numeric", e2="DiscretePBox"),
          function(e1, e2) e2 + e1)

setMethod("-", signature(e1="DiscretePBox", e2="numeric"),
          function(e1, e2) e1 + (-e2))

setMethod("-", signature(e1="numeric", e2="DiscretePBox"),
          function(e1, e2) e1 + (-e2))

setMethod(">=", signature(e1="DiscretePBox", e2="numeric"),
          function(e1, e2) e1@u@points[1] >= e2)

setMethod(">", signature(e1="DiscretePBox", e2="numeric"),
          function(e1, e2) e1@u@points[1] > e2)

setMethod("<=", signature(e1="DiscretePBox", e2="numeric"),
          function(e1, e2) e1@l@points[length(e1)] <= e2)

setMethod("<", signature(e1="DiscretePBox", e2="numeric"),
          function(e1, e2) e1@l@points[length(e1)] < e2)

setMethod(
    "/", signature(e1="numeric", e2="DiscretePBox"),
    function(e1, e2) {
        stopifnot(e1 >= 0, e2 > 0)
        discrete_pbox(rev(e1/e2@l@points), rev(e1/e2@u@points))
    }
)

# apply function on every combination of elements of xs1 and xs2, and sort
sortedfunc = function(func, xs1, xs2) {
    ys1 = rep(xs1, length(xs2))
    ys2 = rep(xs2, rep(length(xs1), length(xs2)))
    sort(func(ys1, ys2))
}

# implementation of Williamson & Downs, Figure 14, page 127
discrete_pbox_convolution = function(func) function(pbox1, pbox2) {
    n = length(pbox1)
    stopifnot(n == length(pbox2))
    ixs = (0:(n-1)) * n
    discrete_pbox(
        upoints=sortedfunc(func, pbox1@u@points, pbox2@u@points)[ixs + 1],
        lpoints=sortedfunc(func, pbox1@l@points, pbox2@l@points)[ixs + n])
}

setMethod("+", "DiscretePBox", function(e1, e2) discrete_pbox_convolution(`+`)(e1, e2))

setMethod("*", "DiscretePBox", function(e1, e2) {
    stopifnot(e1 >= 0, e2 >= 0)
    discrete_pbox_convolution(`*`)(e1, e2)
})

setMethod("-", "DiscretePBox", function(e1, e2) e1 + (-e2))

setMethod("/", "DiscretePBox", function(e1, e2) e1 * (1 / e2))

# implementation of Williamson & Downs, page 120-121 (cases i and ii)
discrete_pbox_frechet = function(func) function(pbox1, pbox2) {
    n = length(pbox1)
    stopifnot(n == length(pbox2))
    discrete_pbox(
        upoints=sapply(1:n, function(i) max(func(pbox1@u@points[1:i], pbox2@u@points[i:1]))),
        lpoints=sapply(1:n, function(i) min(func(pbox1@l@points[i:n], pbox2@l@points[n:i]))))
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

setMethod("plot", "DiscretePBox", function(x) {
    n = length(x)
    cs = rep(2, n)
    uxs = c(rep(x@u@points, cs), x@l@points[length(x)])
    uys = rep((0:n) / n, c(1, cs))
    lxs = c(x@u@points[1], rep(x@l@points, cs))
    lys = rep((0:n) / n, c(cs, 1))
    plot(uxs, uys, col=2, type="l", xlab="x", ylab="cdf")
    lines(lxs, lys, col=1)
    legend("topleft",
           legend=c("upper cdf", "lower cdf"),
           col=c(2, 1),
           lty=1)    
})

test1 = function() {
    pb = discrete_pbox_norm(n=10)
    plot(pb)
    xs = seq(-5,5,0.01)
    lines(xs, sapply(xs, pnorm), type="l")
}

test2 = function() {
    pb = discrete_pbox_unif(min=-2, max=3, n=10)
    plot(pb)
    xs = seq(-4,4,0.01)
    lines(xs, sapply(xs, function(x) punif(x, min=-2, max=3)), type="l")
}

# W&D figure 19
test3 = function() {
    pb1 = discrete_pbox_unif(min=1, max=2, n=40)
    pb2 = discrete_pbox_unif(min=1, max=2, n=40)
    pb3 = pb1 + pb2
    plot(pb3)
}

# W&D figure 21
test4 = function() {
    pb1 = discrete_pbox_unif(min=0, max=1, n=40)
    pb2 = discrete_pbox_unif(min=0, max=1, n=40)
    pb3 = pb1 %fadd% pb2
    plot(pb3)
}

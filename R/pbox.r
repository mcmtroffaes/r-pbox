# a cdf sequence is a sequence of numbers x1 ... xn so that
# cdf(x)=0   for x<x1
# cdf(x)=1/n for x1<=x<x2
# cdf(x)=2/n for x2<=x<x3
# ...
# cdf(x)=1   for xn<=x

check_cdf = function(object) {
    if (length(object@xs) < 1L) "must contain at least one point"
    else if (is.unsorted(object@xs)) "points must be sorted"
    else TRUE
}

setClass("Cdf", representation(xs="numeric"),
         list(xs=0),
         validity=check_cdf)

setMethod("show", "Cdf", function(object) cat("points for cdf =", object@xs, "\n"))

cdf = function(xs) new("Cdf", xs=xs)

setMethod("length", "Cdf", function(x) length(x@xs))

setMethod("mean", "Cdf", function(x) mean(x@xs))

setMethod("median", "Cdf", function(x) median(x@xs))

# calculate the cdf as a function of x according the above
# discrete representation (note: implementation is simple but slow)
eval_cdf = function(cdf) function(x) {
    n = length(cdf)
    for (i in 1:n) { if (x < cdf@xs[i]) return ((i - 1) / n) }
    return(1)
}

# a p-box is represented by two cdf sequences
# u@xs[1] ... u@xs[n]
# l@xs[1] ... l@xs[n]

check_pbox = function(object) {
    if (length(object@u) != length(object@l)) "ucdf and lcdf must have same length"
    else if (any(object@u@xs > object@l@xs)) "ucdf must be left of lcdf"
    else TRUE
}

setClass("PBox", representation(u="Cdf", l="Cdf"),
         validity=check_pbox)

pbox = function(uxs, lxs) new("PBox", u=cdf(uxs), l=cdf(lxs))

setMethod("length", "PBox", function(x) length(x@u@xs))

setMethod("show", "PBox", function(object) {
    cat("points for upper cdf =", object@u@xs, "\n");
    cat("points for lower cdf =", object@l@xs, "\n")
})

# bot = bottom quantile (cut off at 1/n)
# top = top quantile (cut off at 1-1/n)
pbox_norm = function(mean=0, sd=1, bot=0.001, top=0.999, n=200) {
    stopifnot(bot > 0, top < 1)
    ps = c(min(bot, 1/n), seq(0, 1, length.out=n+1)[2:n], max(top, 1-1/n))
    xs = sapply(ps, function(x) qnorm(x, mean=mean, sd=sd))
    pbox(uxs=xs[-(n+1)], lxs=xs[-1])
}

pbox_unif = function(min=0, max=1, n=200) {
    ps = seq(0, 1, length.out=n+1)
    xs = sapply(ps, function(x) qunif(x, min=min, max=max))
    pbox(uxs=xs[-(n+1)], lxs=xs[-1])
}

eval_pbox_u = function(pbox) eval_cdf(pbox@u)

eval_pbox_l = function(pbox) eval_cdf(pbox@l)

setMethod("mean", "PBox", function(x) c(mean(x@u), mean(x@l)))

setMethod("median", "PBox", function(x) c(median(x@u), median(x@l)))

# Williamson & Downs, p. 126-127
pbox_central_moment = function(pbox, k) {
    mu = mean(pbox@u@xs)
    ml = mean(pbox@l@xs)
    xs1 = pbox@l@xs - mu
    xs2 = pbox@u@xs - mu
    xs3 = pbox@l@xs - ml
    xs4 = pbox@u@xs - ml
    mml = pmin(xs1, xs2, xs3, xs4) ** k
    mmu = pmax(xs1, xs2, xs3, xs4) ** k
    c(mean(pmin(mml, mmu)), mean(pmax(mml, mmu)))
}

setMethod("sd", "PBox", function(x) sqrt(pbox_central_moment(x, 2)))

# unary negation
setMethod("-", signature(e1="PBox", e2="missing"),
          function(e1) pbox(-rev(e1@l@xs), -rev(e1@u@xs)))

setMethod("+", signature(e1="PBox", e2="numeric"),
          function(e1, e2) pbox(e1@u@xs+e2, e1@l@xs+e2))

setMethod("+", signature(e1="numeric", e2="PBox"),
          function(e1, e2) e2 + e1)

setMethod("-", signature(e1="PBox", e2="numeric"),
          function(e1, e2) e1 + (-e2))

setMethod("-", signature(e1="numeric", e2="PBox"),
          function(e1, e2) e1 + (-e2))

setMethod(">=", signature(e1="PBox", e2="numeric"),
          function(e1, e2) e1@u@xs[1] >= e2)

setMethod(">", signature(e1="PBox", e2="numeric"),
          function(e1, e2) e1@u@xs[1] > e2)

setMethod("<=", signature(e1="PBox", e2="numeric"),
          function(e1, e2) e1@l@xs[length(e1)] <= e2)

setMethod("<", signature(e1="PBox", e2="numeric"),
          function(e1, e2) e1@l@xs[length(e1)] < e2)

setMethod(
    "/", signature(e1="numeric", e2="PBox"),
    function(e1, e2) {
        stopifnot(e1 >= 0, e2 > 0)
        pbox(rev(e1/e2@l@xs), rev(e1/e2@u@xs))
    }
)

# apply function on every combination of elements of xs1 and xs2, and sort
sortedfunc = function(func, xs1, xs2) {
    ys1 = rep(xs1, length(xs2))
    ys2 = rep(xs2, rep(length(xs1), length(xs2)))
    sort(func(ys1, ys2))
}

# implementation of Williamson & Downs, Figure 14, page 127
pbox_convolution = function(func) function(pbox1, pbox2) {
    n = length(pbox1)
    stopifnot(n == length(pbox2))
    ixs = (0:(n-1)) * n
    pbox(uxs=sortedfunc(func, pbox1@u@xs, pbox2@u@xs)[ixs + 1],
         lxs=sortedfunc(func, pbox1@l@xs, pbox2@l@xs)[ixs + n])
}

setMethod("+", "PBox", function(e1, e2) pbox_convolution(`+`)(e1, e2))

setMethod("*", "PBox", function(e1, e2) {
    stopifnot(e1 >= 0, e2 >= 0)
    pbox_convolution(`*`)(e1, e2)
})

setMethod("-", "PBox", function(e1, e2) e1 + (-e2))

setMethod("/", "PBox", function(e1, e2) e1 * (1 / e2))

# implementation of Williamson & Downs, page 120-121 (cases i and ii)
pbox_frechet = function(func) function(pbox1, pbox2) {
    n = length(pbox1)
    stopifnot(n == length(pbox2))
    pbox(uxs=sapply(1:n, function(i) max(func(pbox1@u@xs[1:i], pbox2@u@xs[i:1]))),
         lxs=sapply(1:n, function(i) min(func(pbox1@l@xs[i:n], pbox2@l@xs[n:i]))))
}

setGeneric("%fadd%", function(e1, e2) standardGeneric("%fadd%"))
setGeneric("%fmul%", function(e1, e2) standardGeneric("%fmul%"))
setGeneric("%fsub%", function(e1, e2) standardGeneric("%fsub%"))
setGeneric("%fdiv%", function(e1, e2) standardGeneric("%fdiv%"))

setMethod("%fadd%", "PBox", function(e1, e2) pbox_frechet(`+`)(e1, e2))

setMethod("%fmul%", "PBox", function(e1, e2) {
    stopifnot(e1 >= 0, e2 >= 0)
    pbox_frechet(`*`)(e1, e2)
})

setMethod("%fsub%", "PBox", function(e1, e2) e1 %fadd% (-e2))

setMethod("%fdiv%", "PBox", function(e1, e2) e1 %fmul% (1 / e2))

setMethod("plot", "PBox", function(x) {
    n = length(x)
    cs = rep(2, n)
    uxs = c(rep(x@u@xs, cs), x@l@xs[length(x)])
    uys = rep((0:n) / n, c(1, cs))
    lxs = c(x@u@xs[1], rep(x@l@xs, cs))
    lys = rep((0:n) / n, c(cs, 1))
    plot(uxs, uys, col=2, type="l", xlab="x", ylab="cdf")
    lines(lxs, lys, col=1)
    legend("topleft",
           legend=c("upper cdf", "lower cdf"),
           col=c(2, 1),
           lty=1)    
})

test1 = function() {
    pb = pbox_norm(n=10)
    plot(pb)
    xs = seq(-5,5,0.01)
    lines(xs, sapply(xs, pnorm), type="l")
}

test2 = function() {
    pb = pbox_unif(min=-2, max=3, n=10)
    plot(pb)
    xs = seq(-4,4,0.01)
    lines(xs, sapply(xs, function(x) punif(x, min=-2, max=3)), type="l")
}

# W&D figure 19
test3 = function() {
    pb1 = pbox_unif(min=1, max=2, n=40)
    pb2 = pbox_unif(min=1, max=2, n=40)
    pb3 = pb1 + pb2
    plot(pb3)
}

# W&D figure 21
test4 = function() {
    pb1 = pbox_unif(min=0, max=1, n=40)
    pb2 = pbox_unif(min=0, max=1, n=40)
    pb3 = pb1 %fadd% pb2
    plot(pb3)
}

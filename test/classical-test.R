
# Get p-values from classic tests
pillai_pval <- function(x, test) {
	# the following code is adapted from summary.manova
	eigs <- Re(eigen(qr.coef(qr(x$SSPE), x$SSPH), symmetric = FALSE)$values)
	stat <- test(eigs, x$df, x$df.residual)
	pval <- pf(stat[2], stat[3], stat[4], lower.tail = FALSE)
	return(pval)
}

library(car)
N <- 145
df <- data.frame(
    y1 = rnorm(N), 
    y2 = rnorm(N), 
	y3 = rnorm(N), 
	y4 = rnorm(N), 
    x = as.factor(sample(1:4, N, replace = T))
)
model = lm(cbind(y1, y2, y3, y4) ~ x + 0, df)
x = linearHypothesis(model, h = matrix(c(1, -1, 0, 0), 1, 4))
x

pillai_pval(x, stats:::Pillai)


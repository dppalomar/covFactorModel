png(file="cov_stat_perform.png",width = 10000,height = 8000,res = 1000)
par(mfcol = c(2, 2))
colors <- c("blue", "green4", "darkmagenta", "red3", "goldenrod4")

load("data/diagonal_Psi.RData")
matplot(index_T/N, PRIAL*100,
        xlab = "T/N", ylab = "PRIAL",
        main = "Data generated with diagonal Psi",
        type = "b", pch = 20, lwd = 2, col = colors)
legend("bottomleft", inset=0.01, legend = colnames(res), pch = 20, col = colors)

load("data/block_diagonal_Psi.RData")
matplot(index_T/N, PRIAL*100,
        xlab = "T/N", ylab = "PRIAL",
        main = "Data generated with block diagonal Psi",
        type = "b", pch = 20, lwd = 2, col = colors)
legend("bottomleft", inset=0.01, legend = colnames(res), pch = 20, col =colors)

load("data/scale_identity_Psi.RData")
matplot(index_T/N, PRIAL*100,
        xlab = "T/N", ylab = "PRIAL",
        main = "Data generated with scaled identity Psi",
        type = "b", pch = 20, lwd = 2, col = colors)
legend("bottomleft", inset=0.01, legend = colnames(res), pch = 20, col = colors)

load("data/full_Psi.RData")
matplot(index_T/N, PRIAL*100,
        xlab = "T/N", ylab = "PRIAL",
        main = "Data generated with full Psi",
        type = "b", pch = 20, lwd = 2, col = colors)
legend("bottomleft", inset=0.01, legend = colnames(res), pch = 20, col = colors)
dev.off()

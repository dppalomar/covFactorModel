png(file="package_compare.png",width = 10000,height = 4000,res = 1000)
par(mfcol = c(1, 2))

load("data/compare_covFactorModel_factanal.RData")
colors <- c("blue", "green4", "darkmagenta", "red3")

matplot(index_T/N, PRIAL,
        xlab = "T/N", ylab = "PRIAL",
        main = "",
        type = "b", pch = 20, lwd = 2, col = colors)
legend("center", inset=0.01, legend = colnames(PRIAL), pch = 20, col = colors)

matplot(index_T/N, time,
        xlab = "T/N", ylab = "Time",
        main = "",
        type = "b", pch = 20, lwd = 2, col = colors)
legend("center", inset=0.01, legend = colnames(time), pch = 20, col =colors)

dev.off()

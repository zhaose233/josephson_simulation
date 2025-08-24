d <- read.csv("out.csv")
dd <- read.csv("out_driving.csv")
# d$delta_phi <- d$phi - c(0,d$phi[1:length(d$phi) - 1])

plot(d$V, d$Ib, type = "l")

plot(dd$t, dd$V, type = "l")

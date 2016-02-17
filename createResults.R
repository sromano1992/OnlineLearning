names = c("nano-nano","disease-disease","drugs-drugs","nano-disease","nano-drugs","drugs-disease","nano-drugs-disease");
vi = c(0.07,0.24,0.78,0.49,1.01,0.15,0.18)
kt = c(0.8620,0.9815,0.9958,0.9965,0.9973,0.9978,0.9980)
chemical = data.frame(names,vi,kt)

names = c("nano-nano","disease-disease","drugs-drugs","chemical-chemical", "nano-disease","nano-drugs","drugs-disease","nano-drugs-disease","disease-chemical","drugs-chemical","drugs-disease-chemical","nano-chemical","nano-disease-chemical","nano-drugs-chemical")
vi = c(0,0,0,0.12,0,0.52,0,0,0.04,0.08,0.01,0.08,0.33,0.12)
kt = c(1,1,0.9968,0.9994,0.9999,0.9971,0.9983,0.9989,0.9995,0.9994,0.9995,0.9994,0.9994,0.9995)
nano_1 = data.frame(names,vi,kt)

vi = c(0,0,0.26,0.34,0.02,0.47,0,0,0.03,0.12,0.08,0.16,0.1,0.12)
kt = c(1,1,0.9963,0.9994,1,0.9956,0.9978,0.9982,0.9995,0.9995,0.9994,0.9994,0.9995,-1)
nano_2 = data.frame(names,vi,kt)

vi = c(0,0,0.26,1.34,0,0.23,0,0,0.09,0.16,0.02,1.01,0.06,0.27)
kt = c(1,1,0.9959,0.9993,1,0.9960,0.9979,0.9983,0.9993,0.9994,0.9995,0.9994,0.9995,0.9994)
nano_3 = data.frame(names,vi,kt)

par(mfrow=c(2,1))
plot(chemical$kt,chemical$vi,xlim = c(0.86,1),ylim=c(0,2.5), col=topo.colors(4)[1],pch=2,xlab = "Kendall Tau",ylab = "VI",lwd=2)
points(nano_1$kt, nano_1$vi, col=topo.colors(4)[2],pch=1,lwd=2)
points(nano_2$kt, nano_2$vi, col=topo.colors(4)[3],pch=1,lwd=2)
points(nano_3$kt, nano_3$vi, col=heat.colors(1),pch=1,lwd=2)
legend("topleft", ncol=1,inset=.01, 
        c("no-chemical","no-nano_1","no-nano_2","no-nano_3"), 
        col=c(topo.colors(4)[1:3],heat.colors(1)), 
        cex=0.8, pch = c(2,1,1,1),bty = "n",
        y.intersp = 0.4)
title("Results")

plot(chemical$kt,chemical$vi,xlim = c(0.995,1),ylim=c(0,2.5), col=topo.colors(4)[1],pch=2,xlab = "Kendall Tau",ylab = "VI",lwd=2)
points(nano_1$kt, nano_1$vi, col=topo.colors(4)[2],pch=1,lwd=2)
points(nano_2$kt, nano_2$vi, col=topo.colors(4)[3],pch=1,lwd=2)
points(nano_3$kt, nano_3$vi, col=heat.colors(1),pch=1,lwd=2)
legend("topleft", ncol=1,inset=.01, 
       c("no-chemical","no-nano_1","no-nano_2","no-nano_3"), 
       col=c(topo.colors(4)[1:3],heat.colors(1)), 
       cex=0.8, pch = c(2,1,1,1),bty = "n",
       y.intersp = 0.4)
title("Results")

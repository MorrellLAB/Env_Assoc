#lengend for figure3-Lei et al., 2018 Enviromental assocaiton
pdf(file = "/Users/lilei/Downloads/legend",width=16.00, height=8.67)
plot(0,0)
legend("topright",c("Ancestral allelic type","Derived allelic type","Wild barley range"),bty='n',lty = c(0,0,3),pch=c(1,20,NA),pt.cex = 1.2, cex = 1.5,col=c(adjustcolor("blue",alpha.f = 0.9),adjustcolor("deeppink1",alpha.f = 0.5),"black"))
dev.off()

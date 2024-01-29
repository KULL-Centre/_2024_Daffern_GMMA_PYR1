
plot_comp = function(dp, legend_pos="topright", ...) {
    plot(0, 0, col=0, xlim=c(1,nrow(dp)), ylim=c(0,1), xlab="Number of substitutions per variant", ylab="Fraction", ...)
    points(dp$n_mut, dp$dens_lib, pch=20, lwd=2, type="b", col=2)
    points(dp$n_mut, dp$frac_act, pch=20, lwd=2, type="b", col=4)
    legend(legend_pos,c("Library composition","Fraction active"),pch=20,col=c(2,4))
}


d2 = read.csv("dia_lib2/composition.csv")
jpeg("library_dia_lib2.jpg", width=8, height=6, units="cm", res=300, quality=95, pointsize=8)
plot_comp(d2, main="DIAZI library 2")
dev.off()

d2r = read.csv("dia_lib2_reref/composition.csv")
jpeg("library_dia_lib2_reref.jpg", width=8, height=6, units="cm", res=300, quality=95, pointsize=8)
plot_comp(d2r, main="DIAZI:D80E:Y81V library 2")
dev.off()


a1 = read.csv("aba_lib1/composition.csv")
a2 = read.csv("aba_lib2/composition.csv")
a3 = read.csv("aba_lib3/composition.csv")
aa = read.csv("aba_lib123/composition.csv")

jpeg("library_aba.jpg", width=16, height=12, units="cm", res=300, quality=95, pointsize=8)
par(mfrow=c(2,2), mar=c(5,4,4,2)+.1, bg="white")

plot_comp(a1, main="ABA library 1")
mtext("A", cex=1.6, font=2, line=1.5, at=0)
plot_comp(a2, main="ABA library 2")
mtext("B", cex=1.6, font=2, line=1.5, at=0)
plot_comp(a3, main="ABA library 3")
mtext("C", cex=1.6, font=2, line=1.5, at=0)
plot_comp(aa, main="ABA library 1+2+3")
mtext("D", cex=1.6, font=2, line=1.5, at=0)

dev.off()


d1 = read.csv("dia_lib1/composition.csv")
d1r = read.csv("dia_lib1_reref/composition.csv")
d3 = read.csv("dia_lib3/composition.csv")
d12 = read.csv("dia_lib12_reref/composition.csv")

jpeg("library_dia.jpg", width=16, height=12, units="cm", res=300, quality=95, pointsize=8)
par(mfrow=c(2,2), mar=c(5,4,4,2)+.1, bg="white")

plot_comp(d1, legend_pos="top", main="DIAZI library 1")
mtext("A", cex=1.6, font=2, line=1.5, at=0)

plot_comp(d1r, main="DIAZI:D80E:Y81V library 1")
mtext("B", cex=1.6, font=2, line=1.5, at=0)

plot_comp(d3, main="DIAZI library 3")
mtext("C", cex=1.6, font=2, line=1.5, at=0)

plot_comp(d12, main="DIAZI:D80E:Y81V library 1+2")
mtext("D", cex=1.6, font=2, line=1.5, at=0)


dev.off()

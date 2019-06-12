plot_barplots <- function(by_type_tables, SV_type, main = "", ymax = NA, legend = F){
    mins <- apply(by_type_tables[[SV_type]][-1,], 2, min)
    maxes <- apply(by_type_tables[[SV_type]][-1,], 2, max)
    bar_sizes <- apply(by_type_tables[[SV_type]][-1,], 2, median)
    ref_ind <- by_type_tables[[SV_type]][1,]

    ymax <- ifelse(is.na(ymax), max(c(unlist(maxes), unlist(ref_ind))), ymax)
    locations <- barplot(bar_sizes, col = rep(rep(c(asex_blue, sex_red), 10), each = 3),
                         ylim = c(0, ymax), main = main, xaxt = "n", cex.axis=2) # ylab = 'Number of called variants'
    # axis(1, tick=F, labels = , cex.axis=2)
    text(locations[seq(2, by = 3, length = 10)],
         par("usr")[3] - 50, pos = 1,
         xpd = TRUE, labels = timemas$labels, cex = 1.5)

    w <- 0.1
    for ( i in 1:length(mins) ){
        x <- locations[i]
        y_min <- mins[i]
        y_max <- maxes[i]
        lines(c(x, x), c(y_min, y_max), xpd=T, lwd = 1.5)
        lines(c(x - w, x + w), c(y_min, y_min), xpd=T, lwd = 1.5)
        lines(c(x - w, x + w), c(y_max, y_max), xpd=T, lwd = 1.5)
        points(x, ref_ind[i], pch = 20, cex = 1.5)
    }

    if (legend){
        legend('topleft', bty = 'n', pch = 20, 'reference individual')
    }
}
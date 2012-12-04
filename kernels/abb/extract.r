
source('overview.r')

plots.tmp = plots

plotsizes = c()
for (plot in plots.tmp) {
    p = psp[which(psp$PLOT == plot),]
    plotsize = unique(p$PLOTSIZE)
    if (length(plotsize) > 1) {
      plots = plots[-which(plots==plot)]
    } else {
      plotsizes = rbind(plotsizes, c(plot, plotsize))
    }
}
write.table(plotsizes, file='plotsizes.csv', sep=',', row.names=F)


for (plot in plots) {

    p = psp[which(psp$PLOT == plot),]

    p_fy = p[which(p$YEAR == min(p$YEAR)),]
    p_fy_trees = unique(p_fy$TREE)

    p = p[which(p$TREE %in% p_fy_trees),]

    d = data.frame(year=p$YEAR, species=as.vector(p$SPEC), dbh=p$DBH, stringsAsFactors=F)
    write.table(na.exclude(d), file=paste(plot, 'csv', sep='.'), sep=',', row.names=F)
}

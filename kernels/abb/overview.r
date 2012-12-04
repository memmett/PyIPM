
load('pspksclean.rdata')

# cycle through plots and print some info
pnames = unique(psp$PLOT)

plots = c()

for (name in pnames) {
  p    = psp[which(psp$PLOT == name), ]
  p_fy = p[which(p$YEAR == min(p$YEAR)), ]
  yrs  = unique(p$YEAR)

  p_fy_trees = unique(p_fy$TREE)

  nyrs = length(yrs)

  if (nyrs < 2) next

  if (('SW' %in% unique(p_fy$SPEC)) & ('AW' %in% unique(p_fy$SPEC))) {
    nsw = length(p_fy$DBH[which(p_fy$SPEC == 'SW')])
    naw = length(p_fy$DBH[which(p_fy$SPEC == 'AW')])
    mdbh = mean(p_fy$DBH[which(p_fy$SPEC == 'SW')], na.rm=T)

    if (nsw > 30 & naw > 10 & mdbh > 160.0) {
      print(c('PLOT', name))

      plots = c(plots, name)

      summary = c("year", "nsw ", "tnsw", "naw ", "tnaw", "mdbh")

      for (yr in yrs) {
        p_yr = p[which(p$YEAR == yr), ]
        nsw = length(na.exclude(p_yr$DBH[which(p_yr$SPEC == 'SW' & p_yr$TREE %in% p_fy_trees)]))
        mdbh = mean(na.exclude(p_yr$DBH[which(p_yr$SPEC == 'SW' & p_yr$TREE %in% p_fy_trees)]))
        nsw2 = length(na.exclude(p_yr$DBH[which(p_yr$SPEC == 'SW')]))
        naw = length(na.exclude(p_yr$DBH[which(p_yr$SPEC == 'AW' & p_yr$TREE %in% p_fy_trees)]))
        naw2 = length(na.exclude(p_yr$DBH[which(p_yr$SPEC == 'AW')]))
	summary = rbind(summary, c(yr, nsw, nsw2, naw, naw2, mdbh))
      }
      print(summary)
    }
  }
}

print(plots)






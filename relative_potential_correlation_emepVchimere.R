# -------------------------------------------------------------------
# CORREALTION BETWEEN RELATIVE POTENTIALS OF CHIMERE AND EMEP sHERPA
# -------------------------------------------------------------------

# Input for this script are the SHERPA output files:
# 1) sherpa_emep_PM25_perSNAP_perPrec.txt 
# 2) sherpa_chimere_PM25_perSNAP_perPrec.txt
# They contain the results of 2 source apportionment runs (EMEP and CHIMERE) for 150 cities
# identifying 200 contributions from 4 regions (city core, functional urban area, national, 
# International), 10 sectors and 5 precursors. For small cities the core and FUA are taken together.

# Output of this script:
# - correlation plots between CHIMERE and EMEP basecase concentrations 
#   ('basecase_concentration_EMEP_vs_CHIMERE_v1/2/3.png')
#
# - for each city a X-Y plot between relative potentials (concentration change of a area-sector-precursor 
#   emission devided by the basecae concentration) of EMEP and CHIMERE. This is done at different
#   aggregation levels:
#     * per Area 
#     * per sector (the 10 SNAP sectors are aggregated to 5: residential, industry, transport, agriculture, other)
#     * per area-sector combination
#     * per area-sector-precursor combination
#   All these correlations are presented as tables and as a plot
#   correlation/<cityname>_relative_potential_correlation.png
#
# - A Histogram of pearson correlation with a classification of the aggreement.
#   'correlations/_rPearson_histogram.png'
#
# - maps of Europe with for the 150 cities the various pearson correlations
#     * per Area 
#     * per sector
#     * per area-sector combination
#     * per area-sector-precursor combination

# clean up and load libraries
rm(list = ls())
library(ggplot2)
library(plyr)
library(gridExtra) # to add tables to graphs
library(ggmap)
# library(scales)
# library(scatterpie)

# working directory
wd <- paste0('D:/handover/Chimere_Emep_Comparison/') 
setwd(wd)

# read output of the SHERPA rusn with CHIMERE-SHERPA and EMEP-sHERPA
emep.file <- 'sherpa_emep_PM25_perSNAP_perPrec.txt'
emep.res <- read.table(emep.file, skip = 5, header = TRUE, sep = ';', quote = "")
chimere.file <- 'sherpa_chimere_PM25_perSNAP_perPrec.txt'
chimere.res <- read.table(chimere.file, skip = 5, header = TRUE, sep = ';', quote = "")

# --- create a dataframe with the emep and chimere data in different columns ---
# remove 'model' column and change headers, add 'emep_' or 'chimere_' before column names
emep.res <- emep.res[, 2:12]
names(emep.res)[6:11] <- paste0('emep_', names(emep.res)[6:11])
chimere.res <- chimere.res[,2:12]
names(chimere.res)[6:11] <- paste0('chimere_', names(chimere.res)[6:11])

# merge chimere and emep data on city, snap, area and precursor
res.df <- merge(emep.res, chimere.res)

# --------------------------------------------------
# Correlation between total base case concentrations
# --------------------------------------------------

total.conc.df <- ddply(res.df, c("target"), summarise,
                       emep_bc_conc = mean(emep_target_conc_basecase),
                       sd_emep = sd(emep_target_conc_basecase),
                       chimere_bc_conc = mean(chimere_target_conc_basecase))

# regression line
lm_bcc <- lm(data = total.conc.df, formula = chimere_bc_conc ~ emep_bc_conc)
summary(lm_bcc)

# x-y scatter plot of the total concentration predicted by EMEP and CHIMERE
# a simple R plot
png("basecase_concentration_EMEP_vs_CHIMERE_v1.png", width = 2*480, height = 2*480)
plot(total.conc.df$emep_bc_conc, total.conc.df$chimere_bc_conc,
     xlim = c(5, 35), ylim = c(5, 35), asp = 1,
     xlab = bquote("EMEP PM"[2.5]~"("*mu*"g/m"^3*")"),
     ylab = bquote("CHIMERE PM"[2.5]~"("*mu*"g/m"^3*")"))
text(total.conc.df$emep_bc_conc, total.conc.df$chimere_bc_conc, labels = total.conc.df$target, pos = 4)
abline(a=0, b=1, col="red")
abline(a=coefficients(lm_bcc)[1], b=coefficients(lm_bcc)[2], col="blue")
dev.off()

# same as the previous plot but with only the big differences labeled
png("basecase_concentration_EMEP_vs_CHIMERE_v2.png", width = 2*480, height = 2*480)
plot(total.conc.df$emep_bc_conc, total.conc.df$chimere_bc_conc,
     xlim = c(5, 35), ylim = c(5, 35), asp = 1,
     xlab = bquote("EMEP PM"[2.5]~"("*mu*"g/m"^3*")"),
     ylab = bquote("CHIMERE PM"[2.5]~"("*mu*"g/m"^3*")"))
big.diffs <- abs(total.conc.df$emep_bc_conc - total.conc.df$chimere_bc_conc) > 6
text(x=total.conc.df$emep_bc_conc[big.diffs], 
     y=total.conc.df$chimere_bc_conc[big.diffs], 
     labels = total.conc.df$target[big.diffs], pos = 2, pch = 19)
#     cex = abs(total.conc.df$emep_bc_conc - total.conc.df$chimere_bc_conc))
abline(a=0, b=1, col="red")
abline(a=coefficients(lm_bcc)[1], b=coefficients(lm_bcc)[2], col="blue")
dev.off()

# same plot with ggplot
tiff("basecase_concentration_EMEP_vs_CHIMERE_v3.tiff",
     width = 480, height = 480, units = "px", pointsize = 24, res = 72*1.5)
p <- ggplot(data = total.conc.df, aes(x=emep_bc_conc, y=chimere_bc_conc)) + geom_point()
p <- p + geom_abline(intercept = 0, slope = 1, color="red") + geom_smooth(method='lm', se = FALSE)
p <- p + xlab(bquote("EMEP PM"[2.5]~"("*mu*"g/m"^3*")")) + ylab(bquote("CHIMERE PM"[2.5]~"("*mu*"g/m"^3*")")) 
p <- p + xlim(4, 35) + ylim(4, 35)
# p <- p + ggtitle(bquote("Total modelled PM"[2.5]~"concentration"))
p
dev.off()

# bias, RMSE and correlation coefficient between CHIMERE and EMEP
bias.chimere.min.emep <- mean(total.conc.df$chimere_bc_conc - total.conc.df$emep_bc_conc)
# 3.790397
rmse.chimere.emep <- sqrt(mean((total.conc.df$chimere_bc_conc - total.conc.df$emep_bc_conc)^2))
# 4.659562
cor.chimere.emep <- cor(total.conc.df$emep_bc_conc, total.conc.df$chimere_bc_conc)
# 0.833478

# re-label the sectors from 10 SNAPs to 5
res.df <- data.frame(res.df, sector = NA, area = NA)
res.df$sector[res.df$snap %in% c(2)] <- 'Residential'
res.df$sector[res.df$snap %in% c(1, 3, 4)] <- 'Industry'
res.df$sector[res.df$snap %in% c(7)] <- 'Road Transport'
res.df$sector[res.df$snap %in% c(10)] <- 'Agriculture'
res.df$sector[res.df$snap %in% c(5, 6, 8, 9)] <- 'Other'
for (i in 1:NROW(res.df)) {
  res.df$area[i] <- unlist(strsplit(toString(res.df$source[i]), '_'))[2]
}

# aggregation levels
aggregation.levels <- c('Area-Sector-Precursor', 'Area-Sector', 'Sector', 'Area', 'Total', 'Base case')

# different aggregations of relative potency
# base case concentration for each city
rp.bc.df <- ddply(res.df, c('target'), summarise,
                   aggregation = 'Base case',
                   emep_rp = 100,
                   chimere_rp = 100,
                   emep_dc = mean(emep_target_conc_basecase),
                   chimere_dc = mean(chimere_target_conc_basecase))
# total DC
rp.tot.df <- ddply(res.df, c('target'), summarise,
                   aggregation = 'Total',
                   emep_rp = sum(emep_relative_potential),
                   chimere_rp = sum(chimere_relative_potential),
                   emep_dc = sum(emep_delta_conc),
                   chimere_dc = sum(chimere_delta_conc))
# per city and area (4)
rp.area.df <- ddply(res.df, c('target', 'area'), summarise,
                    aggregation = 'Area',
                    emep_rp = sum(emep_relative_potential),
                    chimere_rp = sum(chimere_relative_potential),
                    emep_dc = sum(emep_delta_conc),
                    chimere_dc = sum(chimere_delta_conc))
# per city and sector
rp.sector.df <- ddply(res.df, c('target', 'sector'), summarise,
                      aggregation = 'Sector',
                      emep_rp = sum(emep_relative_potential),
                      chimere_rp = sum(chimere_relative_potential),
                      emep_dc = sum(emep_delta_conc),
                      chimere_dc = sum(chimere_delta_conc))
# per city, sector and area
rp.SectorArea.df <- ddply(res.df, c('target', 'sector', 'area'), summarise,
                          aggregation = 'Area-Sector',
                          emep_rp = sum(emep_relative_potential),
                          chimere_rp = sum(chimere_relative_potential),
                          emep_dc = sum(emep_delta_conc),
                          chimere_dc = sum(chimere_delta_conc))
rp.SAP.df <- ddply(res.df, c('target', 'sector', 'area', 'precursor'), summarise,
                   aggregation = 'Area-Sector-Precursor',
                   emep_rp = sum(emep_relative_potential),
                   chimere_rp = sum(chimere_relative_potential),
                   emep_dc = sum(emep_delta_conc),
                   chimere_dc = sum(chimere_delta_conc))

# put aggregations together to plot them all at once
rp.agg.df <- rbind(rp.bc.df[, names(rp.tot.df)],
                   rp.tot.df[, names(rp.tot.df)],
                    rp.area.df[, names(rp.tot.df)],
                    rp.sector.df[, names(rp.tot.df)],
                    rp.SectorArea.df[, names(rp.tot.df)],
                    rp.SAP.df[, names(rp.tot.df)])
rp.agg.df$aggregation <- factor(rp.agg.df$aggregation, 
                                 levels = aggregation.levels,
                                 ordered = TRUE)
rp.agg.df <- rp.agg.df[order(rp.agg.df$target, rp.agg.df$aggregation),]

# calculate correlations at all aggregation levels
rp.corr.df <- ddply(rp.agg.df[rp.agg.df$aggregation != 'Total' & rp.agg.df$target != 'Lefkosia',], 
                    c('target', 'aggregation'), summarize,
                    rmse.rp = sqrt(mean((emep_rp - chimere_rp)^2)),
                    r.pearson.rp = cor(emep_rp, chimere_rp, method = 'pearson'),
                    spearman.rp = cor(emep_rp, chimere_rp, method = 'spearman'))
rp.corr.df$aggregation <- factor(rp.corr.df$aggregation, 
                                 levels = aggregation.levels,
                                 ordered = TRUE)
rp.corr.df <- rp.corr.df[order(rp.corr.df$target, rp.corr.df$aggregation),]

# Are these indicators very different: RMSE, Pearson, Spearman correlation?

# correlation between correlations
plot(x=rp.corr.df$rmse.rp[rp.corr.df$aggregation == "Area-Sector-Precursor"], 
     y=rp.corr.df$r.pearson.rp[rp.corr.df$aggregation == "Area-Sector-Precursor"],
     xlab = "RMSE of relative potentials", 
     ylab = "Pearson correlation of relative potentials")

plot(x=rp.corr.df$rmse.rp[rp.corr.df$aggregation == "Area-Sector-Precursor"], 
     y=rp.corr.df$spearman.rp[rp.corr.df$aggregation == "Area-Sector-Precursor"],
     xlab = "RMSE of relative potentials", 
     ylab = "Spearman correlation of relative potentials")

plot(x=rp.corr.df$r.pearson.rp[rp.corr.df$aggregation == "Area-Sector-Precursor"], 
     y=rp.corr.df$spearman.rp[rp.corr.df$aggregation == "Area-Sector-Precursor"],
     xlab = "Pearson of relative potentials", 
     ylab = "Spearman correlation of relative potentials")

# plot(rp.corr.df$spearman.rp, rp.corr.df$pearson.rp)


# ---------------------------------
#  Correlation plots for each city
# ---------------------------------

# directories for output plots
if (!dir.exists('correlations')) {dir.create('correlations')}

city.list <- unique(rp.tot.df$target)
city.list <- city.list[city.list != 'Lefkosia']
for (city in city.list) {
  rp.agg.city.df <- rp.agg.df[rp.agg.df$target == city, ]
  
  corr.table <- rp.corr.df[rp.corr.df$target==city & rp.corr.df$aggregation != 'Base case', c('aggregation', 'r.pearson.rp')]
  corr.table$r.pearson.rp <- round(corr.table$r.pearson.rp,2)
  
  # plot the correlation between relative potentials
  p <- ggplot(data = rp.agg.city.df, 
              aes(x = emep_rp, y = chimere_rp, col = aggregation)) + geom_point() + geom_smooth(method="lm", fill=NA)
  p <- p + geom_abline(intercept = 0, slope = 1) + geom_smooth(method="lm", fill=NA)
  p <- p + ggtitle(city) + xlab('EMEP relative potential') + ylab('CHIMERE relative potential')
  p <- p + annotation_custom(tableGrob(corr.table, rows=NULL), xmin=60, xmax=75, ymin=0, ymax=25)
  p <- p + xlim(0, 100) + ylim(0, 100)
  png(paste0('correlations/', city, '_7_relative_potential_correlation.png'), height = 480, width = 480*1.25)
  print(p)
  dev.off()
  
  # without the total and base case
  p <- ggplot(data = rp.agg.city.df[rp.agg.city.df$aggregation == "Area-Sector-Precursor",], 
              aes(x = log(emep_rp), y = log(chimere_rp), col = aggregation)) + geom_point() + geom_smooth(method="lm", fill=NA)
  p <- p + geom_abline(intercept = 0, slope = 1) #+ geom_abline(intercept = 0, slope = 0.9) + geom_abline(intercept = 0, slope = 1.1) 
  # p <- p + geom_smooth(method="lm", fill=NA)
  p <- p + ggtitle(city) + xlab('EMEP relative potential') + ylab('CHIMERE relative potential')
  # p <- p + annotation_custom(tableGrob(corr.table, rows=NULL), xmin=60, xmax=75, ymin=0, ymax=25)
  # p <- p + xlim(0, 4) + ylim(0, 4)
  # png(paste0('correlations/', city, '_8_relative_potential_correlation_with_error.png'), height = 480, width = 480*1.25)
  print(p)
  # dev.off()
  
  # -------------------------
  # table with ranked relative potentials per area-sector-precursor combination
  rp.SAP.city.df <- rp.SAP.df[rp.SAP.df$target == city,]
  emep.rank <- - rank(rp.SAP.city.df$emep_rp) + NROW(rp.SAP.city.df) + 1
  chimere.rank <- - rank(rp.SAP.city.df$chimere_rp) + NROW(rp.SAP.city.df) + 1
  rp.SAP.city.df <- data.frame(rp.SAP.city.df, emep.rank, chimere.rank)
  rp.SAP.city.table <- rp.SAP.city.df[order(emep.rank + chimere.rank), 
                                      c("area", "sector", "precursor", "emep_rp", "emep.rank", "chimere_rp", "chimere.rank")]
  rp.SAP.city.table$emep_rp <- round(rp.SAP.city.table$emep_rp, digits=1)
  rp.SAP.city.table$chimere_rp <- round(rp.SAP.city.table$chimere_rp, digits=1)
  rp.SAP.city.table$emep.rank <- round(rp.SAP.city.table$emep.rank, digits=0)
  rp.SAP.city.table$chimere.rank <- round(rp.SAP.city.table$chimere.rank, digits=0)
  top10 <- rp.SAP.city.table$emep.rank <= 10 | rp.SAP.city.table$chimere.rank <= 10
  
  png(paste0('correlations/', city, '_5_area_sector_precursor_table.png'),
      width = 600, height = 23 * (NROW(rp.SAP.city.table) + 1))
  grid.table(rp.SAP.city.table, row = NULL)
  dev.off()

  png(paste0('correlations/', city, '_6_area_sector_precursor_table_top10.png'),
      width = 600, height = 23 * (sum(top10) + 1))
  grid.table(rp.SAP.city.table[top10,], row = NULL)
  dev.off()

  # -------------------------
  # table with ranked relative potentials per area-sector combination
  rp.AS.city.df <- rp.SectorArea.df[rp.SectorArea.df$target == city,]
  emep.rank <- - rank(rp.AS.city.df$emep_rp) + NROW(rp.AS.city.df) + 1
  chimere.rank <- - rank(rp.AS.city.df$chimere_rp) + NROW(rp.AS.city.df) + 1
  rp.AS.city.df <- data.frame(rp.AS.city.df, emep.rank, chimere.rank)
  rp.AS.city.table <- rp.AS.city.df[order(emep.rank + chimere.rank), 
                                      c("area", "sector", "emep_rp", "emep.rank", "chimere_rp", "chimere.rank")]
  rp.AS.city.table$emep_rp <- round(rp.AS.city.table$emep_rp, digits=1)
  rp.AS.city.table$chimere_rp <- round(rp.AS.city.table$chimere_rp, digits=1)
  rp.AS.city.table$emep.rank <- round(rp.AS.city.table$emep.rank, digits=0)
  rp.AS.city.table$chimere.rank <- round(rp.AS.city.table$chimere.rank, digits=0)
  top10 <- rp.AS.city.table$emep.rank <= 10 | rp.AS.city.table$chimere.rank <= 10
  
  png(paste0('correlations/', city, '_3_area_sector_table.png'),
      width = 600, height = 23 * (NROW(rp.AS.city.table) + 1))
  grid.table(rp.AS.city.table, row = NULL)
  dev.off()
  
  png(paste0('correlations/', city, '_4_area_sector_table_top10.png'),
      width = 600, height = 23 * (sum(top10) + 1))
  grid.table(rp.AS.city.table[top10,], row = NULL)
  dev.off()

  # -------------------------
  # table with ranked relative potentials per area
  rp.area.city.df <- rp.area.df[rp.area.df$target == city,]
  emep.rank <- - rank(rp.area.city.df$emep_rp) + NROW(rp.area.city.df) + 1
  chimere.rank <- - rank(rp.area.city.df$chimere_rp) + NROW(rp.area.city.df) + 1
  rp.area.city.df <- data.frame(rp.area.city.df, emep.rank, chimere.rank)
  rp.area.city.table <- rp.area.city.df[order(emep.rank + chimere.rank), 
                                    c("area", "emep_rp", "emep.rank", "chimere_rp", "chimere.rank")]
  rp.area.city.table$emep_rp <- round(rp.area.city.table$emep_rp, digits=1)
  rp.area.city.table$chimere_rp <- round(rp.area.city.table$chimere_rp, digits=1)
  rp.area.city.table$emep.rank <- round(rp.area.city.table$emep.rank, digits=0)
  rp.area.city.table$chimere.rank <- round(rp.area.city.table$chimere.rank, digits=0)
  
  png(paste0('correlations/', city, '_1_area_table.png'),
      width = 400, height = 23 * (NROW(rp.area.city.table) + 1))
  grid.table(rp.area.city.table, row = NULL)
  dev.off()

  # -------------------------
  # table with ranked relative potentials per sector
  rp.sector.city.df <- rp.sector.df[rp.sector.df$target == city,]
  emep.rank <- - rank(rp.sector.city.df$emep_rp) + NROW(rp.sector.city.df) + 1
  chimere.rank <- - rank(rp.sector.city.df$chimere_rp) + NROW(rp.sector.city.df) + 1
  rp.sector.city.df <- data.frame(rp.sector.city.df, emep.rank, chimere.rank)
  rp.sector.city.table <- rp.sector.city.df[order(emep.rank + chimere.rank), 
                                        c("sector", "emep_rp", "emep.rank", "chimere_rp", "chimere.rank")]
  rp.sector.city.table$emep_rp <- round(rp.sector.city.table$emep_rp, digits=1)
  rp.sector.city.table$chimere_rp <- round(rp.sector.city.table$chimere_rp, digits=1)
  rp.sector.city.table$emep.rank <- round(rp.sector.city.table$emep.rank, digits=0)
  rp.sector.city.table$chimere.rank <- round(rp.sector.city.table$chimere.rank, digits=0)
  
  png(paste0('correlations/', city, '_2_sector_table.png'),
      width = 400, height = 23 * (NROW(rp.sector.city.table) + 1))
  grid.table(rp.sector.city.table, row = NULL)
  dev.off()
  

  # plot the correlation between delta concentrations
  p <- ggplot(data = rp.agg.city.df, aes(x = emep_dc, y = chimere_dc, col = aggregation)) + geom_point() + geom_smooth(method="lm", fill=NA)
  p <- p + geom_abline(intercept = 0, slope = 1) + geom_smooth(method="lm", fill=NA)
  p <- p + ggtitle(city) + xlab('EMEP delta concentration') + ylab('CHIMERE delta concentration')
  png(paste0('correlations/', city, '_8_delta_concentration_correlation.png'), height = 480, width = 480*1.25)
  print(p)
  dev.off()
  
}

# --------------------------------
# histogram of pearson correlation
# --------------------------------

score.colors <- c('firebrick4', 'firebrick1', 'gold', 'darkolivegreen2', 'darkolivegreen4')
score.breaks <- c(0.3, 0.7, 0.85, 0.9, 0.95, 1)
png('correlations/_rPearson_histogram.png')
pearson.hist <- hist(rp.corr.df$r.pearson.rp[rp.corr.df$aggregation == 'Area-Sector-Precursor'],
                     freq = FALSE,
                     breaks = score.breaks,
                     col = score.colors,
                     main = 'Pearson correlatoin between relative potentials for 150 cities',
                     xlab = "Pearson correlation")
percentages <- paste0(round(pearson.hist$density*c(0.4, 0.15, 0.05, 0.05, 0.05) * 100), '%')
text(x = (score.breaks[1:5] + score.breaks[2:6])/2, y = 0.5, labels = percentages)
legend("topleft", c('Very good', 'good', 'fair', 'bad', 'very bad'),
       pch = 15, col = rev(score.colors))
dev.off()

write.table(rp.corr.df, 'RP_corrs.txt', row.names = FALSE, quote = FALSE, sep = '\t')


# ------------------------------------------
# maps of correlations per aggregation level
# ------------------------------------------

city.file <- "city_list_fua150.txt"
city.df <- read.table(city.file, header = TRUE, sep = ';', quote = '')
map.boundaries <- c(min(city.df$lon), min(city.df$lat), max(city.df$lon), max(city.df$lat))
map.boundaries <- c(left = -10, bottom = 33, right = 34, top = 61)

# get an EU map from google
eumap <- get_map(source = "google", maptype = "roadmap", location = map.boundaries)

names(rp.corr.df)[1] <- 'cityname'
rp.corr.df <- merge(rp.corr.df, city.df[, c('cityname', 'lat', 'lon')])

# colours
score.colors <- c('firebrick4', 'firebrick1', 'gold', 'darkolivegreen2', 'darkolivegreen4')
score.breaks <- c(0.3, 0.7, 0.85, 0.9, 0.95, 1)

# create an output folder for the maps
if(!dir.exists("maps")) {dir.create("maps")}

# loop over aggregation levels
for (aggregation in aggregation.levels[!(aggregation.levels %in% c('Base case', 'Total'))]) {
  # select the data with the desired aggregation
  map.data <- rp.corr.df[rp.corr.df$aggregation == aggregation, ]
  # corr.quantiles <- as.numeric(quantile(map.data$rmse.rp, seq(0, 1, 0.25)) - min(map.data$rmse.rp)) / (max(map.data$rmse.rp) - min(map.data$rmse.rp))

  m1 <- ggmap(eumap) + ggtitle(paste("Pearson correlation between Emep and Chimere relative potential\naggregation:", aggregation))
  m1 <- m1 + geom_point(aes(x = lon, y = lat, colour = r.pearson.rp), 
                        data = rp.corr.df[rp.corr.df$aggregation == aggregation,],
                        size = 6)
  m1 <- m1 + scale_colour_gradientn(colours = score.colors, 
                                    values = rescale((score.breaks[1:5]+score.breaks[2:6])/2, to = c(0, 1))) # 
  m1 <- m1 + theme(text = element_text(size=20))
  png(file.path("maps", paste0('rPearson_Emep_Chimere_RP_', aggregation, '.png')),
      width = 2*480, height = 2*480)
  print(m1)
  dev.off()
}



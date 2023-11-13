# Analysis to compare the efficacy of ant repellents.
library(survival);library(survminer);library(grid);library(ggplot2);library(graphics)
setwd("/Users/tkay/Desktop/Work/Repellent")


# Read in and organize data
LasDF <- read.csv("Lasius.csv", header = TRUE, check.names = FALSE)
LasDF <- LasDF[,order(colMeans(LasDF))]
LasDF_FORCOL <- LasDF[colnames(LasDF) %in% c("Control", "2-HydroxyBA", "2-BromoBA", "2-EthoxyBA", "2-MethoxyBA")]
LasDF_trimmed <- LasDF[colnames(LasDF)%in%colnames(WasDF)]
LasDF_trimmed <- LasDF_trimmed[c(5,3,2,4,1)]
LasCol <- rainbow(5, alpha = 0.4)
LasCol_All <- rainbow(18, alpha = 0.4)
WasDF <- read.csv("Wasmannia.csv", header = TRUE, check.names = FALSE)
WasCol <- LasCol[colnames(LasDF_FORCOL) %in% colnames(WasDF)]
ParDF <- read.csv("Paratrechina.csv", header = TRUE, check.names = FALSE)
ParCol <- LasCol[colnames(LasDF_FORCOL) %in% colnames(ParDF)]
SolDF <- read.csv("Solenopsis.csv", header = TRUE, check.names = FALSE)
SolCol <- LasCol[colnames(LasDF_FORCOL) %in% colnames(SolDF)]
TapDF <- read.csv("Tapinoma.csv", header = TRUE, check.names = FALSE)
TapCol <- LasCol[colnames(LasDF_FORCOL) %in% colnames(TapDF)]

# Lasius boxplot (Fig. 1)
jpeg('Las_boxplot.jpg', width=6000, height=3000, unit='px')
par(mfrow=c(1,1), mar = c(1,1,1,1), mai = c(12,5,1,1), mgp = c(5, 5, 0), family = 'serif', bty="n")
boxplot(LasDF, col="darkgray", las=2, cex=3, ylab ="Time (seconds)", yaxt = "n", xaxt = "n",
        cex.axis = 5, cex.lab = 10, boxlwd = 8, lwd = 5, medlwd = 15)
axis(2, at = c(0,1200), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
axis(1, at = c(1:18), labels = c("citrus", "control", "3-hydroxyBA", "4-hydroxyBA", "hexanoic acid", "cis-2-hexenol", "AITC",
                                 "2-phenoxyBA", "trans-2-hexenol", "2-butoxyBA", "2-methylBA", "2-propoxyBA", "2-aminoBA",
                                 "BA", "2-hydroxyBA", "2-bromoBA", "2-ethoxyBA", "2-methoxyBA"),
     tick = TRUE, lwd = 8, cex.axis = 10, las = 2)
points(1:ncol(LasDF), colMeans(LasDF), col="red", pch=18, cex=12)
dev.off()

# Comparative boxplot (Fig. 2)
jpeg('Comparative_boxplot.jpg', width=3000, height=4000, unit='px')
par(mfrow=c(3,2), mar = c(2,2,2,2), mai = c(3,5,3,3), mgp = c(15, 5, 0), family = 'serif', bty="n")
boxplot(LasDF_trimmed, col=WasCol, las=2, cex=3, xaxt = "n", ylab ="Time (seconds)", yaxt = "n",
        cex.axis = 5, cex.lab = 10, boxlwd = 8, lwd = 10, medlwd = 20, main = expression(italic("Lasius niger")), cex.main = 10)
axis(2, at = c(0,1200), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
boxplot(WasDF, col=WasCol, las=2, cex=3, xaxt = "n", ylab ="Time (seconds)", yaxt = "n",
        cex.axis = 5, cex.lab = 10, boxlwd = 8, lwd = 10, medlwd = 20, main = expression(italic("Wasmannia auropunctata")), cex.main = 10)
axis(2, at = c(0,8000), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
boxplot(ParDF, col=ParCol, las=2, cex=3, xaxt = "n", ylab ="Time (seconds)", yaxt = "n",
        cex.axis = 5, cex.lab = 10, boxlwd = 8, lwd = 10, medlwd = 20, main = expression(italic("Paratrechina longicornis")), cex.main = 10)
axis(2, at = c(0,8000), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
boxplot(SolDF, col=SolCol, las=2, cex=3, xaxt = "n", ylab ="Time (seconds)", yaxt = "n",
        cex.axis = 5, cex.lab = 10, boxlwd = 8, lwd = 10, medlwd = 20, main = expression(italic("Solenopsis invicta")), cex.main = 10)
axis(2, at = c(0,1200), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
boxplot(TapDF, col=SolCol, las=2, cex=3, xaxt = "n", ylab ="Time (seconds)", yaxt = "n",
        cex.axis = 5, cex.lab = 10, boxlwd = 8, lwd = 10, medlwd = 20, main = expression(italic("Tapinoma magnum")), cex.main = 10)
axis(2, at = c(0,1200), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend(x = -0.1, y= 1, legend = colnames(WasDF),
       fill = WasCol, col = WasCol, cex = 12, lwd = 5, bty = "n")
dev.off()

# Survival analysis of all Lasius data
LasSurvDF <- data.frame(time = as.vector(as.matrix(LasDF)), chemical = rep(names(LasDF), each = nrow(LasDF)))
LasSurvDF$status <- LasSurvDF$time
LasSurvDF$status[LasSurvDF$status<1200] <- 2
LasSurvDF$status[LasSurvDF$status == 1200] <- 1
LasSurvDF$chemical <- factor(LasSurvDF$chemical, levels = rev(colnames(LasDF)))

LasFit <- survfit(Surv(time, status) ~ chemical, data = LasSurvDF)
LasFit.coxph <- coxph(Surv(time, status) ~ chemical, data = LasSurvDF)

jpeg('LasSurvCurves.jpg', width=6000, height=3000, unit='px')
ggsurvplot(LasFit,
           ggtheme = theme_classic2(base_size=100, base_family = "serif"),
           palette = rev(LasCol_All),
           size = 10,
           surv.scale = "percent",
           censor = FALSE) + 
  ylab("Percentage within perimeter") +
  xlab("Time (s)")
dev.off()
jpeg('LasSurvStats.jpg', width=3000, height=1500, unit='px')
a <- ggforest(LasFit.coxph, data = LasSurvDF,
              noDigits = 2,
              cpositions=c(0.01, 0.08, 0.25),
              fontsize = 3,
              main = "")
a + theme(plot.margin = unit(c(2,2,2,2), "cm"), element_line(size = 10))
dev.off()

# Contrasting only HydroxyBA with the negative control, Lasius
LasSurvDF_H <- LasSurvDF[LasSurvDF$chemical == "Control" | LasSurvDF$chemical == "2-HydroxyBA",]
LasSurvDF_H$chemical <- as.factor(as.character(LasSurvDF_H$chemical))
LasFit_H <- survfit(Surv(time, status) ~ chemical, data = LasSurvDF_H)
LasFit_H.coxph <- coxph(Surv(time, status) ~ chemical, data = LasSurvDF_H)

# Moving the carbon atom to positions 3 and 4
LasSurvDF_34 <- LasSurvDF[LasSurvDF$chemical == "Control" | LasSurvDF$chemical == "3-HydroxyBA" | LasSurvDF$chemical == "4-HydroxyBA",]
LasSurvDF_34$chemical <- as.factor(as.character(LasSurvDF_34$chemical))
LasSurvDF_34$chemical <- factor(LasSurvDF_34$chemical , levels = c("Control","3-HydroxyBA","4-HydroxyBA"))
LasFit_34 <- survfit(Surv(time, status) ~ chemical, data = LasSurvDF_34)
LasFit_34.coxph <- coxph(Surv(time, status) ~ chemical, data = LasSurvDF_34)

# Contrasting alternatives to the -hydroxy group
LasSurvDF_small <- LasSurvDF[LasSurvDF$chemical == "2-MethoxyBA" | LasSurvDF$chemical == "2-HydroxyBA" | LasSurvDF$chemical == "2-BromoBA" | LasSurvDF$chemical == "2-MethylBA" | LasSurvDF$chemical == "2-AminoBA" | LasSurvDF$chemical == "BA",]
LasSurvDF_small$chemical <- as.factor(as.character(LasSurvDF_small$chemical))
LasSurvDF_small$chemical <- factor(LasSurvDF_small$chemical , levels = c("2-MethoxyBA", "2-HydroxyBA","2-BromoBA","2-MethylBA", "2-AminoBA", "BA"))
LasFit_small <- survfit(Surv(time, status) ~ chemical, data = LasSurvDF_small)
LasFit_small.coxph <- coxph(Surv(time, status) ~ chemical, data = LasSurvDF_small)

# Varying carbon chain length
LasSurvDF_cc <- LasSurvDF[LasSurvDF$chemical == "2-MethoxyBA" | LasSurvDF$chemical == "2-EthoxyBA" | LasSurvDF$chemical == "2-PropoxyBA" | LasSurvDF$chemical == "2-ButoxyBA" | LasSurvDF$chemical == "2-PhenoxyBA",]
LasSurvDF_cc$chemical <- as.factor(as.character(LasSurvDF_cc$chemical))
LasSurvDF_cc$chemical <- factor(LasSurvDF_cc$chemical , levels = c("2-MethoxyBA", "2-EthoxyBA","2-PropoxyBA","2-ButoxyBA", "2-PhenoxyBA"))
LasFit_cc <- survfit(Surv(time, status) ~ chemical, data = LasSurvDF_cc)
LasFit_cc.coxph <- coxph(Surv(time, status) ~ chemical, data = LasSurvDF_cc)

# More contrasting chemical structures
LasSurvDF_other <- LasSurvDF[LasSurvDF$chemical == "Control" | LasSurvDF$chemical == "Cis-2-hexenol" | LasSurvDF$chemical == "Trans-2-hexenol" | LasSurvDF$chemical == "Hexanoic acid",]
LasSurvDF_other$chemical <- as.factor(as.character(LasSurvDF_other$chemical))
LasSurvDF_other$chemical <- factor(LasSurvDF_other$chemical , levels = c("Control","Cis-2-hexenol", "Trans-2-hexenol", "Hexanoic acid"))
LasFit_other <- survfit(Surv(time, status) ~ chemical, data = LasSurvDF_other)
LasFit_other.coxph <- coxph(Surv(time, status) ~ chemical, data = LasSurvDF_other)

# Other repellents
LasSurvDF_diff <- LasSurvDF[LasSurvDF$chemical == "Control" | LasSurvDF$chemical == "AITC" | LasSurvDF$chemical == "Citrus",]
LasSurvDF_diff$chemical <- as.factor(as.character(LasSurvDF_diff$chemical))
LasSurvDF_diff$chemical <- factor(LasSurvDF_diff$chemical , levels = c("Control","AITC","Citrus"))
LasFit_diff <- survfit(Surv(time, status) ~ chemical, data = LasSurvDF_diff)
LasFit_diff.coxph <- coxph(Surv(time, status) ~ chemical, data = LasSurvDF_diff)


# Comparative section
##########################################

# Wasmannia
WasSurvDF <- data.frame(time = as.vector(as.matrix(WasDF)), chemical = rep(names(WasDF), each = nrow(WasDF)))
WasSurvDF$status <- WasSurvDF$time
WasSurvDF$status[WasSurvDF$status<9000] <- 2
WasSurvDF$status[WasSurvDF$status == 9000] <- 1

WasSurvDF_cont <- WasSurvDF
WasSurvDF_cont$chemical <- as.factor(as.character(WasSurvDF_cont$chemical))
WasSurvDF_cont$chemical <- factor(WasSurvDF_cont$chemical , levels = c("Control","2-HydroxyBA","2-MethoxyBA","2-BromoBA", "2-EthoxyBA"))
WasFit_cont <- survfit(Surv(time, status) ~ chemical, data = WasSurvDF_cont)
WasFit_cont.coxph <- coxph(Surv(time, status) ~ chemical, data = WasSurvDF_cont)
jpeg('WasSurvStats.jpg', width=3000, height=1500, unit='px')
a <- ggforest(WasFit_cont.coxph, data = WasSurvDF_cont,
              noDigits = 2,
              cpositions=c(0.01, 0.08, 0.25),
              fontsize = 3,
              main = "")
a + theme(plot.margin = unit(c(2,2,2,2), "cm"), element_line(size = 10))
dev.off()
jpeg('WasSurvCurves.jpg', width=6000, height=3000, unit='px')
WAS_SC <- ggsurvplot(WasFit_cont,
                     ggtheme = theme_classic2(base_size=100, base_family = "serif"),
                     palette = rev(WasCol),
                     size = 10,
                     surv.scale = "percent",
                     censor = FALSE) + 
  ylab("Percentage within perimeter") +
  xlab("Time (s)")
WAS_SC
dev.off()

# Contrast with hydroxyBA
WasSurvDF_H <- WasSurvDF[WasSurvDF$chemical == "2-HydroxyBA" | WasSurvDF$chemical == "2-MethoxyBA" | WasSurvDF$chemical == "2-EthoxyBA" | WasSurvDF$chemical == "2-BromoBA",]
WasSurvDF_H$chemical <- as.factor(as.character(WasSurvDF_H$chemical))
WasSurvDF_H$chemical <- factor(WasSurvDF_H$chemical , levels = c("2-HydroxyBA","2-MethoxyBA","2-BromoBA", "2-EthoxyBA"))
WasFit_H <- survfit(Surv(time, status) ~ chemical, data = WasSurvDF_H)
WasFit_H.coxph <- coxph(Surv(time, status) ~ chemical, data = WasSurvDF_H)

# Contrast with bromoBA
WasSurvDF_B <- WasSurvDF[WasSurvDF$chemical == "2-BromoBA" | WasSurvDF$chemical == "2-MethoxyBA" | WasSurvDF$chemical == "2-EthoxyBA" ,]
WasSurvDF_B$chemical <- as.factor(as.character(WasSurvDF_B$chemical))
WasSurvDF_B$chemical <- factor(WasSurvDF_B$chemical , levels = c("2-BromoBA","2-MethoxyBA", "2-EthoxyBA"))
WasFit_B <- survfit(Surv(time, status) ~ chemical, data = WasSurvDF_B)
WasFit_B.coxph <- coxph(Surv(time, status) ~ chemical, data = WasSurvDF_B)

# Contrast methoxyBA ethoxyBA
WasSurvDF_ME <- WasSurvDF[WasSurvDF$chemical == "2-MethoxyBA" | WasSurvDF$chemical == "2-EthoxyBA",]
WasSurvDF_ME$chemical <- as.factor(as.character(WasSurvDF_ME$chemical))
WasSurvDF_ME$chemical <- factor(WasSurvDF_ME$chemical , levels = c("2-MethoxyBA", "2-EthoxyBA"))
WasFit_ME <- survfit(Surv(time, status) ~ chemical, data = WasSurvDF_ME)
WasFit_ME.coxph <- coxph(Surv(time, status) ~ chemical, data = WasSurvDF_ME)

# Paratrechina
ParSurvDF <- data.frame(time = as.vector(as.matrix(ParDF)), chemical = rep(names(ParDF), each = nrow(ParDF)))
ParSurvDF$status <- ParSurvDF$time
ParSurvDF$status[ParSurvDF$status<9000] <- 2
ParSurvDF$status[ParSurvDF$status == 9000] <- 1

ParSurvDF_cont <- ParSurvDF
ParSurvDF_cont$chemical <- as.factor(as.character(ParSurvDF_cont$chemical))
ParSurvDF_cont$chemical <- factor(ParSurvDF_cont$chemical , levels = c("Control","2-HydroxyBA","2-MethoxyBA","2-BromoBA", "2-EthoxyBA"))
ParFit_cont <- survfit(Surv(time, status) ~ chemical, data = ParSurvDF_cont)
ParFit_cont.coxph <- coxph(Surv(time, status) ~ chemical, data = ParSurvDF_cont)
jpeg('ParSurvStats.jpg', width=3000, height=1500, unit='px')
a <- ggforest(ParFit_cont.coxph, data = ParSurvDF_cont,
              noDigits = 2,
              cpositions=c(0.01, 0.08, 0.25),
              fontsize = 3,
              main = "")
a + theme(plot.margin = unit(c(2,2,2,2), "cm"), element_line(size = 10))
dev.off()
jpeg('ParSurvCurves.jpg', width=6000, height=3000, unit='px')
Par_SC <- ggsurvplot(ParFit_cont,
                     ggtheme = theme_classic2(base_size=100, base_family = "serif"),
                     palette = rev(ParCol),
                     size = 10,
                     surv.scale = "percent",
                     censor = FALSE) + 
  ylab("Percentage within perimeter") +
  xlab("Time (s)")
Par_SC
dev.off()

# Contrasting with hydroxyBA
ParSurvDF_H <- ParSurvDF[ParSurvDF$chemical == "2-HydroxyBA" | ParSurvDF$chemical == "2-MethoxyBA" | ParSurvDF$chemical == "2-EthoxyBA" | ParSurvDF$chemical == "2-BromoBA",]
ParSurvDF_H$chemical <- as.factor(as.character(ParSurvDF_H$chemical))
ParSurvDF_H$chemical <- factor(ParSurvDF_H$chemical , levels = c("2-HydroxyBA","2-MethoxyBA","2-BromoBA", "2-EthoxyBA"))
ParFit_H <- survfit(Surv(time, status) ~ chemical, data = ParSurvDF_H)
ParFit_H.coxph <- coxph(Surv(time, status) ~ chemical, data = ParSurvDF_H)

# Solenopsis
SolSurvDF <- data.frame(time = as.vector(as.matrix(SolDF)), chemical = rep(names(SolDF), each = nrow(SolDF)))
SolSurvDF$status <- SolSurvDF$time
SolSurvDF$status[SolSurvDF$status<1200] <- 2
SolSurvDF$status[SolSurvDF$status == 1200] <- 1

SolSurvDF_cont <- SolSurvDF
SolSurvDF_cont$chemical <- as.factor(as.character(SolSurvDF_cont$chemical))
SolSurvDF_cont$chemical <- factor(SolSurvDF_cont$chemical , levels = c("Control","2-HydroxyBA","2-MethoxyBA","2-BromoBA", "2-EthoxyBA"))
SolFit_cont <- survfit(Surv(time, status) ~ chemical, data = SolSurvDF_cont)
SolFit_cont.coxph <- coxph(Surv(time, status) ~ chemical, data = SolSurvDF_cont)
jpeg('SolSurvStats.jpg', width=3000, height=1500, unit='px')
a <- ggforest(SolFit_cont.coxph, data = SolSurvDF_cont,
              noDigits = 2,
              cpositions=c(0.01, 0.08, 0.25),
              fontsize = 3,
              main = "")
a + theme(plot.margin = unit(c(2,2,2,2), "cm"), element_line(size = 10))
dev.off()
jpeg('SolSurvCurves.jpg', width=6000, height=3000, unit='px')
Sol_SC <- ggsurvplot(SolFit_cont,
                     ggtheme = theme_classic2(base_size=100, base_family = "serif"),
                     palette = rev(SolCol),
                     size = 10,
                     surv.scale = "percent",
                     censor = FALSE) + 
  ylab("Percentage within perimeter") +
  xlab("Time (s)")
Sol_SC
dev.off()

# Contrast with bromoBA
SolSurvDF_B <- SolSurvDF[SolSurvDF$chemical == "2-BromoBA" | SolSurvDF$chemical == "2-MethoxyBA" | SolSurvDF$chemical == "2-EthoxyBA" ,]
SolSurvDF_B$chemical <- as.factor(as.character(SolSurvDF_B$chemical))
SolSurvDF_B$chemical <- factor(SolSurvDF_B$chemical , levels = c("2-BromoBA","2-MethoxyBA", "2-EthoxyBA"))
SolFit_B <- survfit(Surv(time, status) ~ chemical, data = SolSurvDF_B)
SolFit_B.coxph <- coxph(Surv(time, status) ~ chemical, data = SolSurvDF_B)

# Tapinoma
TapSurvDF <- data.frame(time = as.vector(as.matrix(TapDF)), chemical = rep(names(TapDF), each = nrow(TapDF)))
TapSurvDF$status <- TapSurvDF$time
TapSurvDF$status[TapSurvDF$status<1200] <- 2
TapSurvDF$status[TapSurvDF$status == 1200] <- 1

TapSurvDF_cont <- TapSurvDF
TapSurvDF_cont$chemical <- as.factor(as.character(TapSurvDF_cont$chemical))
TapSurvDF_cont$chemical <- factor(TapSurvDF_cont$chemical , levels = c("Control","2-HydroxyBA","2-MethoxyBA","2-BromoBA", "2-EthoxyBA"))
TapFit_cont <- survfit(Surv(time, status) ~ chemical, data = TapSurvDF_cont)
TapFit_cont.coxph <- coxph(Surv(time, status) ~ chemical, data = TapSurvDF_cont)
jpeg('TapSurvStats.jpg', width=3000, height=1500, unit='px')
a <- ggforest(TapFit_cont.coxph, data = TapSurvDF_cont,
              noDigits = 2,
              cpositions=c(0.01, 0.08, 0.25),
              fontsize = 3,
              main = "")
a + theme(plot.margin = unit(c(2,2,2,2), "cm"), element_line(size = 10))
dev.off()
jpeg('TapSurvCurves.jpg', width=6000, height=3000, unit='px')
Tap_SC <- ggsurvplot(TapFit_cont,
                     ggtheme = theme_classic2(base_size=100, base_family = "serif"),
                     palette = rev(TapCol),
                     size = 10,
                     surv.scale = "percent",
                     censor = FALSE) + 
  ylab("Percentage within perimeter") +
  xlab("Time (s)")
Tap_SC
dev.off()

# Contrasting with hydroxyBA
TapSurvDF_H <- TapSurvDF[TapSurvDF$chemical == "2-HydroxyBA" | TapSurvDF$chemical == "2-MethoxyBA" | TapSurvDF$chemical == "2-EthoxyBA" | TapSurvDF$chemical == "2-BromoBA",]
TapSurvDF_H$chemical <- as.factor(as.character(TapSurvDF_H$chemical))
TapSurvDF_H$chemical <- factor(TapSurvDF_H$chemical , levels = c("2-HydroxyBA","2-MethoxyBA","2-BromoBA", "2-EthoxyBA"))
TapFit_H <- survfit(Surv(time, status) ~ chemical, data = TapSurvDF_H)
TapFit_H.coxph <- coxph(Surv(time, status) ~ chemical, data = TapSurvDF_H)

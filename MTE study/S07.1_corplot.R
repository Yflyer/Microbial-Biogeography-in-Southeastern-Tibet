library(corrplot)
library(lme4)
library(rsq)
library(car) # mixed model anova
#################### save dir #######################
folder = paste0('S07')
dir.create(folder)
########
Env = read.csv('S01/MTE_ENV.csv',row.names = 1)
All.var = c('Temp.','pH','TN','SOC','Water_content','NH4')

cor_matrix <- cor(Env[,All.var])
cor_mtest <- cor.mtest(Env[,All.var], conf.level = .95)

pdf(file = paste0(folder,"/Env_corrplot.pdf"))
corrplot(cor_matrix, method = "color",
         type = "upper",
         order="hclust",
         addCoef.col = "white",
         tl.col="black",
         number.cex = 0.8,
         p.mat = cor_mtest$p,
         sig.level = 0.05,
         insig = "blank",
         #pch.cex = 0.8,
         #pch.col = "black",
         number.digits = 3)

dev.off()
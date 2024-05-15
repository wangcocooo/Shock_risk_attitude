df_AP <- readxl::read_xlsx('Country.xlsx', sheet = "AP")
df_DS <- readxl::read_xlsx('Country.xlsx', sheet = "DS")

library(dplyr)
library(data.table)
library(RColorBrewer)

df_AP <- data.matrix(df_AP)
df_AP <- df_AP[,-1]
rownames(df_AP) <- colnames(df_AP)

df_DS <- data.matrix(df_DS)
df_DS <- df_DS[,-1]
rownames(df_DS) <- colnames(df_DS)

heatmap.2(df_AP, symm = T, trace = 'none', key = F, Rowv = NA)
heatmap.2(df_AP, symm = TRUE, Rowv = NA, scale = 'none',dendrogram='none',
              xlab = "Countries", ylab =  "Countries",
              main = "AP aversion",
              col = colorRampPalette(brewer.pal(8,"Blues"))(4),
          key = F, trace = 'none')

heatmap.2(df_DS, symm = TRUE, Rowv = NA, scale = 'none',
        xlab = "Countries", ylab =  "Countries",
        main = "DS aversion",
        col = colorRampPalette(brewer.pal(8,"Blues"))(4),
        key = F, trace = 'none')

melt_AP <- melt(df_AP, na.rm = T)

p_AP <- ggplot(data = melt_AP, aes(x = Var1, y = Var2, fill = factor(value)))+
  geom_tile(lty = 3, col = 'black')+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+
  xlab('Countries') + ylab('Countries')+
  scale_fill_manual(name = 'Prob',
                    values = c("grey100", "grey80", "grey40",'grey3'),
                    labels = c('Prob < 1%','1% <= Prob < 5%','5% <= Prob < 10%','10% <= Prob'))+
  coord_fixed()

melt_DS <- melt(df_DS, na.rm = T)
p_DS <- ggplot(data = melt_DS, aes(x = Var1, y = Var2, fill = factor(value)))+
  geom_tile(lty = 3, col = 'black')+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+
  xlab('Countries') + ylab('Countries')+
  scale_fill_manual(name = 'Prob',
                    values = c("grey100", "grey80", "grey40",'grey3'),
                    labels = c('Prob < 1%','1% <= Prob < 5%','5% <= Prob < 10%','10% <= Prob'))+
  coord_fixed()


jpeg('heatmap.jpeg', width = 10, height = 6, units = 'in', res = 200)
ggpubr::ggarrange(p_AP,p_DS, ncol = 2, nrow = 1, common.legend = T, legend = 'bottom',
                  font.label = list(size = 12), labels = c('(a) AP aversion','(b) DS aversion'),
                  vjust = 2)
dev.off()  
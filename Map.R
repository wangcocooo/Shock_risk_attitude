df_coun <- readxl::read_xlsx('Country attitudes.xlsx')

library(maps)
library(ggplot2)
library(dplyr)

world <- map_data('world')
eu <- subset(world, region %in% c('Belgium','Bulgaria','Cyprus','Czech Republic','Denmark','Germany','Greece',
                                  'Spain','Estonia','France','Croatia','Hungary','Ireland','Italy',
                                  'Lithuania','Luxembourg','Latvia','Malta','Netherlands','Austria',
                                  'Poland','Portugal','Romania','Finland','Sweden','Slovakia','Slovenia'))

eu_plot <- merge(eu, df_coun, by.x = "region", by.y = "Country", all.x = T)%>%
  arrange(group, order)


AP <- ggplot(eu_plot, aes(x = long, y = lat, group = group, fill = AP*2)) +
  scale_fill_gradient(low = 'grey90',high = 'black', name = 'AP')+
  geom_polygon(color = 'white', size =0.4, lty = 1)+
  theme_minimal()+
  xlab('Longitude')+
  ylab('Latitude')

DS <- ggplot(eu_plot, aes(x = long, y = lat, group = group, fill = DS*-6)) +
  scale_fill_gradient(low = 'grey90',high = 'black', name = 'DS')+
  geom_polygon(color = 'white', size =0.4, lty = 1)+
  theme_minimal()+
  xlab('Longitude')+
  ylab('Latitude')

jpeg('country attitudes.jpeg', width = 10, height = 3.5, units = 'in', res = 200)
ggpubr::ggarrange(AP,DS, ncol = 2, nrow = 1, font.label = list(size = 12), legend = 'right',labels = c('(a) AP aversion','(b) DS aversion'),
                  vjust = -.1, hjust = 0)+
  theme(plot.margin = margin(0.7,0,0,0, "cm"))
dev.off() 

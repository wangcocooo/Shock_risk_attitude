coef <- readxl::read_xlsx('/Users/cocowang/Documents/Research Practice/Results/sure_entire_untrimmed.xlsx', sheet = 'R_data')

library(ggplot2)
library(geomtextpath)

df <- with(coef,
                data.frame(value = c(AP, DS),
                           variable = factor(rep(c("AP","DS"),
                                                 each = NROW(coef))),
                           Year = Year))

p1 <- ggplot2::ggplot(df, aes(Year, value, lty = variable)) + geom_line()+
  geom_textvline(label = '2008 shock',xintercept = 2008, col = 'grey')+
  geom_textvline(label = '2011 shock',xintercept = 2011, col = 'grey')+
  theme_bw() + theme(legend.position="bottom")+
  xlab('Year')+
  ylab('Risk aversion coefficients')+labs(lty="Risk aversion type")
p1

jpeg('risk aversion.jpeg', width = 7.5, height = 4, units = 'in', res = 200)
ggpubr::ggarrange(p1, ncol = 1, nrow = 1, font.label = list(size = 10))
dev.off()  

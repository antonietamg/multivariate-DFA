plotdfa <- ggplot(dplot, aes(x=sx, y=fx)) + 
  #stat_smooth(method='lm', formula = y~poly(x,1), col="steelblue1")+
  geom_point(col="firebrick2")+
  geom_point(shape = 1, colour = "black")+
  labs(title="your_title", x="log(s)", y="log(F(s))") #change the title in labs(title="")

plotdfa1 <- plotdfa + theme(plot.title = element_text(hjust = 0.5, size=15, face="bold"))

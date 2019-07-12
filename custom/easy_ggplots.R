## Easier ggplots

# ggplot axis angle adjust
gg_axis_adjust_angle <- function(x=0, y=0) {
	theme(
		axis.text.x = element_text(angle = x, hjust = 1),
		axis.text.y = element_text(angle = y, vjust = 1))
}


## Mimic base R plots

gg_base_legend <- function() {
	theme(legend.position = c(1,1), legend.justification = c(1,1),
		legend.background = element_rect(fill=NA, linetype='solid', color='grey45'),
		legend.key = element_rect(fill=NA, linetype='blank'))
}


set.seed(42)
plotDF <- data.frame(
	value = c(rnorm(100,2), rnorm(100,-1)), 
	tag = rep(c('a','b'), each = 100))

## base plot reference
d1 = density(plotDF$value)
d2 = density(plotDF$value[which(plotDF$tag == 'a')])
d3 = density(plotDF$value[which(plotDF$tag == 'b')])
ylim = range(c(d1$y, d2$y, d3$y))
plot(d1, main = 'density plots', col = 'grey25', ylim = ylim)
lines(d2, col = 'green')
lines(d3, col = 'blue')
legend('topright', legend=c('all','a','b'), 
	col=c('grey25','blue','green'), lty=1, xjust=1)

## ggplot mimic
ggplot(plotDF, aes(x=value)) +
  geom_line(stat='density', col='grey25') +
  geom_line(aes(colour=tag), stat='density') +
  gg_base_legend() + 
  scale_x_continuous(breaks=seq(-4,4,1))


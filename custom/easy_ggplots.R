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


# set.seed(42)
# plotDF <- data.frame(
# 	value = c(rnorm(100,2), rnorm(100,-1)), 
# 	tag = rep(c('a','b'), each = 100))

# ## base plot reference
# d1 = density(plotDF$value)
# d2 = density(plotDF$value[which(plotDF$tag == 'a')])
# d3 = density(plotDF$value[which(plotDF$tag == 'b')])
# ylim = range(c(d1$y, d2$y, d3$y))
# plot(d1, main = 'density plots', col = 'grey25', ylim = ylim)
# lines(d2, col = 'green')
# lines(d3, col = 'blue')
# legend('topright', legend=c('all','a','b'), 
# 	col=c('grey25','blue','green'), lty=1, xjust=1)

## ggplot mimic
ggplot(plotDF, aes(x=value)) +
  geom_line(stat='density', col='grey25') +
  geom_line(aes(colour=tag), stat='density') +
  # gg_base_legend() + 
  scale_x_continuous(breaks=seq(-4,4,1))


##  Themes!  ##
### source : https://www.shanelynn.ie/themes-and-colours-for-r-ggplots-with-ggthemr/

# reference
theme_grey2 <- function (base_size = 11, base_family = "") 
{
	half_line <- base_size/2
	theme(
		line = element_line(colour = "black", size = 0.5, linetype = 1, lineend = "butt"), 
		rect = element_rect(fill = "white", colour = "black", size = 0.5, linetype = 1), 
		text = element_text(family = base_family, face = "plain", colour = "black", size = base_size, 
			lineheight = 0.9, hjust = 0.5, vjust = 0.5, angle = 0, margin = margin(), debug = FALSE), 

		# Axis
		axis.line = element_line(), 
		axis.line.x = element_blank(), 
		axis.line.y = element_blank(), 
		axis.text = element_text(size = rel(0.8), colour = "grey30"), 
		axis.text.x = element_text(margin = margin(t = 0.8 * half_line/2), vjust = 1), 
		axis.text.y = element_text(margin = margin(r = 0.8 * half_line/2), hjust = 1), 
		axis.ticks = element_line(colour = "grey20"), 
		axis.ticks.length = unit(half_line/2, "pt"), 
		axis.title.x = element_text(margin = margin(t = 0.8 * half_line, b = 0.8 * half_line/2)), 
		axis.title.y = element_text(angle = 90, margin = margin(r = 0.8 * half_line, l = 0.8 * half_line/2)), 

		# Legends
		legend.background = element_rect(colour = NA), 
		legend.margin = unit(0.2, "cm"), 
		legend.key = element_rect(fill = "grey95", colour = "white"), 
		legend.key.size = unit(1.2, "lines"), 
		legend.key.height = NULL, 
		legend.key.width = NULL, 
		legend.text = element_text(size = rel(0.8)), 
		legend.text.align = NULL, 
		legend.title = element_text(hjust = 0), 
		legend.title.align = NULL, 
		legend.position = "right", 
		legend.direction = NULL, 
		legend.justification = "center", 
		legend.box = NULL, 

		# Panel
		panel.background = element_rect(fill = "grey92", colour = NA), 
		panel.border = element_blank(), 
		panel.grid.major = element_line(colour = "white"), 
		panel.grid.minor = element_line(colour = "white", size = 0.25), 
		panel.margin = unit(half_line, "pt"), 
		panel.margin.x = NULL, 
		panel.margin.y = NULL, 
		panel.ontop = FALSE, 

		strip.background = element_rect(fill = "grey85", colour = NA), 
		strip.text = element_text(colour = "grey10", size = rel(0.8)), 
		strip.text.x = element_text(margin = margin(t = half_line, b = half_line)), 
		strip.text.y = element_text(angle = -90, margin = margin(l = half_line, r = half_line)), 
		strip.switch.pad.grid = unit(0.1, "cm"), 
		strip.switch.pad.wrap = unit(0.1, "cm"), 

		plot.background = element_rect(colour = "white"), 
		plot.title = element_text(size = rel(1.2), margin = margin(b = half_line * 1.2)), 
		plot.margin = margin(half_line, half_line, half_line, half_line), 
		complete = TRUE)
}

# set theme
theme_set(theme_grey2())
theme_get()



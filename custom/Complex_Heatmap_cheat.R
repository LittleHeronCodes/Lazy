## ComplexHeatmap grob components

mtitle = 'This is a Title'
hp = Heatmap(plotMat, col=colors, name = 'mut', column_title = mtitle, row_title = 'Row', 
	cluster_rows = FALSE, show_row_dend = FALSE, show_row_names = TRUE,
	cluster_columns = FALSE, show_column_dend = FALSE, show_column_names = TRUE,
	top_annotation = cell_ha)

# Title
grobHeight(textGrob(mtitle, gp = gpar(fontsize = 14))) + unit(5, "mm") + unit(1,'mm')

# dendrobram
unit(10, "mm") + unit(.5, "mm")

# annotation size
ht_opt$simple_anno_size
cell_ha@height

# column names
max_text_width(colnames(plotMat), gp = gpar(fontsize = 12))

max_text_width(row(plotMat), gp = gpar(fontsize = 12))



hopt <- list(
	cluster_rows = TRUE, show_row_dend = FALSE, show_row_names = TRUE,
	cluster_columns = TRUE, show_column_dend = FALSE, show_column_names = TRUE
	)


add_text_cheatmap <- function(text, col='grey15') {
	cell_fun = function(j,i,x,y,w,h,col) {
		grid.text(text[i,j], x, y, gp=gpar(col=col, fontface='bold'))		
	}
	return(cell_fun)
}

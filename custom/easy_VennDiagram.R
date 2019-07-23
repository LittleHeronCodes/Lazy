
vennSimple = function(vennLs, mtitle = '', fname = 'venn.png', 
	colors = NULL, scale = FALSE, h=800, w=800) {
	require(RColorBrewer)
	require(VennDiagram)

	if(is.null(colors)) colors = RColorBrewer::brewer.pal(8, "Set2")[1:length(vennLs)]

	vennLs = lapply(vennLs, unique)
	vennLs = lapply(vennLs, function(v) v[which(!is.na(v))] )

	venn.diagram(vennLs, filename = fname, na = 'remove', height = h, width = w, 
		main = mtitle, euler.d = scale, scaled = scale, alpha = 0.5, col = "transparent", 
		fill = colors, cex = 0.3, cat.cex = 0.3, main.cex = 0.4, margin = 0.1, 
		imagetype = "png")
}


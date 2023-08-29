custom_enrichR_plot = function(df, showTerms = 20, numChar = 40, y = "gene.count", orderBy = "Adjusted.p.value", 
          xlab = NULL, ylab = NULL, title = NULL) 
{
  if (!is.data.frame(df)) {
    stop("Input df is malformed - must be a data.frame object.")
  }
  if (nrow(df) == 0 | ncol(df) == 0) {
    stop("Input df is empty.")
  }
  if (!is.numeric(numChar)) {
    stop(paste0("numChar '", numChar, "' is invalid."))
  }
  df = df[c(1:showTerms),] 
  shortName <- paste(substr(df$Term, 1, numChar), ifelse(nchar(df$Term) > 
                                                           numChar, "...", ""), sep = "")
  names(shortName) <- df$Term.name
  if (any(duplicated(shortName))) {
    warning("There are duplicated trimmed names in the plot, consider increasing the 'numChar' setting.")
  }
  #df$Ratio <- df$Significant/df$Annotated
  if (orderBy == "Combined.Score") {
    fill <- "Combined.Score"
  }
  else {
    fill <- "Adjusted.p.value"
  }

  map <- aes_string(x = "Term.name", y = y, fill = fill)
  if (is.null(xlab)) {
    xlab <- "Enriched terms"
  }
  if (is.null(ylab)) {
    ylab <- "Gene count"
  }
  if (is.null(title)) {
    title <- "Enrichment Analysis by EnrichR"
  }
  p <- ggplot(df, map) + geom_bar(stat = "identity") + coord_flip() + 
    theme_bw() + scale_x_discrete(labels = rev(shortName), 
                                  limits = rev(df$Term))
  if (orderBy == "Combined.Score") {
    p <- p + scale_fill_continuous(low = "blue", high = "red") + 
      guides(fill = guide_colorbar(title = "Combined Score", 
                                   reverse = FALSE))
  }
  else {
    p <- p + scale_fill_continuous(low = "red", high = "blue") + 
      guides(fill = guide_colorbar(title = "P value", reverse = TRUE))
  }
  p <- p + theme(axis.text.x = element_text(colour = "black", 
                                            vjust = 1), axis.text.y = element_text(colour = "black", 
                                                                                   hjust = 1), axis.title = element_text(color = "black", 
                                                                                                                         margin = margin(10, 5, 0, 0)), axis.title.y = element_text(angle = 90))
  p <- p + xlab(xlab) + ylab(ylab) + ggtitle(title)
  return(p)
}


plot.weights <- function(m, exlabs, group = "Variable") {
  
  ### Get weights from model fit object 
  weights <- rbind(data.frame(value = m$pos.weights, name = names(m$pos.weights)),
                   data.frame(value = m$neg.weights * -1, name = names(m$neg.weights)))
  weights$name <- weights$name |> plyr::mapvalues(from = exlabs$var, to = exlabs$label, warn_missing = FALSE)
  
  ### Set colors based on psi values
  if (m$pos.psi > abs(m$neg.psi)) {
    color <- ifelse(weights$value <= 0, "grey75", "grey35")
  } else {
    color <- ifelse(weights$value > 0, "grey75", "grey35")
  }
  
  ### Plot
  g <- ggplot(weights, aes(x = reorder(name, abs(value)), y = value)) + 
    theme_bw() + 
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          text = element_text(size=15),
          axis.title.y =element_text(size=16),
          axis.title.x =element_text(size=12)) + 
    #axis.text=element_text(size=rel(1.2))) + 
    geom_bar(stat = "identity",
             show.legend = FALSE,
             fill = color) + 
    scale_y_continuous(breaks=seq(-1,1,0.25), limits = c(-1.0, 1.0)) +
    xlab(group) + 
    ylab("Negative Weights     Positive Weights") + 
    coord_flip()
  
  return(g)
}

                 


  
  

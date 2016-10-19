#
# probability density functions
#

pdf.sif<-function(sif){
  return(density(x = sif$V2))
}

sif.density <- function(sif){
  return(density(x = sif$V2))
}

#
#Plotting functions
#

plot_density<-function(...){
  dots = list(...)
  colores = rainbow(n = length(seq_along(dots)))

  plot(x = dots[[1]]$x,
       y = dots[[1]]$y,
       xlab = "MI",
       ylab = "Frequency",
       col = colores[1],
       type = "l")

  for(i in seq_along(dots)){
    if(i > 1){
      points(x = dots[[i]]$x,
             y = dots[[i]]$y,
             type = "l",
             col = colores[i])
    }
  }
}

plot_density_log<-function(...){
  dots = list(...)
  colores = rainbow(n = length(seq_along(dots)))

  plot(x = dots[[1]]$x,
       y = dots[[1]]$y,
       xlab = "MI",
       ylab = "log(Frequency)",
       col = colores[1],
       type = "l",
       log = "y")

  for(i in seq_along(dots)){
    if(i > 1){
      points(x = dots[[i]]$x,
             y = dots[[i]]$y,
             type = "l",
             col = colores[i])
    }
  }
}

plot_density_zoom<-function(..., zoom = 10){
  dots = list(...)
  colores = rainbow(n = length(seq_along(dots)))
  zm = round(length(dots[[1]]$y)/zoom)

  plot(x = dots[[1]]$x[1:zm],
       y = dots[[1]]$y[1:zm],
       xlab = "MI",
       ylab = "Frequency",
       col = colores[1],
       type = "l",
       #log = "y"
  )

  for(i in seq_along(dots)){
    if(i > 1){
      points(x = dots[[i]]$x[1:zm],
             y = dots[[i]]$y[1:zm],
             type = "l",
             col = colores[i])
    }
  }
}

#
# Plot (write out version)
#
plot_density.wo<-function(..., file = "plot.pdf"){
  dots = list(...)
  colores = rainbow(n = length(seq_along(dots)))
  pdf(file = file, paper = "a4r")
  plot(x = dots[[1]]$x,
       y = dots[[1]]$y,
       xlab = "MI",
       ylab = "Frequency",
       col = colores[1],
       type = "l")

  for(i in seq_along(dots)){
    if(i > 1){
      points(x = dots[[i]]$x,
             y = dots[[i]]$y,
             type = "l",
             col = colores[i])
    }
  }
  dev.off()
}

plot_density_log.wo<-function(..., file = "plot.pdf"){
  dots = list(...)
  colores = rainbow(n = length(seq_along(dots)))

  pdf(file = file, paper = "a4r")
  plot(x = dots[[1]]$x,
       y = dots[[1]]$y,
       xlab = "MI",
       ylab = "log(Frequency)",
       col = colores[1],
       type = "l",
       log = "y")

  for(i in seq_along(dots)){
    if(i > 1){
      points(x = dots[[i]]$x,
             y = dots[[i]]$y,
             type = "l",
             col = colores[i])
    }
  }
  dev.off()
}

plot_density_zoom.wo<-function(..., zoom = 10, file = "plot.pdf"){
  dots = list(...)
  colores = rainbow(n = length(seq_along(dots)))
  zm = round(length(dots[[1]]$y)/zoom)
  pdf(file = file, paper = "a4r")
  plot(x = dots[[1]]$x[1:zm],
       y = dots[[1]]$y[1:zm],
       xlab = "MI",
       ylab = "Frequency",
       col = colores[1],
       type = "l",
       #log = "y"
  )

  for(i in seq_along(dots)){
    if(i > 1){
      points(x = dots[[i]]$x[1:zm],
             y = dots[[i]]$y[1:zm],
             type = "l",
             col = colores[i])
    }
  }
  dev.off()
}

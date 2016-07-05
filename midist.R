#LIBRARIES
library("data.table")

#
# Load index
#
load.index<-function(index){
  fread(input = index, data.table = FALSE)
}

#
#Subgraphing functions
#
intra.chromosomic<-function(sif, index){
  index<-as.data.frame(index, stringsAsFactors = FALSE)
  return(subset(sif, index$V1 ==index$V2 &
                  !is.na(index$V1) &
                  !is.na(index$V2)))
}

inter.chromosomic <-function(sif, index){
  index<-as.data.frame(index, stringsAsFactors = FALSE)
  return(subset(sif, index$V1 !=index$V2 &
                  !is.na(index$V1) &
                  !is.na(index$V2)))
}

selector<-function(sif, criteria){
  #example criteria index$V1 == "19"
  #example criteria index$V1 == "19" & index$V2 == "20"
  return(subset(sif, criteria))
}

#
# probability density functions
#

pdf.sif<-function(sif){
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


#
# Comparing via Mann Whitney
#

U_test<- function(pdf1, pdf2){

  test<-wilcox.test(x = pdf1$y,
                    y = pdf2$y,
                    alternative = "two.sided",
                    paired = FALSE,
                    exact = FALSE)
  return(test)
}

# PENDING
# Shuffle interactions test
#

shuffle_interactions_test<-function(sif, index){
  #shuffle interactions
  dex_shuffle <- sample(as.data.frame(index))
  #shuffled networks
  intranw <- intra(sif = sif, index = dex_shuffle)
  internw <- inter(sif = sif, index = dex_shuffle)
  #density
  d.intra <- density(intranw$V2)
  d.inter <- density(internw$V2)
  #KW
  return(U_test(d.intra, d.inter))

}

shuffle_repeat <- function(sif, index, n){
  L<-list()
  for(i in 1:n){
    name = paste(i)
    k <- shuffle_interactions_test(sif, index)
    L <- c(L, i = as.list(k))
  }
  return(L)
}

require("data.table")
require("parallel")
enfermos_sif <-fread(input = "/mnt/e/jespinal/sif/Enfermos1_SinDuplicados.sif", data.table = FALSE)
sanos_sif <-fread(input = "/mnt/e/jespinal/sif/Sanos1_SinDuplicados.sif", data.table = FALSE)

#histogramas
hist_enfermos <- hist(enfermos_sif$V2, breaks = 200, plot = FALSE)
hist_sanos <- hist(sanos_sif$V2, breaks = 200, plot = FALSE)

##plots histogramas
plot(x = hist_sanos$mids, y = log(hist_sanos$count), type = "l")
lines(x = hist_enfermos$mids, y = log(hist_enfermos$count), type = "l", col = "red")

#TODO: Check whether this is the right way to estimate the distribution.
#density plots
density_enfermos<-density(enfermos_sif$V2)
density_sanos<-density(sanos_sif$V2)

##plot density
plot(x = density_enfermos, type = "l", col = "red")
lines(x = density_sanos, type = "l", col = "blue")

plot(x = density_sanos, type = "l", col = "green")
#compare densities?
ks.test(x = density_enfermos$y, y = density_sanos$y)

# #try with KernSmooth
# bk_enfermos<-bkde(enfermos_sif$V2)
# bk_sanos<-bkde(sanos_sif$V2)
# plot(bk_enfermos)
# plot(bk_sanos)
# ks.test(x = bk_enfermos$y, y = bk_sanos$y)

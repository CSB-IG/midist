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

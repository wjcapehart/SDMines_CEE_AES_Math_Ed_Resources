

library(package = "tidyverse")
library(package  = "lubridate")
library(package = "readxl")


Landsat_7_Samples = read_excel(path  = "./Landsat_7_Rapid_City_ETM_2001-08-16_1725_Classes.xlsx") 

Landsat_7_Samples$Class=as.factor( Landsat_7_Samples$Class )
Landsat_7_Samples$Color=as.factor( Landsat_7_Samples$Color )

colormap = unique(data_frame(Class = Landsat_7_Samples$Class,
                             Color = Landsat_7_Samples$Color))

Landsat_7_Samples_ColorMap = array(data     = colormap$Color,
                                   dimnames = list(Class = colormap$Class))

remove(colormap)

Landsat_7_Samples = subset(Landsat_7_Samples, select=-Color)   


# top of atmos radiance.


 Landsat_7_Samples$ETM1 = 0.766 * Landsat_7_Samples$ETM1 - 2.28583
 Landsat_7_Samples$ETM2 = 1.448 * Landsat_7_Samples$ETM2 - 4.28819
 Landsat_7_Samples$ETM3 = 1.044 * Landsat_7_Samples$ETM3 - 2.21398
 Landsat_7_Samples$ETM4 = 0.876 * Landsat_7_Samples$ETM4 - 2.38602
 Landsat_7_Samples$ETM5 = 0.120 * Landsat_7_Samples$ETM5 - 0.49035
 Landsat_7_Samples$ETM7 = 0.066 * Landsat_7_Samples$ETM7 - 0.21555
 
 Landsat_7_Bands_Metadata = array(data      =      c("0.45-0.52",
                                                     "0.52-0.60",
                                                     "0.63-0.69",
                                                     "0.77-0.90",
                                                     "1.55-1.75",
                                                     "2.09-2.35"),
                         dimnames = list(channel = c("ETM1",
                                                     "ETM2",
                                                     "ETM3",
                                                     "ETM4",
                                                     "ETM5",
                                                     "ETM7")))

 
    SUN_AZIMUTH = 136.0658436 
    SUN_ELEVATION = 56.5251637  
    SUN_ZENITH  = (90 - SUN_ELEVATION ) * pi / 360.0
    
    
    DATE = as.POSIXct(x = "2001-08-16 17:24:56 UTC", tz = "UTC") 
    DD   = (decimal_date(DATE) - year(DATE)) * 2 * pi / 365.25
    
    dmean_over_d_squared = 1.000110 + 0.034221*cos(DD)   + 0.001280*sin(DD) +
                                      0.000719*cos(2*DD) + 0.000077*sin(2*DD)
    
    
    Landsat_7_Time_Metadata = data.frame(date = DATE,
                                                      sun_elevation = SUN_ELEVATION,
                                                      sun_azimuth   = SUN_AZIMUTH,
                                                      sun_zenith    = SUN_ZENITH,
                                                      dmean_over_d_squared = dmean_over_d_squared)
      
    
    
 Landsat_7_Samples$ETM1 = pi * Landsat_7_Samples$ETM1 / 1997.0 / dmean_over_d_squared / cos(SUN_ZENITH)
 Landsat_7_Samples$ETM2 = pi * Landsat_7_Samples$ETM2 / 1812.0 / dmean_over_d_squared / cos(SUN_ZENITH)
 Landsat_7_Samples$ETM3 = pi * Landsat_7_Samples$ETM3 / 1533.0 / dmean_over_d_squared / cos(SUN_ZENITH)
 Landsat_7_Samples$ETM4 = pi * Landsat_7_Samples$ETM4/ 1039.0 / dmean_over_d_squared / cos(SUN_ZENITH)
 Landsat_7_Samples$ETM5 = pi * Landsat_7_Samples$ETM5 /  230.8 / dmean_over_d_squared / cos(SUN_ZENITH)
 Landsat_7_Samples$ETM7 = pi * Landsat_7_Samples$ETM7 /   84.9 / dmean_over_d_squared / cos(SUN_ZENITH)

 remove(DATE)
 remove(DD)
 remove(SUN_ELEVATION)
 remove(SUN_AZIMUTH)
 remove(SUN_ZENITH)
 remove(dmean_over_d_squared)
 
 
 
 
 
 







save(Landsat_7_Samples,
     Landsat_7_Bands_Metadata,
     Landsat_7_Time_Metadata,
     Landsat_7_Samples_ColorMap,
     file = "../Landsat_7_Rapid_City_ETM_2001-08-16_1725_Classes.Rdata") 




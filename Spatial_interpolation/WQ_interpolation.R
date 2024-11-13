library(tidyr)
library(dplyr)
library(purrr)
library(raster)
library(automap)

###-------- import Enterococcus -----#########
setwd("~/Spatial_interpolation")
Entero0<-read.csv("WQD_Enterococcus.csv")
Entero1<- Entero0 %>%
  separate(ActivityStartDate, c("Month", "Day", "Year"), "/")
Entero<-data.frame(cbind(Entero1$MonitoringLocationName, Entero1$ActivityLocation.LatitudeMeasure, Entero1$ActivityLocation.LongitudeMeasure, Entero1$Month, Entero1$Year, Entero1$ResultMeasureValue))
names(Entero)<-c("Site", "Lat", "Lon", "Month", "Year", "Entero")
Entero<-Entero %>% drop_na()
Entero<-Entero[Entero$Entero != 0, ]
# Identify and remove sites lacking sufficient data points
total.samples <- group_by(Entero, Site) %>%
  dplyr::summarise(n = n())
ggplot(total.samples)+
  aes(x = n) +
  geom_histogram(bins = 30, fill = "#0c4c8a") +
  theme_minimal()
Entero.data.filt1 <- filter(Entero, Site %in% total.samples[total.samples$n > 40,]$Site)
Entero.data.filt <- filter(Entero.data.filt1, Year > 10 & Year < 22)

Entero.data.filt$Entero<-as.numeric(Entero.data.filt$Entero)
Entero.data.filt <- Entero.data.filt[!is.na(as.numeric(Entero.data.filt$Entero)), ]   
n_distinct(Entero.data.filt$Site)
# Summarize data for each site
Entero.summ <- group_by(Entero.data.filt, Site) %>%
  dplyr::summarise(across(1:3, first),
                   across(5, list(mean = function(x) mean(x, na.rm = T), 
                                  max = function(x) max(x, na.rm = T), 
                                  min = function(x) min(x, na.rm = T))))
Entero.month <- group_by(Entero.data.filt, Site, Month) %>%
  dplyr::summarise(monthly_range = (max(Entero, na.rm=TRUE) - min(Entero, na.rm=TRUE)))
Entero.month2 <- group_by(Entero.month, Site) %>%
  dplyr::summarise(Entero_mean_monthly_range = mean(monthly_range, na.rm = T))

Entero.year <- group_by(Entero.data.filt, Site, Year) %>%
  dplyr::summarise(yearly_range = (max(Entero, na.rm=TRUE) - min(Entero, na.rm=TRUE)))
Entero.year2 <- group_by(Entero.year, Site) %>%
  dplyr::summarise(Entero_mean_yearly_range = mean(yearly_range, na.rm = T))

Entero.site <-list(Entero.summ, Entero.month2, Entero.year2) %>% reduce(inner_join, by = "Site")
Entero.site$Lat<-as.numeric(Entero.site$Lat)
Entero.site$Lon<-as.numeric(Entero.site$Lon)
Entero.site$Month <- NULL
write.csv(Entero.site, "Entero_STX.csv")
Entero.site<-read.csv("Entero_STX.csv")[,-1]

Entero.data.filt.year=Entero.data.filt %>% 
  select(Site, Year) %>% 
  unique() %>% 
  group_by(Year) %>% 
  summarize(Enteron = n())
Entero.data.filt.month=Entero.data.filt %>% 
  select(Site, Month) %>% 
  unique() %>% 
  group_by(Month) %>% 
  summarize(Enteron = n())


###-------- import Dissolved Oxygen -----#########
setwd("~/Spatial_interpolation")
DO0<-read.csv("WQD_DO.csv")
DO1<- DO0 %>%
  separate(ActivityStartDate, c("Month", "Day", "Year"), "/")
DO<-data.frame(cbind(DO1$MonitoringLocationName, DO1$ActivityLocation.LatitudeMeasure, DO1$ActivityLocation.LongitudeMeasure, DO1$Month, DO1$Year, DO1$ResultMeasureValue))
names(DO)<-c("Site", "Lat", "Lon", "Month", "Year", "DO")
DO<-DO %>% drop_na()
DO <- DO[DO$DO != 0, ]
# Identify and remove sites lacking sufficient data points
total.samples <- group_by(DO, Site) %>%
  dplyr::summarise(n = n())
ggplot(total.samples)+
  aes(x = n) +
  geom_histogram(bins = 30, fill = "#0c4c8a") +
  theme_minimal()
DO.data.filt1 <- filter(DO, Site %in% total.samples[total.samples$n > 50,]$Site)
DO.data.filt <- filter(DO.data.filt1, Year > 10 & Year < 23)

DO.data.filt$DO<-as.numeric(DO.data.filt$DO)
DO.data.filt <- DO.data.filt[!is.na(as.numeric(DO.data.filt$DO)), ]   
n_distinct(DO.data.filt$Site)
# Summarize data for each site
DO.summ <- group_by(DO.data.filt, Site) %>%
  dplyr::summarise(across(1:3, first),
                   across(5, list(mean = function(x) mean(x, na.rm = T), 
                                  max = function(x) max(x, na.rm = T), 
                                  min = function(x) min(x, na.rm = T))))
DO.month <- group_by(DO.data.filt, Site, Month) %>%
  dplyr::summarise(monthly_range = (max(DO, na.rm=TRUE) - min(DO, na.rm=TRUE)))
DO.month2 <- group_by(DO.month, Site) %>%
  dplyr::summarise(DO_mean_monthly_range = mean(monthly_range, na.rm = T))

DO.year <- group_by(DO.data.filt, Site, Year) %>%
  dplyr::summarise(yearly_range = (max(DO, na.rm=TRUE) - min(DO, na.rm=TRUE)))
DO.year2 <- group_by(DO.year, Site) %>%
  dplyr::summarise(DO_mean_yearly_range = mean(yearly_range, na.rm = T))

DO.site <-list(DO.summ, DO.month2, DO.year2) %>% reduce(inner_join, by = "Site")
DO.site$Lat<-as.numeric(DO.site$Lat)
DO.site$Lon<-as.numeric(DO.site$Lon)
DO.site$Month <- NULL
write.csv(DO.site, "DO_STX.csv")
DO.site<-read.csv("DO_STX.csv")[,-1]

DO.data.filt.year=DO.data.filt %>% 
  select(Site, Year) %>% 
  unique() %>% 
  group_by(Year) %>% 
  summarize(DOn = n())
DO.data.filt.month=DO.data.filt %>% 
  select(Site, Month) %>% 
  unique() %>% 
  group_by(Month) %>% 
  summarize(DOn = n())


###-------- import E. coli -----#########
setwd("~/Spatial_interpolation")
Ecoli0<-read.csv("WQD_Fecal_Coliform.csv")
Ecoli1<- Ecoli0 %>%
  separate(ActivityStartDate, c("Month", "Day", "Year"), "/")
Ecoli2<-data.frame(cbind(Ecoli1$MonitoringLocationName, Ecoli1$ActivityLocation.LatitudeMeasure, Ecoli1$ActivityLocation.LongitudeMeasure, Ecoli1$Month, Ecoli1$Year, Ecoli1$ResultMeasureValue))
names(Ecoli2)<-c("Site", "Lat", "Lon", "Month", "Year", "Ecoli")
Ecoli2 <- Ecoli2[Ecoli2$Ecoli != 0, ]
Ecoli2 <- with(Ecoli2, Ecoli2[!(Ecoli == "" | is.na(Ecoli)), ])
Ecoli2 <- with(Ecoli2, Ecoli2[!(Ecoli == "<1" | is.na(Ecoli)), ])
Ecoli2 <- with(Ecoli2, Ecoli2[!(Ecoli == "<10" | is.na(Ecoli)), ])
Ecoli<-Ecoli2 %>% drop_na()
total.samples <- group_by(Ecoli, Site) %>%
  dplyr::summarise(n = n())
ggplot(total.samples)+
  aes(x = n) +
  geom_histogram(bins = 30, fill = "#0c4c8a") +
  theme_minimal()
Ecoli.data.filt1 <- filter(Ecoli, Site %in% total.samples[total.samples$n > 5,]$Site)
Ecoli.data.filt <- filter(Ecoli.data.filt1, Year > 00 & Year < 14)

Ecoli.data.filt$Ecoli<-as.numeric(Ecoli.data.filt$Ecoli)
Ecoli.data.filt <- Ecoli.data.filt[!is.na(as.numeric(Ecoli.data.filt$Ecoli)), ]   
n_distinct(Ecoli.data.filt$Site)
# Summarize data for each site
Ecoli.summ <- group_by(Ecoli.data.filt, Site) %>%
  dplyr::summarise(across(1:3, first),
                   across(5, list(mean = function(x) mean(x, na.rm = T), 
                                  max = function(x) max(x, na.rm = T), 
                                  min = function(x) min(x, na.rm = T))))
Ecoli.month <- group_by(Ecoli.data.filt, Site, Month) %>%
  dplyr::summarise(monthly_range = (max(Ecoli, na.rm=TRUE) - min(Ecoli, na.rm=TRUE)))
Ecoli.month2 <- group_by(Ecoli.month, Site) %>%
  dplyr::summarise(Ecoli_mean_monthly_range = mean(monthly_range, na.rm = T))

Ecoli.year <- group_by(Ecoli.data.filt, Site, Year) %>%
  dplyr::summarise(yearly_range = (max(Ecoli, na.rm=TRUE) - min(Ecoli, na.rm=TRUE)))
Ecoli.year2 <- group_by(Ecoli.year, Site) %>%
  dplyr::summarise(Ecoli_mean_yearly_range = mean(yearly_range, na.rm = T))

Ecoli.site <-list(Ecoli.summ, Ecoli.month2, Ecoli.year2) %>% reduce(inner_join, by = "Site")
Ecoli.site$Lat<-as.numeric(Ecoli.site$Lat)
Ecoli.site$Lon<-as.numeric(Ecoli.site$Lon)
Ecoli.site$Month <- NULL
write.csv(Ecoli.site, "Ecoli_STX.csv")
Ecoli.site<-read.csv("Ecoli_STX.csv")[,-1]

Ecoli.data.filt.year=Ecoli.data.filt %>% 
  select(Site, Year) %>% 
  unique() %>% 
  group_by(Year) %>% 
  summarize(Ecolin = n())
Ecoli.data.filt.month=Ecoli.data.filt %>% 
  select(Site, Month) %>% 
  unique() %>% 
  group_by(Month) %>% 
  summarize(Ecolin = n())


###-------- import Nitrogen -----#########
setwd("~/Spatial_interpolation")
Nitrogen0<-read.csv("WQD_Nitrogen.csv")
Nitrogen1<- Nitrogen0 %>%
  separate(ActivityStartDate, c("Month", "Day", "Year"), "/")
Nitrogen2<-data.frame(cbind(Nitrogen1$MonitoringLocationName, Nitrogen1$ActivityLocation.LatitudeMeasure, Nitrogen1$ActivityLocation.LongitudeMeasure, Nitrogen1$Month, Nitrogen1$Year, Nitrogen1$ResultMeasureValue))
names(Nitrogen2)<-c("Site", "Lat", "Lon", "Month", "Year", "Nitrogen")
Nitrogen2 <- Nitrogen2[Nitrogen2$Nitrogen != 0, ]
Nitrogen2 <- with(Nitrogen2, Nitrogen2[!(Nitrogen == "" | is.na(Nitrogen)), ])
Nitrogen2 <- with(Nitrogen2, Nitrogen2[!(Nitrogen == "<1" | is.na(Nitrogen)), ])
Nitrogen2 <- with(Nitrogen2, Nitrogen2[!(Nitrogen == "<10" | is.na(Nitrogen)), ])
Nitrogen<-Nitrogen2 %>% drop_na()
total.samples <- group_by(Nitrogen, Site) %>%
  dplyr::summarise(n = n())
ggplot(total.samples)+
  aes(x = n) +
  geom_histogram(bins = 30, fill = "#0c4c8a") +
  theme_minimal()
Nitrogen.data.filt1 <- filter(Nitrogen, Site %in% total.samples[total.samples$n > 10,]$Site)
Nitrogen.data.filt <- filter(Nitrogen.data.filt1, Year > 17 & Year < 23)

Nitrogen.data.filt$Nitrogen<-as.numeric(Nitrogen.data.filt$Nitrogen)
Nitrogen.data.filt <- Nitrogen.data.filt[!is.na(as.numeric(Nitrogen.data.filt$Nitrogen)), ]   
n_distinct(Nitrogen.data.filt$Site)
# Summarize data for each site
Nitrogen.summ <- group_by(Nitrogen.data.filt, Site) %>%
  dplyr::summarise(across(1:3, first),
                   across(5, list(mean = function(x) mean(x, na.rm = T), 
                                  max = function(x) max(x, na.rm = T), 
                                  min = function(x) min(x, na.rm = T))))
Nitrogen.month <- group_by(Nitrogen.data.filt, Site, Month) %>%
  dplyr::summarise(monthly_range = (max(Nitrogen, na.rm=TRUE) - min(Nitrogen, na.rm=TRUE)))
Nitrogen.month2 <- group_by(Nitrogen.month, Site) %>%
  dplyr::summarise(Nitrogen_mean_monthly_range = mean(monthly_range, na.rm = T))

Nitrogen.year <- group_by(Nitrogen.data.filt, Site, Year) %>%
  dplyr::summarise(yearly_range = (max(Nitrogen, na.rm=TRUE) - min(Nitrogen, na.rm=TRUE)))
Nitrogen.year2 <- group_by(Nitrogen.year, Site) %>%
  dplyr::summarise(Nitrogen_mean_yearly_range = mean(yearly_range, na.rm = T))

Nitrogen.site <-list(Nitrogen.summ, Nitrogen.month2, Nitrogen.year2) %>% reduce(inner_join, by = "Site")
Nitrogen.site$Lat<-as.numeric(Nitrogen.site$Lat)
Nitrogen.site$Lon<-as.numeric(Nitrogen.site$Lon)
Nitrogen.site$Month <- NULL
write.csv(Nitrogen.site, "Nitrogen_STX.csv")
Nitrogen.site<-read.csv("Nitrogen_STX.csv")[,-1]

Nitrogen.data.filt.year=Nitrogen.data.filt %>% 
  select(Site, Year) %>% 
  unique() %>% 
  group_by(Year) %>% 
  summarize(Nitrogenn = n())
Nitrogen.data.filt.month=Nitrogen.data.filt %>% 
  select(Site, Month) %>% 
  unique() %>% 
  group_by(Month) %>% 
  summarize(Nitrogenn = n())


###-------- import pH -----#########

setwd("~/Spatial_interpolation")
pH0<-read.csv("WQD_pH.csv")
pH1<- pH0 %>%
  separate(ActivityStartDate, c("Month", "Day", "Year"), "/")
pH2<-data.frame(cbind(pH1$MonitoringLocationName, pH1$ActivityLocation.LatitudeMeasure, pH1$ActivityLocation.LongitudeMeasure, pH1$Month, pH1$Year, pH1$ResultMeasureValue))
names(pH2)<-c("Site", "Lat", "Lon", "Month", "Year", "pH")
pH2 <- pH2[pH2$pH != 0, ]
pH2 <- with(pH2, pH2[!(pH == "" | is.na(pH)), ])
pH2 <- with(pH2, pH2[!(pH == "<1" | is.na(pH)), ])
pH2 <- with(pH2, pH2[!(pH == "<10" | is.na(pH)), ])
pH2$pH <- as.numeric(pH2$pH)
pH2 <- with(pH2, pH2[ (pH < 14 | is.na(pH)), ])
#pH2 <- with(pH2, pH2[!(pH < 0 | is.na(pH)), ])
pH<-pH2 %>% drop_na()
total.samples <- group_by(pH, Site) %>%
  dplyr::summarise(n = n())
ggplot(total.samples)+
  aes(x = n) +
  geom_histogram(bins = 30, fill = "#0c4c8a") +
  theme_minimal()
pH.data.filt1 <- filter(pH, Site %in% total.samples[total.samples$n > 50,]$Site)
pH.data.filt <- filter(pH.data.filt1, Year > 11 & Year < 23)

pH.data.filt$pH<-as.numeric(pH.data.filt$pH)
pH.data.filt <- pH.data.filt[!is.na(as.numeric(pH.data.filt$pH)), ]   
n_distinct(pH.data.filt$Site)
# Summarize data for each site
pH.summ <- group_by(pH.data.filt, Site) %>%
  dplyr::summarise(across(1:3, first),
                   across(5, list(mean = function(x) mean(x, na.rm = T), 
                                  max = function(x) max(x, na.rm = T), 
                                  min = function(x) min(x, na.rm = T))))
pH.month <- group_by(pH.data.filt, Site, Month) %>%
  dplyr::summarise(monthly_range = (max(pH, na.rm=TRUE) - min(pH, na.rm=TRUE)))
pH.month2 <- group_by(pH.month, Site) %>%
  dplyr::summarise(pH_mean_monthly_range = mean(monthly_range, na.rm = T))

pH.year <- group_by(pH.data.filt, Site, Year) %>%
  dplyr::summarise(yearly_range = (max(pH, na.rm=TRUE) - min(pH, na.rm=TRUE)))
pH.year2 <- group_by(pH.year, Site) %>%
  dplyr::summarise(pH_mean_yearly_range = mean(yearly_range, na.rm = T))

pH.site <-list(pH.summ, pH.month2, pH.year2) %>% reduce(inner_join, by = "Site")
pH.site$Lat<-as.numeric(pH.site$Lat)
pH.site$Lon<-as.numeric(pH.site$Lon)
pH.site$Month <- NULL
write.csv(pH.site, "pH_STX.csv")
pH.site<-read.csv("pH_STX.csv")[,-1]

pH.data.filt.year=pH.data.filt %>% 
  select(Site, Year) %>% 
  unique() %>% 
  group_by(Year) %>% 
  summarize(pHn = n())
pH.data.filt.month=pH.data.filt %>% 
  select(Site, Month) %>% 
  unique() %>% 
  group_by(Month) %>% 
  summarize(pHn = n())


###-------- import Phosphorus -----#########

setwd("~/Spatial_interpolation")
Phosphorus0<-read.csv("WQD_Phosphorus.csv")
Phosphorus1<- Phosphorus0 %>%
  separate(ActivityStartDate, c("Month", "Day", "Year"), "/")
Phosphorus2<-data.frame(cbind(Phosphorus1$MonitoringLocationName, Phosphorus1$ActivityLocation.LatitudeMeasure, Phosphorus1$ActivityLocation.LongitudeMeasure, Phosphorus1$Month, Phosphorus1$Year, Phosphorus1$ResultMeasureValue))
names(Phosphorus2)<-c("Site", "Lat", "Lon", "Month", "Year", "Phosphorus")
Phosphorus2 <- Phosphorus2[Phosphorus2$Phosphorus != 0, ]
Phosphorus2 <- with(Phosphorus2, Phosphorus2[!(Phosphorus == "" | is.na(Phosphorus)), ])
Phosphorus2 <- with(Phosphorus2, Phosphorus2[!(Phosphorus == "<1" | is.na(Phosphorus)), ])
Phosphorus2 <- with(Phosphorus2, Phosphorus2[!(Phosphorus == "<10" | is.na(Phosphorus)), ])
Phosphorus<-Phosphorus2 %>% drop_na()
total.samples <- group_by(Phosphorus, Site) %>%
  dplyr::summarise(n = n())
ggplot(total.samples)+
  aes(x = n) +
  geom_histogram(bins = 30, fill = "#0c4c8a") +
  theme_minimal()
Phosphorus.data.filt1 <- filter(Phosphorus, Site %in% total.samples[total.samples$n > 30,]$Site)
Phosphorus.data.filt <- filter(Phosphorus.data.filt1, Year > 12 & Year < 23)

Phosphorus.data.filt$Phosphorus<-as.numeric(Phosphorus.data.filt$Phosphorus)
Phosphorus.data.filt <- Phosphorus.data.filt[!is.na(as.numeric(Phosphorus.data.filt$Phosphorus)), ]   
n_distinct(Phosphorus.data.filt$Site)
# Summarize data for each site
Phosphorus.summ <- group_by(Phosphorus.data.filt, Site) %>%
  dplyr::summarise(across(1:3, first),
                   across(5, list(mean = function(x) mean(x, na.rm = T), 
                                  max = function(x) max(x, na.rm = T), 
                                  min = function(x) min(x, na.rm = T))))
Phosphorus.month <- group_by(Phosphorus.data.filt, Site, Month) %>%
  dplyr::summarise(monthly_range = (max(Phosphorus, na.rm=TRUE) - min(Phosphorus, na.rm=TRUE)))
Phosphorus.month2 <- group_by(Phosphorus.month, Site) %>%
  dplyr::summarise(Phosphorus_mean_monthly_range = mean(monthly_range, na.rm = T))

Phosphorus.year <- group_by(Phosphorus.data.filt, Site, Year) %>%
  dplyr::summarise(yearly_range = (max(Phosphorus, na.rm=TRUE) - min(Phosphorus, na.rm=TRUE)))
Phosphorus.year2 <- group_by(Phosphorus.year, Site) %>%
  dplyr::summarise(Phosphorus_mean_yearly_range = mean(yearly_range, na.rm = T))

Phosphorus.site <-list(Phosphorus.summ, Phosphorus.month2, Phosphorus.year2) %>% reduce(inner_join, by = "Site")
Phosphorus.site$Lat<-as.numeric(Phosphorus.site$Lat)
Phosphorus.site$Lon<-as.numeric(Phosphorus.site$Lon)
Phosphorus.site$Month <- NULL
write.csv(Phosphorus.site, "Phosphorus_STX.csv")
Phosphorus.site<-read.csv("Phosphorus_STX.csv")[,-1]

Phosphorus.data.filt.year=Phosphorus.data.filt %>% 
  select(Site, Year) %>% 
  unique() %>% 
  group_by(Year) %>% 
  summarize(Phosphorusn = n())
Phosphorus.data.filt.month=Phosphorus.data.filt %>% 
  select(Site, Month) %>% 
  unique() %>% 
  group_by(Month) %>% 
  summarize(Phosphorusn = n())



###-------- import Secchi -----#########

setwd("~/Spatial_interpolation")
Secchi0<-read.csv("WQD_Secchi_Depth.csv")
Secchi1<- Secchi0 %>%
  separate(ActivityStartDate, c("Month", "Day", "Year"), "/")
Secchi2<-data.frame(cbind(Secchi1$MonitoringLocationName, Secchi1$ActivityLocation.LatitudeMeasure, Secchi1$ActivityLocation.LongitudeMeasure, Secchi1$Month, Secchi1$Year, Secchi1$ResultMeasureValue))
names(Secchi2)<-c("Site", "Lat", "Lon", "Month", "Year", "Secchi")
Secchi2 <- Secchi2[Secchi2$Secchi != 0, ]
Secchi2 <- with(Secchi2, Secchi2[!(Secchi == "" | is.na(Secchi)), ])
Secchi2 <- with(Secchi2, Secchi2[!(Secchi == "<1" | is.na(Secchi)), ])
Secchi2 <- with(Secchi2, Secchi2[!(Secchi == "<10" | is.na(Secchi)), ])
Secchi<-Secchi2 %>% drop_na()
total.samples <- group_by(Secchi, Site) %>%
  dplyr::summarise(n = n())
ggplot(total.samples)+
  aes(x = n) +
  geom_histogram(bins = 30, fill = "#0c4c8a") +
  theme_minimal()
Secchi.data.filt1 <- filter(Secchi, Site %in% total.samples[total.samples$n > 30,]$Site)
Secchi.data.filt <- filter(Secchi.data.filt1, Year > 12 & Year < 23)

Secchi.data.filt$Secchi<-as.numeric(Secchi.data.filt$Secchi)
Secchi.data.filt <- Secchi.data.filt[!is.na(as.numeric(Secchi.data.filt$Secchi)), ]   
n_distinct(Secchi.data.filt$Site)
# Summarize data for each site
Secchi.summ <- group_by(Secchi.data.filt, Site) %>%
  dplyr::summarise(across(1:3, first),
                   across(5, list(mean = function(x) mean(x, na.rm = T), 
                                  max = function(x) max(x, na.rm = T), 
                                  min = function(x) min(x, na.rm = T))))
Secchi.month <- group_by(Secchi.data.filt, Site, Month) %>%
  dplyr::summarise(monthly_range = (max(Secchi, na.rm=TRUE) - min(Secchi, na.rm=TRUE)))
Secchi.month2 <- group_by(Secchi.month, Site) %>%
  dplyr::summarise(Secchi_mean_monthly_range = mean(monthly_range, na.rm = T))

Secchi.year <- group_by(Secchi.data.filt, Site, Year) %>%
  dplyr::summarise(yearly_range = (max(Secchi, na.rm=TRUE) - min(Secchi, na.rm=TRUE)))
Secchi.year2 <- group_by(Secchi.year, Site) %>%
  dplyr::summarise(Secchi_mean_yearly_range = mean(yearly_range, na.rm = T))

Secchi.site <-list(Secchi.summ, Secchi.month2, Secchi.year2) %>% reduce(inner_join, by = "Site")
Secchi.site$Lat<-as.numeric(Secchi.site$Lat)
Secchi.site$Lon<-as.numeric(Secchi.site$Lon)
Secchi.site$Month <- NULL
write.csv(Secchi.site, "Secchi_STX.csv")
Secchi.site<-read.csv("Secchi_STX.csv")[,-1]

Secchi.data.filt.year=Secchi.data.filt %>% 
  select(Site, Year) %>% 
  unique() %>% 
  group_by(Year) %>% 
  summarize(Secchin = n())
Secchi.data.filt.month=Secchi.data.filt %>% 
  select(Site, Month) %>% 
  unique() %>% 
  group_by(Month) %>% 
  summarize(Secchin = n())


###-------- import Temperature -----#########

setwd("~/Spatial_interpolation")
Temp0<-read.csv("WQD_Temp.csv")
Temp1<- Temp0 %>%
  separate(ActivityStartDate, c("Month", "Day", "Year"), "/")
Temp2<-data.frame(cbind(Temp1$MonitoringLocationName, Temp1$ActivityLocation.LatitudeMeasure, Temp1$ActivityLocation.LongitudeMeasure, Temp1$Month, Temp1$Year, Temp1$ResultMeasureValue))
names(Temp2)<-c("Site", "Lat", "Lon", "Month", "Year", "Temp")
Temp2 <- Temp2[Temp2$Temp != 0, ]
Temp2 <- with(Temp2, Temp2[!(Temp == "" | is.na(Temp)), ])
Temp2 <- with(Temp2, Temp2[!(Temp == "<1" | is.na(Temp)), ])
Temp2 <- with(Temp2, Temp2[!(Temp == "<10" | is.na(Temp)), ])
Temp<-Temp2 %>% drop_na()
total.samples <- group_by(Temp, Site) %>%
  dplyr::summarise(n = n())
ggplot(total.samples)+
  aes(x = n) +
  geom_histogram(bins = 30, fill = "#0c4c8a") +
  theme_minimal()
Temp.data.filt1 <- filter(Temp, Site %in% total.samples[total.samples$n > 50,]$Site)
Temp.data.filt <- filter(Temp.data.filt1, Year > 00 & Year < 23)

Temp.data.filt$Temp<-as.numeric(Temp.data.filt$Temp)
Temp.data.filt <- Temp.data.filt[!is.na(as.numeric(Temp.data.filt$Temp)), ]   
n_distinct(Temp.data.filt$Site)
# Summarize data for each site
Temp.summ <- group_by(Temp.data.filt, Site) %>%
  dplyr::summarise(across(1:3, first),
                   across(5, list(mean = function(x) mean(x, na.rm = T), 
                                  max = function(x) max(x, na.rm = T), 
                                  min = function(x) min(x, na.rm = T))))
Temp.month <- group_by(Temp.data.filt, Site, Month) %>%
  dplyr::summarise(monthly_range = (max(Temp, na.rm=TRUE) - min(Temp, na.rm=TRUE)))
Temp.month2 <- group_by(Temp.month, Site) %>%
  dplyr::summarise(Temp_mean_monthly_range = mean(monthly_range, na.rm = T))

Temp.year <- group_by(Temp.data.filt, Site, Year) %>%
  dplyr::summarise(yearly_range = (max(Temp, na.rm=TRUE) - min(Temp, na.rm=TRUE)))
Temp.year2 <- group_by(Temp.year, Site) %>%
  dplyr::summarise(Temp_mean_yearly_range = mean(yearly_range, na.rm = T))

Temp.site <-list(Temp.summ, Temp.month2, Temp.year2) %>% reduce(inner_join, by = "Site")
Temp.site$Lat<-as.numeric(Temp.site$Lat)
Temp.site$Lon<-as.numeric(Temp.site$Lon)
Temp.site$Month <- NULL
write.csv(Temp.site, "Temp_STX.csv")
Temp.site<-read.csv("Temp_STX.csv")[,-1]

Temp.data.filt.year=Temp.data.filt %>% 
  select(Site, Year) %>% 
  unique() %>% 
  group_by(Year) %>% 
  summarize(Tempn = n())
Temp.data.filt.month=Temp.data.filt %>% 
  select(Site, Month) %>% 
  unique() %>% 
  group_by(Month) %>% 
  summarize(Tempn = n())


# Observations by month
month_n = full_join(Entero.data.filt.month, full_join(DO.data.filt.month, full_join(Ecoli.data.filt.month, full_join(Nitrogen.data.filt.month, full_join(pH.data.filt.month, full_join(Phosphorus.data.filt.month, full_join(Secchi.data.filt.month, Temp.data.filt.month)))))))
write.csv(month_n, "Month_n.csv")
month_n=read.csv("Month_n.csv")

# Observations by year
year_n = full_join(Entero.data.filt.year, full_join(DO.data.filt.year, full_join(Ecoli.data.filt.year, full_join(Nitrogen.data.filt.year, full_join(pH.data.filt.year, full_join(Phosphorus.data.filt.year, full_join(Secchi.data.filt.year, Temp.data.filt.year)))))))
write.csv(year_n, "Year_n.csv")
year_n=read.csv("Year_n.csv")

# Import shoreline map
library(maptools)
setwd("~/Documents/SERC/") 
if (!rgeosStatus()) gpclibPermit()
gshhs.f.b <- "gshhg-bin-2.3.6/gshhs_f.b"
sf1 <- getRgshhsMap(gshhs.f.b, xlim = c(-65,-64.5), ylim = c(17.6,17.9)) %>%
  fortify()

# Example plot
ggplot()+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.title = element_blank(), plot.title = element_text(hjust = 0.5))+
  geom_polygon(data = sf1, aes(x=long, y = lat, group = group), fill = "grey70", color='black', lwd = 0.1)+
  geom_point(data=Entero.site, aes(x=Lon, y=Lat, color = Entero_mean))+
  scale_color_viridis_c(na.value = NA)+
  coord_fixed(xlim = c(-65,-64.5), ylim = c(17.6,17.9), expand = 0)



###-------- Spatial interpolation -----#########

# Create convex hull around sampling sites
testpts<-Entero.site[,2:3]
maphull.index<-chull(testpts)
maphull.index <- as.integer(maphull.index)
map.hull <- Entero.site[maphull.index, 3:2]
ggplot()+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.title = element_blank(), plot.title = element_text(hjust = 0.5))+
  geom_polygon(data = map.hull, aes(x=Lon, y = Lat), fill = NA, color='black', lwd = 0.5)+
  geom_polygon(data = sf1, aes(x=long, y = lat, group = group), fill = "grey70", color='black', lwd = 0.1)+
  geom_point(data=Entero.site, aes(x=Lon, y=Lat))+
  scale_color_viridis_c(na.value = NA)+
  coord_fixed(xlim = c(-65,-64.5), ylim = c(17.6,17.9), expand = 0)

# Stretch the hull by generating a circle of points that surround the vertices
n <- 6 # number of points you want on the circle
pts.circle <- t(sapply(1:n,function(r)c(cos(2*r*pi/n),sin(2*r*pi/n))))
buffer.data <- list()
for(i in 1:nrow(map.hull)){
  # The higher the denominator, the smaller the circle/buffer
  buffer.data[[i]] <- data.frame(x = map.hull$Lon[i] + pts.circle[,1]/20, 
                                 y = map.hull$Lat[i] + pts.circle[,2]/20)
}
buffer.data.long <- do.call(rbind, buffer.data)
# Generate a new concave hull using the outside of the circles as the reference
map.hull.stretch <- data.frame(concaveman::concaveman(as.matrix(buffer.data.long), concavity = 1, length_threshold = 0.3))
#map.hull.stretch <- read.table('map_hull_stretch.txt', sep = " ", header = FALSE)
map.hull.stretch <- as.data.frame(map.hull.stretch)
ggplot()+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.title = element_blank(), plot.title = element_text(hjust = 0.5))+
  geom_polygon(data = map.hull.stretch, aes(x=V1, y = V2), fill = NA, color='black', lwd = 0.5)+
  geom_polygon(data = sf1, aes(x=long, y = lat, group = group), fill = "grey70", color='black', lwd = 0.1)+
  geom_point(data=Entero.site, aes(x=Lon, y=Lat))+
  scale_color_viridis_c(na.value = NA)+
  coord_fixed(xlim = c(-65,-64.5), ylim = c(17.6,17.9), expand = 0)

# Convert hull to Spatial object
library(raster)
hull.poly <- SpatialPolygons(list(Polygons(list(Polygon(cbind(map.hull.stretch$V1, map.hull.stretch$V2))), ID=1)), proj4string = CRS("+proj=longlat +datum=NAD83"))
hull.poly.df <- SpatialPolygonsDataFrame(hull.poly, data=data.frame(ID=1))
hull.grid <- raster(hull.poly, res = 1/250)



###-------- Select coral sampling sites for extracting data -----#########
setwd("~/Data")
site.locations <- read.csv('sitesSamples.csv', header = T)

# View collection sites on the map
ggplot()+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.title = element_blank(), plot.title = element_text(hjust = 0.5))+
  geom_polygon(data = sf1, aes(x=long, y = lat, group = group), fill = "grey70", color='black', lwd = 0.1)+
  geom_point(data = site.locations, aes(x = Longitude, y = Latitude), color = 'red', size = 0.5)+
  scale_fill_viridis_c(na.value = NA)+
  coord_fixed(xlim = c(-65,-64.5), ylim = c(17.6,17.9), expand = 0)

site.coords <- SpatialPoints(site.locations[,c('Longitude', 'Latitude')], proj4string = CRS("+proj=longlat +datum=NAD83"))
site.locations.shifted <- vector()
for (each.site in 1:nrow(site.locations)){
  # geodesic distances
  distances <- fields::rdist.earth(coordinates(site.coords[each.site]), coordinates(krige.mask), miles=FALSE)
  # coordinates of the pixel that is at shorter distance 
  shifted.coordinates<-coordinates(krige.mask)[which.min(distances),]
  site.locations.shifted <- rbind(site.locations.shifted, shifted.coordinates)
}
new.locations <- data.frame(site.locations.shifted) %>%
  cbind(site.locations, .)



###-------- Convert sampling data to Spatial object -----#########

# Choose one variable at a time to send through kriging loop
Temp.only<-Temp.site[,2:8]
wq.sp <- SpatialPoints(Temp.only[,2:1], proj4string = CRS("+proj=longlat +datum=NAD83"))
wq.spdf <- SpatialPointsDataFrame(wq.sp, as.data.frame(Temp.only))

DO.only<-DO.site[,2:8]
wq.sp <- SpatialPoints(DO.only[,2:1], proj4string = CRS("+proj=longlat +datum=NAD83"))
wq.spdf <- SpatialPointsDataFrame(wq.sp, as.data.frame(DO.only))

Ecoli.only<-Ecoli.site[,2:8]
Ecoli.only$Ecoli_min=NULL
wq.sp <- SpatialPoints(Ecoli.only[,2:1], proj4string = CRS("+proj=longlat +datum=NAD83"))
wq.spdf <- SpatialPointsDataFrame(wq.sp, as.data.frame(Ecoli.only))

Entero.only<-Entero.site[,2:9]
wq.sp <- SpatialPoints(Entero.only[,2:1], proj4string = CRS("+proj=longlat +datum=NAD83"))
wq.spdf <- SpatialPointsDataFrame(wq.sp, as.data.frame(Entero.only))

Nitrogen.only<-Nitrogen.site[,2:8]
wq.sp <- SpatialPoints(Nitrogen.only[,2:1], proj4string = CRS("+proj=longlat +datum=NAD83"))
wq.spdf <- SpatialPointsDataFrame(wq.sp, as.data.frame(Nitrogen.only))

pH.only<-pH.site[,2:8]
wq.sp <- SpatialPoints(pH.only[,2:1], proj4string = CRS("+proj=longlat +datum=NAD83"))
wq.spdf <- SpatialPointsDataFrame(wq.sp, as.data.frame(pH.only))

Phosphorus.only<-Phosphorus.site[,2:8]
wq.sp <- SpatialPoints(Phosphorus.only[,2:1], proj4string = CRS("+proj=longlat +datum=NAD83"))
wq.spdf <- SpatialPointsDataFrame(wq.sp, as.data.frame(Phosphorus.only))

Secchi.only<-Secchi.site[,2:8]
wq.sp <- SpatialPoints(Secchi.only[,2:1], proj4string = CRS("+proj=longlat +datum=NAD83"))
wq.spdf <- SpatialPointsDataFrame(wq.sp, as.data.frame(Secchi.only))



###-------- Kriging loop -----#########

setwd("~/Spatial_interpolation/Output")
library(automap)

data.out = data.frame(new.locations[,1:3])
for(c in names(wq.spdf)) { 
  interp.var <- c
  test <- na.omit(wq.spdf@data[,c("Lat", "Lon", c)])
  #make test into a new spatial dataframe
  wq.sp.new <- SpatialPoints(test[,2:1], proj4string = CRS("+proj=utm +datum=NAD83"))
  wq.spdf.new <- SpatialPointsDataFrame(wq.sp.new, as.data.frame(test))
  
  # Auto-fit a model variogram
  variogram<-automap::autofitVariogram(get(interp.var)~1, wq.spdf.new)
  
  plot(variogram)
  spatial.points <- SpatialPoints(coords = coordinates(hull.grid), 
                                  proj4string = CRS("+proj=utm +datum=NAD83"))
  spatial.grid <- SpatialPixels(points = spatial.points)
  kriging_result <- automap::autoKrige(get(interp.var)~1, wq.spdf.new, spatial.grid)
  # save variogram
  pdf(file = paste0("variogram_", interp.var, ".pdf"))
  plot(variogram)
  dev.off()
  
  # Create Rasterbrick object and mask to St. Croix polygon
  krige.brick <- brick(kriging_result$krige_output)
  krige.mask <- mask(krige.brick, hull.poly)
  names(krige.mask) <- c('prediction', 'variance', 'stdev')
  # save prediction map
  pdf(file = paste0("prediction_", interp.var, ".pdf"))
  plot(krige.mask)
  dev.off()
  
  krige.predict <- cbind(as.data.frame(krige.mask$prediction), as.data.frame(spatial.grid))
  ggpredict <- ggplot()+
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.title = element_blank(), plot.title = element_text(hjust = 0.5))+
    geom_tile(data = krige.predict, aes(x = x, y = y, fill = prediction, col = prediction))+
    geom_polygon(data = sf1, aes(x=long, y = lat, group = group), fill = "grey70", color='black', lwd = 0.1)+
    scale_fill_viridis_c(na.value = NA)+
    scale_color_viridis_c(na.value = NA)+
    coord_fixed(xlim = c(-65,-64.5), ylim = c(17.6,17.9), expand = 0)
  # # save ggplot prediction map
  ggsave(file = paste0("ggpredict_", interp.var, ".pdf"))
  plot(ggpredict)
  dev.off()
  
  # save a raster
  krig.fixed<-data.frame(cbind(krige.predict$x, krige.predict$y, krige.predict$prediction))
  r = rasterFromXYZ(krig.fixed)
  plot(r)
  writeRaster(r, file = paste0(interp.var), format="raster", overwrite=TRUE)
  
  #extract data
  extract.data <- data.frame(raster::extract(krige.mask, new.locations[,4:5])[,1])
  colnames(extract.data) <- interp.var
  data.out <- cbind(data.out, extract.data)
  data.out$Lat=NULL
  data.out$Lon=NULL
  
}

write.csv(data.out, "STXenv_allvars.csv")








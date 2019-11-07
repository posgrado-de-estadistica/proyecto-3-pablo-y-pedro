# Unión tablas

library(RDS)
library(dplyr)
library(sp)
library(stringr)

datos_encabezado <- readRDS("Datos/datos_encabezado.Rds")
datos_tabla1 <- readRDS("Datos/datos_tabla1.Rds")

datos_encabezado <- datos_encabezado %>% 
    mutate(ch1= as.numeric(substr(LAT, 1, 2)),
           chm1 = as.numeric(substr(LAT, 6,8)),
           ch2= as.numeric(substr(LON, 1, 2)),
           chm2 = as.numeric(substr(LON, 6,8))) %>% 
    mutate(chm1= ifelse(is.na(chm1),59,chm1)) %>% 
    mutate(LAT= ch1+ chm1/60,
           LON= (ch2+ chm2/60)*-1,
           ALT= as.numeric(str_extract(datos_encabezado$ALT,"\\d+"))
           ) %>%
    dplyr::select(-ch1, -ch2, -chm1,- chm2)


datos_tabla1 <- datos_tabla1 %>% 
    rename("ELEMENTOS"= "V1",
           "PERIODO1"= "V2",
           "PERIODO2"= "V3",
           "ENE"= "V4",
           "FEB"= "V5",
           "MAR"= "V6",
           "ABR"= "V7",
           "MAY"= "V8",
           "JUN"= "V9",
           "JUL"= "V10",
           "AGO"= "V11",
           "SET"= "V12",
           "OCT"= "V13",
           "NOV"= "V14",
           "DIC"= "V15",
           "COD_ESTACION"= "archivo") %>% 
    mutate(COD_ESTACION= gsub("_.*","",COD_ESTACION))

datos_finales <- datos_encabezado %>% left_join(datos_tabla1)

saveRDS(datos_finales, file="Datos/datos_finales.Rds")

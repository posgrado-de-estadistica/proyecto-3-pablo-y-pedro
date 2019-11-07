# Extracción PDF

library(dplyr)
library(pdftools)
library(readr)    
library(data.table)
library(stringr)
library(RDS)

# Path pdf

file_vector <- list.files(path = "PDF/")
file_vector %>% head()

# Prueba

pdf_text(paste0("PDF/",file_vector[[2]])) %>% 
    strsplit(split = "\n")


# Tabla 1

tablas <- vector("list", NROW(file_vector))

for (i in seq_along(file_vector)) {
    leer <- pdf_text(paste0("PDF/",file_vector[[i]])) %>% 
        strsplit(split = "\n")
    validos <-grepl("[12]\\d{3}\\s{5}",leer[[1]])
    
    a<-paste(leer[[1]][validos])
    b<- str_trim(a)
    b2<-  gsub(".*(VEL.)","VELOCIDAD",b)
    b2<-  gsub(".*(TANQUE)","EVAPORACION",b2)
    b2<-  gsub(".*(SOLAR)","BRILLO",b2)
    c<- strsplit(b2,"[[:blank:]]")
    c2<-lapply(c, function(x) x[!x %in% ""])
    c3 <- rbindlist(lapply(c2, function(x) data.table((t(x)))),fill = TRUE) %>% 
        dplyr::select(V1:V15) %>% mutate(archivo= !!file_vector[[i]])
    tablas[[i]] <- c3
}

datos <- bind_rows(tablas)  

saveRDS(datos, file = "Datos/datos_tabla1.Rds")

# Tabla 2

tablas2 <- vector("list", NROW(file_vector))

for (i in seq_along(file_vector)) {
    leer <- pdf_text(paste0("PDF/",file_vector[[i]])) %>% 
        strsplit(split = "\n")
    validos <-grepl("(PREDOMINANTE)",leer[[1]])
    
    a<-c(paste(leer[[1]][which(validos)]),paste(leer[[1]][which(validos)+1]))
    b<- str_trim(a)
    c <- substring(b,40,200)
    c2<- strsplit(c,"[[:blank:]]")
    c3<-lapply(c2, function(x) x[!x %in% ""])
    c4 <- rbindlist(lapply(c3, function(x) data.table((t(x)))),fill = TRUE) %>% 
        dplyr::select(V1:V12) %>% mutate(archivo= !!file_vector[[i]])
    tablas2[[i]] <- c4
}

datos2 <- bind_rows(tablas2)  

saveRDS(datos2, file = "Datos/datos_tabla2.Rds")

# LON-LAT

#LATITUD<- substring(b, 68, 75) 
#LONGITUD<- substring(b, 90, 98)
#ALTITUD<- substring(b, 114, 124)


tablas3 <- vector("list", NROW(file_vector))

for (i in seq_along(file_vector)) {
    leer <- pdf_text(paste0("PDF/",file_vector[[i]])) %>% 
        strsplit(split = "\n")
    validos <-leer[[1]][[5]]
    
    b<- str_trim(validos)
    
    c <- data.frame(
    COD_ESTACION= gsub("_.*","",file_vector[[i]]),
    LAT= gsub(".*(Latitud)","",b) %>% substring(., 3, 14) %>% str_trim(.),
    LON= gsub(".*(Longitud)","",b)%>% substring(., 3, 15)%>% str_trim(.),
    ALT= gsub(".*(Altitud)","",b)%>% substring(., 3, 15)%>% str_trim(.)
    )
    tablas3[[i]] <- c
}

datos3 <- bind_rows(tablas3)  

saveRDS(datos3, file = "Datos/datos_encabezado.Rds")

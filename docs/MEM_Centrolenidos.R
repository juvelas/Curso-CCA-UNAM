

library(raster)
library(sp)
library(maptools)
library(rgdal)
library(dismo)
library(XML)
library(maps)
library(letsR)
library(plyr)
library(groupdata2)
library(rgeos)
library(classInt)
library(Metrics)
library(RColorBrewer)
library(rcompanion)
library(plotKML)
library(lmtest)
library(spdep)
library(lme4)
library(reshape2)
library(usdm)


# cargar shapefile con ranges de todas las especies
crs.geo <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")  # geographical, datum WGS84

mapas <-readOGR("./Shapefiles/Centrolenidae_IUCN.shp")
projection(mapas)<-crs.geo


# generar un mapa de riqueza y una matriz de presencia-ausencia (PAM) usando los mapas
PAM_centro <-lets.presab(mapas, xmn=-100.4, xmx=-37.5, ymn=-29.3, ymx=19.6, resol=1,
                         remove.cells=TRUE, remove.sp=TRUE, show.matrix=FALSE,
                         crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"), cover=0, presence=NULL,
                         origin=NULL, seasonal=NULL, count=FALSE)
#
plot(PAM_centro)

#
PAM_centro

#
str(PAM_centro)

# extraer raster
rr <- PAM_centro$Richness_Raster

# cargar un mapa de países
am <-readOGR("./Shapefiles/america_completo.shp")
projection(am)<-crs.geo

# graficar
plot(rr, main="Riqueza de especies")
plot(am, add=T)


# Subir capas bioclimáticas actuales y futuras
bios.act <- stack(list.files(path="./Bioclimaticos/actual/", pattern = "tif", full.names = T))
bios.fut <- stack(list.files(path="./Bioclimaticos/rcp85_2070/bc85bi70/", pattern = "tif", full.names = T))

projection(bios.act)<-crs.geo
projection(bios.fut)<-crs.geo


# Extraer datos de variables ambientales para cada celda con datos de riqueza
PAM_centro_bios <- lets.addvar(PAM_centro, bios.act, fun=mean,  onlyvar = TRUE)
PAM_centro_biosF <- lets.addvar(PAM_centro, bios.fut, fun=mean,  onlyvar = TRUE)

# Preparar los conjuntos de datos de riqueza de especies y variables bioclimáticas
PAM_DF <- as.data.frame(PAM_centro$Presence_and_Absence_Matrix)
a <- rowSums(PAM_DF[,c(-1:-2)])
DF_SR <- PAM_DF[,c(1:2)]
DF_SR["Riqueza"] <- a

DF_SR_Bios <- round(cbind(DF_SR, PAM_centro_bios),1)
DF_SR_BiosFut <- round(cbind(DF_SR, PAM_centro_biosF),1)

# Construir el modelo de regresión (calibrar, ajustar, etc.)

# Generamos una partición aleatoria de los datos
partitions <- partition(DF_SR_Bios, p=0.7)
train <- partitions[[1]]
test <- partitions[[2]]
trainDF <- as.data.frame(train)
testDF <- as.data.frame(test)

# Ajustar el modelo (Ordinary Least Squares)
m.ols <- lm(Riqueza~., data=trainDF[,c(-1,-2)])
summary(m.ols)
# Generar la predicción
pred.ols <- predict(m.ols, testDF)

par(mfrow=c(1,1))
hist(testDF$Riqueza, main="Observaciones")
hist(pred.ols, main="Predicciones (fuera de la muestra)")


###


# loop para generar 100 particiones aleatorias y hacer validación cruzada
nsim <- 100
mfit_ols <- list() # guarda cada modelo ajustado
rsq_ols <- list() # guarda cada R2 calculado
rmse_ols <- list() # guarda la raíz del error cuadrático medio RMSE
pred_ols <- list() # predicciones fuera de la muestra

for (i in 1:nsim){
#
partitions <- partition(DF_SR_Bios, p=0.7)
train <- partitions[[1]]
test <- partitions[[2]]
trainDF <- as.data.frame(train)
testDF <- as.data.frame(test)

# model fit
m.ols <- lm(Riqueza~., data=trainDF[,c(-1:-2)])
mfit_ols[[i]] <- m.ols

# predictions
pred.ols <- predict(m.ols, testDF[,c(-1:-2)])
pred_ols[[i]] <- pred.glm

# metrics
rss <- sum((pred.ols - testDF$Riqueza) ^ 2)  ## suma de cuadrados residuales
tss <- sum((testDF$Riqueza - mean(testDF$Riqueza)) ^ 2)  ## suma total de cuadrados
rsq_ols[[i]] <- 1 - rss/tss # R2
rmse_ols[[i]] <- rmse(testDF$Riqueza, pred.ols) #  raíz del error cuadrático medio RMSE
}

# Histograma de los 100 valores de R2 y RMSE
bb <- cbind(as.numeric(rsq_ols), as.numeric(rmse_ols))
colnames(bb) <- c("R2", "RMSE")
v.ols <- as.data.frame(bb)

par(mfrow=c(1,2))
hist(v.ols$R2, main="R2")
hist(v.ols$RMSE, main="RMSE")


# GLM
nsim <- 100
# to run models across a sample of 100 spatial folds
mfit_glm <- list() # model fitted
rsq_glm <- list() # R squared
rmse_glm <- list() # Root Mean Squared Error
pred_glm <- list() # predicciones fuera de la muestra

for (i in 1:nsim){
#
partitions <- partition(DF_SR_Bios, p=0.7)
train <- partitions[[1]]
test <- partitions[[2]]
trainDF <- as.data.frame(train)
testDF <- as.data.frame(test)

# model fit
m.glm <- glm(Riqueza~., data = trainDF[,c(-1:-2)], family=poisson())
mfit_glm[[i]] <- m.glm

# predictions
pred.glm <- predict(m.glm, newdata=testDF[,c(-1:-2)], type="response", se.fit=FALSE)
pred_glm[[i]] <- pred.glm

# metrics
rss <- sum((pred.glm - testDF$Riqueza) ^ 2)  ## suma de cuadrados residuales
tss <- sum((testDF$Riqueza - mean(testDF$Riqueza)) ^ 2)  ## suma total de cuadrados
rsq_glm[[i]] <- 1 - rss/tss # R2
rmse_glm[[i]] <- rmse(testDF$Riqueza, pred.glm) #  raíz del error cuadrático medio RMSE
}


# Histograma de los 100 valores de R2 y RMSE
bb <- cbind(as.numeric(rsq_glm), as.numeric(rmse_glm))
colnames(bb) <- c("R2", "RMSE")
v.glm <- as.data.frame(bb)


par(mfrow=c(1,2))
hist(v.ols$R2, main="R2 OLS")
hist(v.ols$RMSE, main="RMSE OLS")

par(mfrow=c(1,2))
hist(v.glm$R2, main="R2 GLM")
hist(v.glm$RMSE, main="RMSE GLM")

#  Proyectar a futuro
options(scipen=999)


m.ols <- lm(Riqueza~., data=DF_SR_Bios[,c(-1:-2)])
m.glm <- glm(Riqueza~., data=DF_SR_Bios[,c(-1:-2)], family=poisson())

pred.ols.act <- round(predict(m.ols, newdata=DF_SR_Bios[,c(-1:-2)], type="response", se.fit=FALSE),1)
pred.glm.act <- round(predict(m.glm, newdata=DF_SR_Bios[,c(-1:-2)], type="response", se.fit=FALSE),1)

pred.ols.fut <- round(predict(m.ols, newdata=DF_SR_BiosFut[,c(-1:-2)], type="response", se.fit=FALSE),1)
pred.glm.fut <- round(predict(m.glm, newdata=DF_SR_BiosFut[,c(-1:-2)], type="response", se.fit=FALSE),1)

# Comparar
pp.ols <- pred.ols.fut - pred.ols.act # negativo = pérdidas, positivo = ganancias
pp.glm <- pred.glm.fut - pred.glm.act # negativo = pérdidas, positivo = ganancias

# pegar las predicciones en el dataframe completo
DF_SR_Bios["pred_ols_act"] <- pred.ols.act
DF_SR_Bios["pred_ols_fut"] <- pred.ols.fut

DF_SR_Bios["pred_glm_act"] <- pred.glm.act
DF_SR_Bios["pred_glm_fut"] <- pred.glm.fut


DF_SR_Bios_sp <- DF_SR_Bios
DF_SR_Bios_sp <- plyr::rename(DF_SR_Bios_sp, replace=c("Longitude(x)"="x", "Latitude(y)"="y"))
coordinates(DF_SR_Bios_sp)<- ~ x + y

class(DF_SR_Bios_sp)

obs.r <- vect2rast(DF_SR_Bios_sp, fname = "Riqueza", cell.size=1)
obs.r  <- raster(obs.r)
projection(obs.r)<-crs.geo

ols.act.r <- vect2rast(DF_SR_Bios_sp, fname = "pred_ols_act", cell.size=1)
ols.act.r  <- raster(ols.act.r)
projection(ols.act.r)<-crs.geo

ols.fut.r <- vect2rast(DF_SR_Bios_sp, fname = "pred_ols_fut", cell.size=1)
ols.fut.r  <- raster(ols.fut.r)
projection(ols.fut.r)<-crs.geo

glm.act.r <- vect2rast(DF_SR_Bios_sp, fname = "pred_glm_act", cell.size=1)
glm.act.r  <- raster(glm.act.r)
projection(glm.act.r)<-crs.geo

glm.fut.r <- vect2rast(DF_SR_Bios_sp, fname = "pred_glm_fut", cell.size=1)
glm.fut.r  <- raster(glm.fut.r)
projection(glm.fut.r)<-crs.geo

dif.ols <- ols.fut.r - ols.act.r # negativo = pérdidas, positivo = ganancias
plot(dif.ols)

dif.glm <- glm.fut.r - glm.act.r # negativo = pérdidas, positivo = ganancias
plot(dif.glm)

perd.riq.ols <- reclassify(dif.ols, c(0, Inf, 0))
perd.riq.glm <- reclassify(dif.glm, c(0, Inf, 0))

gan.riq.ols <- reclassify(dif.ols, c(-Inf, 0, 0))
gan.riq.glm <- reclassify(dif.glm, c(-Inf, 0, 0))

# mapas
library(viridis)

par(mfrow=c(2,3))
plot(ols.act.r, main="OLS predicción", col=viridis(30))
plot(am, add=T, lwd=0.5)

plot(perd.riq.ols, main="OLS pérdidas", col=viridis(30))
plot(am, add=T, lwd=0.5)

plot(gan.riq.ols, main="OLS ganancias", col=viridis(30))
plot(am, add=T, lwd=0.5)

plot(glm.act.r, main="GLM predicción", col=viridis(30))
plot(am, add=T, lwd=0.5)

plot(perd.riq.glm, main="GLM pérdidas", col=viridis(30))
plot(am, add=T, lwd=0.5)

plot(gan.riq.glm, main="GLM ganancias", col=viridis(30))
plot(am, add=T, lwd=0.5)

# FIN






#

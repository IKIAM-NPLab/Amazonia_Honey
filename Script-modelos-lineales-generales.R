#Taller de análisis multivariado
#Estudiante: Jose A. Molina Mora
#Práctica: Modelos lineales generales

#librerias
library("psych")
library("corrplot")

#A. datos
datos<-read.csv("Datos_presion.csv")

#B. Estadisticas y graficos generales
summary(datos)
multi.hist(datos, global = FALSE)
corrplot.mixed(cor(datos, method = "pearson"))

#C. Modelo simple: Presion y Peso

colnames(datos)
simple<-lm(Presion_sistolica ~ Peso, data = datos) # lm(Y ~ X, data = )
names(simple)
summary(simple)

# Existe al menos un factor (peso, 0.000665) que se relaciona linealmente con la presion (p-value: 0.0006654)
# Adjusted R-squared:  0.2521: El peso explica un 25.21% de la variabilidad de la presion
# Curva ajuste: Presion = 66.6 + 0.96*Peso

# D. Graficas del modelo simple
pairs.panels(datos[,c("Peso","Presion_sistolica")])
plot(datos$Peso, datos$Presion_sistolica, main = "Regresion simple", pch=8, col="gray")
abline(simple, col="red")

# E. Prediccion
predict(simple, newdata = data.frame(Peso=c(80,70,25,0)))

# F. Modelo lineal multiple (con todos los factores)
multiple_todos<-lm(Presion_sistolica ~ . , data =datos)
names(multiple_todos)
summary(multiple_todos)
# Existe al menos un factor (peso, Años de vivir en region urbana) que se relaciona linealmente con la presion (p-value: 0.001768)
# Adjusted R-squared:  0.3839: El modelo explica un 38.39% de la variabilidad de la presion

# G. Modelo reducido
multiple_reducido<-lm(Presion_sistolica ~ Peso + Anos_region_urbana, data =datos)
summary(multiple_reducido)
# Existe al menos un factor (peso, Años de vivir en region urbana) que se relaciona linealmente con la presion (p-value: 0.00)
# Adjusted R-squared:  0.3886: El modelo explica un 38.86% de la variabilidad de la presion
# El modelo es: Presion = 50.31 + 1.35*Peso + -0.57*Anos_region_urbana

# H. Modelo con interaccion
multiple_interaccion<-lm(Presion_sistolica ~ Peso + Anos_region_urbana + Edad*Peso, data =datos)
summary(multiple_interaccion)
# Existe al menos un factor que se relaciona linealmente con la presion (p-value: 0.000595)
# Adjusted R-squared:  0.3625: El modelo explica un 36.25% de la variabilidad de la presion
# El R2 disminuye y varios factores se excluyen.

# I. Modelos multiples: Comparacion

# El mejor modelo es el reducido con los 2 factores: Presion = 50.31 + 1.35*Peso + -0.57*Anos_region_urbana

# J. Supuestos

par(mfrow=c(2,2))
plot(multiple_reducido)
par(mfrow=c(1,1))


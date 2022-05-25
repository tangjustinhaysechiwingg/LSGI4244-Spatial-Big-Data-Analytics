library()

install.packages("rgdal", repos = "http://cran.r-project.org")
install.packages("maps", repos = "http://cran.r-project.org")
install.packages("mapproj", repos = "http://cran.r-project.org")

library(rgdal) # shapefile operation
library(maps) # map operation
library(mapproj) # map projection

setwd("C:/Users/20016345d/Documents/R/Lab1/data")
getwd()

columbus.poly = readOGR(".", "columbus")
 class(columbus.poly)
 columbus.df = as.data.frame(columbus.poly) 

hist(columbus.df$CRIME,
main="The number of crime in Columbus in 1980",
xlab = "Number of Crime (per Thousand Households)",
ylab = "Number of Neighbourhood",
xlim = c(0,70),
ylim = c(0,10),
breaks = 20, 
col = "green",
border ="blue",
density = 50, 
angle = 30
)

 hist(columbus.df$HOVAL,
main="The Housing Value in Columbus in 1980",
xlab = "The Housing Value (in $1000)",
ylab = "Number of Neighbourhood",
xlim = c(0,100),
ylim = c(0,10),
breaks = 20, 
col = "red",
border ="green",
density = 60, 
angle = 30
)

 hist(columbus.df$INC,
main="The Household Income in Columbus in 1980",
xlab = "The Household Income (in $1000)",
ylab = "Number of Neighbourhood",
xlim = c(0,35),
ylim = c(0,7),
breaks = 20, 
col = "blue",
border ="purple",
density = 60, 
angle = 30
)


boxplot(columbus.df$CRIME,
main="The number of Crime in Columbus in 1980", 
xlab = "Number of Crime (per Thousand Households)",
notch = TRUE, 
varwidth = TRUE,
col = "green",
border ="blue",
horizontal = TRUE)

boxplot(columbus.df$HOVAL,
main="The Housing Value in Columbus in 1980", 
xlab = "The Housing Value (in $1000)",
notch = TRUE, 
varwidth = TRUE,
col = "red",
border ="purple",
horizontal = TRUE)

boxplot(columbus.df$INC,
main="The Household Income in Columbus in 1980", 
xlab = "The Household Income (in $1000)",
notch = TRUE, 
varwidth = TRUE,
col = "blue",
border ="black",
horizontal = TRUE)



qqnorm(columbus.df$CRIME,col = "blue",main="The QQ-Plot of the number of Crime in Columbus in 1980",ylab = "Number of Crime (per Thousand Households)",pch = 1, frame = FALSE)
qqline(columbus.df$CRIME,col = "red")


qqnorm(columbus.df$HOVAL,col = "red",main="The QQ-Plot of the Housing Value in Columbus in 1980",ylab = "The Housing Value (in $1000)")
qqline(columbus.df$HOVAL,col = "purple")


qqnorm(columbus.df$INC,col = "blue",main="The QQ-Plot of the Household Income in Columbus in 1980",ylab = "The Household Income (in $1000)")
qqline(columbus.df$INC,col = "red")


summary(columbus.df$CRIME)
summary(columbus.df$HOVAL)
summary(columbus.df$INC)



plot(columbus.df$INC,columbus.df$HOVAL,columbus.df$CRIME, xlab="Number of Crime", ylab = "The Housing Value",)
plot(columbus.df[,c(10,8,9)],pch=8, main="The Scatterplot Matrix for CRIME, HOVAL and INC Variables", col=rainbow(20),panel = panel.smooth)
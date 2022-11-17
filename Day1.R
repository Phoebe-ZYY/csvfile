data("women")
data("mtcars")

str(mtcars)
summary(mtcars)

mean(women$height)
mean(women$weight)

hist(mtcars$mpg,      
     main="Histogram for mpg",
     xlab="MPG",      
		 border="blue",      
		 col="green", 
		 las=1,      
		 breaks=5)

var(women$height)
sd(women$height)

getwd()
4*6
15/3
log(10)
log10(10)
4^3
5 < 8
5 == 8
5 != 8
1:250

"a" == "A"

x<-7
x+3
y <- 4
x + y
x <- x+3
x
x <- "Hello"
x

x <- c(1:10)
x <- 1:10
numeric_vector <- c(1, 2, 3)
character_vector <- c("a", "b", "c")

x<- 1:10
x*2
x < 10
y <- 5
x*y
x + y
z <- x + y
z
x <- x * y
x

data(mtcars)
mtcars
head(mtcars)
tail(mtcars)
str(mtcars)


scale(women$height)
women$height

pnorm(-0.68)
qnorm(0.9)

country <- read.csv("C:/me/course/Statistics/country.csv")
continent <- read.csv("C:/me/course/Statistics/continent.csv")
gdp <- read.csv("C:/me/course/Statistics/GDP.csv")

countryxlsx <- read_xlsx("C:/me/course/Statistics/country.xlsx")

head(gdp)
head(continent)
continentgdp <- merge(continent, gdp)
head(continentgdp)

#filter
oc_continent <- filter(continentgdp, Continent == "OC")
filter(continentgdp, Continent == "OC" & X2011 > 23000)
#select
select(continentgdp, X2006:X2012)
select(continentgdp, -(ISO2:Country))

#group
contiData <- group_by(continentgdp, Continent)
summarize(contiData, count = n(), GDP2012 = mean(X2012, na.rm = TRUE))

#convert to long format
#from wide to long
continentgdp_l <- reshape(continentgdp, varying = 4:13, direction = "long", sep = "")
head(continentgdp_l)

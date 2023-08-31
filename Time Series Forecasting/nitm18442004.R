#####################################
#installing packages
#####################################
install.packages("ggplot2")
install.packages("dplyr")
install.packages("forecast")
install.packages("tseries")
install.packages("ggseas")
install.packages("fpp2")

library(ggplot2)
library(ggseas)
library(dplyr)
library(forecast)
library(tseries)
library(fpp2)

#####################################
#Setting the Directory
#####################################
# getting the current directory
getwd()
# set the working directory to the folder containing the CSV file
setwd("/Users/sahelidutta/Documents/Statistics for Data Analytics/TABA/dataset")

#####################################
#Loading the data
#####################################

# read in the CSV file
monthly_data <- read.csv("nitm18442004.csv")
monthly_data
class(monthly_data)

# Convert the data to a time series object
ts_data <- ts(monthly_data$x, start = c(1844, 1), frequency = 12)
head(ts_data)
class(ts_data)

start(ts_data)
end(ts_data)
frequency(ts_data)

#####################################
#Priliminary Analysis of Time Series
#####################################

# Monthly data summary statistics
monthly_summary <- summary(monthly_data$x)
monthly_summary

# Create a time plot of the data
par(mfrow =c(1,1))
plot(ts_data, main = "Time plot of monthly average temperature in Armagh")

# Augmented Dickey-Fuller Test to check the data is stationary or not
adf.test(ts_data, alternative = "stationary")
#p value is less than 0.05, the data is stationary

# ADF test with differencing (Data is still stationary)
adf.test(diff(ts_data), alternative = "stationary")

# Visualizing based on century
plot.ts(ts_data, xlim=c(1844,1900), main = "Time plot of monthly average temperature from 1844 to 1900")
plot.ts(ts_data, xlim=c(1901,2004), main = "Time plot of monthly average temperature from 1901 to 2004")

# Checking Seasonality(Time series is already trend-stationary)
ggseasonplot(ts_data)+ggtitle("Seasonal Plot in monthly basis")+ylab("Frquency")

# Checking seaasonality in another subseries plot 
ggsubseriesplot(ts_data)+ggtitle("Seasonal Plot in monthly basis")+ylab("Frquency")

# Returning the seasonal cycle of the time series ts_data
cycle(ts_data)

# Compute mean and standard deviation of seasonal cycle
seasonal_mean <- mean(cycle(ts_data))
seasonal_sd <- sd(cycle(ts_data))

# Print out results
cat("Mean seasonal cycle:", seasonal_mean, "\n")
cat("Standard deviation of seasonal cycle:", seasonal_sd, "\n")

# Seasonal effects across 12 months (more fluctuation in june, july,august,spetember) 
boxplot(ts_data~cycle(ts_data))

# Calculate outliers
outliers <- boxplot.stats(ts_data)$out
outliers

# Create a histogram of the data
hist(ts_data, main = "Histogram of monthly average temperature in Armagh",
     xlab = "Temperature", ylab = "Frequency")

# Decompose trend, seasonal, remainder element in a plot
decomp = stl(ts_data, s.window = "periodic")
deseasonal <- seasadj(decomp)
plot(decomp)

# Plotting Correlogram (ACF, PACF)
par(mfrow =c(1,2))
acf(ts_data, lag.max = 25)
pacf(ts_data, lag.max = 25)

acf(ts_data, lag.max = 25, plot = FALSE)
pacf(ts_data, lag.max = 25, plot = FALSE)

#######################################
# Splitting the dataset
#######################################
# Split the data into training and test sets
train <- window(ts_data, end = c(2003, 12))
test <- window(ts_data, start = c(2004, 1))

#######################################
# Seasonal Naive for forcasting
# y_t = y_(t-s) +e_t
#######################################
fit <- snaive(train) # Residual sd: 1.705 
print(summary(fit))
checkresiduals(fit)

# Make forecasts for the test data
forecasts <- forecast(fit, h = length(test))
# Print the forecasts
print(summary(forecasts))
# Checking the residuals 
checkresiduals(forecasts)
# Evaluate the forecasts against the actual data for 2004
round(accuracy(forecasts, test),3)

# Plot the forecasts
par(mfrow =c(1,1))
plot(forecasts, main = "Seasonal Naive Forecasts")
# Add the actual data for comparison
lines(test, col = "red")
# Add a legend
legend("topleft", legend = c("Forecasts", "Actual", "Test"), col = c("black", "blue", "red"), lty = c(1,1))

#######################################
# Exponential Smoothning Model
#######################################
# Apply Holt's method with additive seasonality
fit_exp <- HoltWinters(train, seasonal = "additive")
print(summary(fit_exp))
checkresiduals(fit_exp)

# Make forecasts for the test data
forecasts_exp <- forecast(fit_exp, h = length(test))
print(summary(forecasts_exp))
# Checking the residuals 
checkresiduals(forecasts_exp)
# Evaluate the forecasts against the actual data for 2004
round(accuracy(forecasts_exp, test),3)

# Extract the residual standard error (RSE) value
rse_exp <- sqrt(mean(residuals(fit_exp)^2))
rse_exp # Residual SD: 1.250518

# Plot the forecasts
par(mfrow =c(1,1))
plot(forecasts_exp, xlab = "Year", ylab = "", main = "Holt Winters Forecast")
# Add the actual test data to the plot
lines(test, col = "red")
# Add a legend
legend("topleft", legend = c("Forecasts", "Actual", "Test"), col = c("black", "blue", "red"), lty = 1)

#########################################
# Automated Exponential Smoothning Model
########################################
# Fit an automated exponential smoothing model to the training data
fit_autoexp <- ets(train)
# Print the model summary
print(summary(fit_autoexp))
checkresiduals(fit_autoexp)

# Make forecasts for the test data
forecasts_autoexp <- forecast(fit_autoexp, h = length(test))
# Print the forecasts
print(summary(forecasts_autoexp))
# Check the residuals
checkresiduals(forecasts_autoexp)
# Evaluate the forecasts against the actual data for 2004
round(accuracy(forecasts_autoexp, test), 3)

# Extract the residual standard error (RSE) value
rse_autoexp <- sqrt(mean(residuals(fit_autoexp)^2))
rse_autoexp # Residual SD: 1.211239

# Plot the forecasts
par(mfrow =c(1,1))
plot(forecasts_autoexp, xlab = "Year", ylab = "", main = "Automated Exponential Model Forecast")
lines(test, col = "red")
legend("topleft", legend = c("Forecast", "Actual","Test"), col = c("black", "blue","red"), lty = 1)

#########################################
# Fitting an Sarima Model
########################################
# Fitting SARIMA models with different seasonal orders
# (P,D,Q)=(1,1,1), (p,d,q)=(1,0,0)
fit_sarima1 <- Arima(train, order=c(1,0,0), seasonal=list(order=c(1,1,1), period=12))
# Print the model summaries
print(summary(fit_sarima1))
# Check the residuals
checkresiduals(fit_sarima1)
# Plot the residuals
tsdisplay(residuals(fit_sarima1), lag.max=50, main="SARIMA(1,0,0)(1,1,1)[12] residual")
#Normal probability plot of the residuals.
par(mfrow =c(1,2))
qqnorm(fit_sarima1$residuals)
qqline(fit_sarima1$residuals)

# (P,D,Q)=(1,1,1), (p,d,q)=(2,0,0)
fit_sarima2 <- Arima(train, order=c(2,0,0), seasonal=list(order=c(1,1,1), period=12))
# Print the model summaries
print(summary(fit_sarima2))
# Check the residuals
checkresiduals(fit_sarima2)
# Plot the residuals
tsdisplay(residuals(fit_sarima2), lag.max=50, main="SARIMA(2,0,0)(1,1,1)[12] residual")
#Normal probability plot of the residuals.
par(mfrow =c(1,2))
qqnorm(fit_sarima2$residuals)
qqline(fit_sarima2$residuals)

# (P,D,Q)=(1,1,1), (p,d,q)=(2,0,2)
fit_sarima3 <- Arima(train, order=c(2,0,2), seasonal=list(order=c(1,1,1), period=12))
# Print the model summaries
print(summary(fit_sarima3))
# Check the residuals
checkresiduals(fit_sarima3)
# Plot the residuals
tsdisplay(residuals(fit_sarima3), lag.max=50, main="SARIMA(2,0,2)(1,1,1)[12] residual")
#Normal probability plot of the residuals.
par(mfrow =c(1,2))
qqnorm(fit_sarima3$residuals)
qqline(fit_sarima3$residuals)


# Fit an ARIMA model automatically
# (P,D,Q)=(2,1,0), (p,d,q)=(1,0,1)[12]
fit_autoarima_sarima3 <- auto.arima(train, seasonal=TRUE, stepwise=FALSE, approximation=FALSE)
# Print the model summary
print(summary(fit_autoarima_sarima3))
# Check the residuals
checkresiduals(fit_autoarima_sarima3)
# Plot the residuals
tsdisplay(residuals(fit_autoarima_sarima3), lag.max=50, main="Auto-ARIMA residual")
# Normal probability plot of the residuals.
par(mfrow =c(1,2))
qqnorm(fit_autoarima_sarima3$residuals)
qqline(fit_autoarima_sarima3$residuals)

######################################
# Best Model SARIMA(2,0,2)(1,1,1)[12]
######################################
# Evaluate the forecasts against the actual data for 2004
forecasts_sarima <- forecast(fit_sarima3, h=length(test))
# Print the forecasts
print(summary(forecasts_sarima))
# Check the residuals
checkresiduals(forecasts_sarima)
# Plot the residuals
tsdisplay(residuals(forecasts_sarima), lag.max=50, main="SARIMA(2,0,2)(1,1,1)[12] residual")
# Evaluate the forecasts against the actual data for 2004
accuracy <- round(accuracy(forecasts_sarima, test), 3)
print(accuracy)

# Plot the forecasts
par(mfrow =c(1,1))
plot(forecasts_sarima, xlab="Year", ylab="", main="SARIMA Forecast")
# Add the actual test data to the plot
lines(test, col="red")
# Add a legend
legend("topleft", legend = c("Forecasts", "Actual", "Test"), col = c("black", "blue", "red"), lty = 1)


# Extract the residual standard error (RSE) value
sarima_exp <- sqrt(mean(residuals(fit_sarima3)^2))
sarima_exp # Residual SD: 1.189311

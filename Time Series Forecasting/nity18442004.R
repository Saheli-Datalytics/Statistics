#####################################
#installing packages
#####################################
install.packages("ggplot2")
install.packages("dplyr")
install.packages("forecast")
install.packages("tseries")
install.packages("ggseas")
install.packages("fpp2")
install.packages("ggfortify")

library(ggplot2)
library(ggseas)
library(dplyr)
library(forecast)
library(tseries)
library(fpp2)
library(ggfortify)

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
yearly_data <- read.csv("nity18442004.csv")
yearly_data
class(yearly_data)

# Convert the data to a time series object
ts_data <- ts(yearly_data, start = 1844, end = 2004, frequency = 1)
head(ts_data)
class(ts_data)

start(ts_data)
end(ts_data)
frequency(ts_data)

#####################################
#Priliminary Analysis of Time Series
#####################################

# Yearly data summary statistics
yearly_summary <- summary(yearly_data$x)
yearly_summary

# Create a time plot of the data
par(mfrow =c(1,1))
plot(ts_data, main = "Time plot of yearly average temperature in Armagh")

# Difference Needed 
ndiffs(ts_data)

# Augmented Dickey-Fuller Test to check the data is stationary or not
adf.test(ts_data, alternative = "stationary")

# Differencing the data to make it stationary
diff_ts_data <- diff(ts_data)
diff_ts_data
par(mfrow =c(1,1))
plot(diff_ts_data, main = "Time plot of differenced yearly average temperature in Armagh")

# Augmented Dickey-Fuller Test to check the data is stationary or not
adf.test(diff_ts_data, alternative = "stationary")

# Visualizing based on century
plot.ts(diff_ts_data, xlim=c(1844,1900), main = "Time plot of yearly average temperature from 1844 to 1900")
plot.ts(diff_ts_data, xlim=c(1901,2004), main = "Time plot of yearly average temperature from 1901 to 2004")

# Returning the seasonal cycle of the time series ts_data
cycle(diff_ts_data)

# Compute mean and standard deviation of seasonal cycle
seasonal_mean <- mean(cycle(diff_ts_data))
seasonal_sd <- sd(cycle(diff_ts_data))

# Print out results
cat("Mean seasonal cycle:", seasonal_mean, "\n")
cat("Standard deviation of seasonal cycle:", seasonal_sd, "\n")

# A box plot can show the distribution of the values in the diff_ts_data time series and outliers
boxplot(diff_ts_data, main = "Box plot of differenced yearly average temperature in Armagh")

# Create a histogram of the data
par(mfrow =c(1,1))
hist(diff_ts_data, main = "Histogram of yearly average temperature in Armagh",
     xlab = "Temperature", ylab = "Frequency")

# Plotting Correlogram (ACF, PACF)
par(mfrow =c(1,2))
acf(diff_ts_data, lag.max = 25)
pacf(diff_ts_data, lag.max = 25)

acf(diff_ts_data, lag.max = 25, plot = FALSE)
pacf(diff_ts_data, lag.max = 25, plot = FALSE)

#######################################
# Hold Out
#######################################
# Split the data into training and test sets
train <- window(diff_ts_data, end = c(2003))
test <- window(diff_ts_data, start = c(2004))

#######################################
# Naive for forcasting
# y_t = y_(t-s) +e_t
#######################################
fit_naive <- naive(train) # Residual sd: 1.07 
print(summary(fit_naive))
checkresiduals(fit_naive)

# Make forecasts for the test data
forecasts_naive <- forecast(fit_naive, h = length(test))
# Print the forecasts
print(summary(forecasts_naive))
checkresiduals(forecasts_naive)
# Evaluate the forecasts against the actual data for 2004
round(accuracy(forecasts_naive, test),3)

# Plot the forecasts
par(mfrow =c(1,1))
plot(forecasts_naive, main = "Naive Forecasts")
# Add the actual data for comparison
lines(test, col = "red")
# Add a legend
legend("topleft", legend = c("Forecasts", "Actual", "Test"), col = c("black", "blue", "red"), lty = c(1,1,1))

#######################################
# Exponential Smoothning Model
#######################################
fit_exp <- ets(train, model = "ZZZ") # Residual sd: 0.61
print(summary(fit_exp))
checkresiduals(fit_exp)

# Make forecasts for the test data
forecasts_exp <- forecast(fit_exp, h = length(test))
print(summary(forecasts_exp))
checkresiduals(forecasts_exp)
# Evaluate the forecasts against the actual data for 2004
round(accuracy(forecasts_exp, test),3)

# Extract the residual standard error (RSE) value
rse_exp <- sqrt(mean(residuals(fit_exp)^2))
rse_exp # Residual SD:  0.61

# Plot the forecasts
par(mfrow =c(1,1))
plot(forecasts_exp, xlab = "Year", ylab = "", main = "Exponential Smoothning Forecast")
# Add the actual test data to the plot
lines(test, col = "red")
# Add a legend
legend("topleft", legend = c("Forecasts", "Actual", "Test"), col = c("black", "blue", "red"), lty = 1)

#########################################
# Fitting an arima Model
########################################
# Fitting ARIMA models with different orders
# (p,d,q)=(1,1,1)
fit_arima1 <- arima(train, order = c(1,1,1))
# Print the model summaries
print(summary(fit_arima1))
# Check the residuals
checkresiduals(fit_arima1)
# Plot the residuals
tsdisplay(residuals(fit_arima1), lag.max=50, main="ARIMA(1,1,1) residual")
#Normal probability plot of the residuals.
par(mfrow =c(1,2))
qqnorm(fit_arima1$residuals)
qqline(fit_arima1$residuals)

# (p,d,q)=(1,1,2)
fit_arima2 <- arima(train, order = c(1,1,2))
# Print the model summaries
print(summary(fit_arima2))
# Check the residuals
checkresiduals(fit_arima2)
# Plot the residuals
tsdisplay(residuals(fit_arima2), lag.max=50, main="ARIMA(1,1,2) residual")
#Normal probability plot of the residuals.
par(mfrow =c(1,2))
qqnorm(fit_arima2$residuals)
qqline(fit_arima2$residuals)

# Auto-ARIMA model selection
fit_autoarima <- auto.arima(train)
# Print the model summaries
print(summary(fit_autoarima))
# Check the residuals
checkresiduals(fit_autoarima)
# Plot the residuals
tsdisplay(residuals(fit_autoarima), lag.max=50, main="Auto-ARIMA residual")
#Normal probability plot of the residuals.
par(mfrow =c(1,2))
qqnorm(fit_autoarima$residuals)
qqline(fit_autoarima$residuals)
# Evaluating Model fit
Box.test(fit_autoarima$residuals, type = "Ljung-Box")

###########################
# Best Model ARIMA(0,0,1)
###########################
# Evaluate the forecasts against the actual data for 2004
forecasts_arima <- forecast(fit_autoarima, h=length(test))
# Print the forecasts
print(summary(forecasts_arima))
checkresiduals(forecasts_arima)
# Plot the forecasts
par(mfrow =c(1,1))
plot(forecasts_arima, xlab="Year", ylab="", main="ARIMA Forecast")
# Add the actual test data to the plot
lines(test, col="red")
# Add a legend
legend("topleft", legend = c("Forecasts", "Actual", "Test"), col = c("black", "blue", "red"), lty = 1)

# Evaluate the forecasts against the actual data for 2004
accuracy <- round(accuracy(forecasts_arima, test), 3)
print(accuracy)

# Extract the residual standard error (RSE) value
arima_exp <- sqrt(mean(residuals(fit_autoarima)^2))
arima_exp # Residual SD: 0.46

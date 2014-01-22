# To execute the simulation, type following commands in R commander:
#	source('simSpurious.R')
#	execSim(10, 1)			# For the case of 10-min integ, Poisson distribution
#
#
allanvar1 <- function(x){
	temp <- diff(x, 1, 2)
	return(0.5* (temp %*% temp) / length(temp))
}

simSpur <- function(ntime, spurAmp, randomProcess){
	if(randomProcess == 1){
		nscale <- 16; spur <- cumsum(rpois(ntime, nscale)/nscale - 1)* spurAmp* 8.0		# Random-walk Poisson spurious
	} else {
		spur <- cumsum(rnorm(ntime, mean=0, sd=1.0))* spurAmp* 2.0		# Random-walk Gaussian spurious
	}
	
	spec <- rnorm(ntime, mean=1.0, sd=sqrt(2.0e-6))
	spch <- rnorm(ntime, mean=1.0, sd=sqrt(2.0e-6)) + spur
	
	scanpattern <- rep(c(1,-1), (ntime/2))	# Scan patrtern
	spchSpec <- (spch %*% scanpattern) / (ntime/2)
	freeSpec <- (spec %*% scanpattern) / (ntime/2)
	spurSpec <- (spur %*% scanpattern) / (ntime/2)
	
	return(list(freeAV=allanvar1(spec), spchAV=allanvar1(spch), freeSpec=freeSpec, spchSpec=spchSpec, spurSpec=spurSpec))
}

execSim <- function(integMin, randomProcess){
	# integMin	: Integration time in [min]
	# randomProcess : Switch to select Random Process, 1 -> Poisson, 2 -> Gauss
	randFn <- c("Poisson", "Gauss")
	ntime <- 2*integMin*60	# Number of 0.5-sec integration periods
	spurPower <- 1.0e-4 * 2^(0:8)
	num_trial <- 256
	spurAV <- spchAV <- freeAV <- maxSPAV <- numeric(0)
	#-------- Loop for various spurious power
	for(power_index in 1:length(spurPower)){
		freeAV_tmp <- spchAV_tmp <- spurAV_tmp <- numeric(num_trial)
		cat(sprintf("Trying Power = %e\n", spurPower[power_index]))
		#-------- Loop for many trials
		for(trial_index in 1:num_trial){
			result <- simSpur(ntime, spurPower[power_index], randomProcess)
			freeAV_tmp[trial_index] <- result$freeAV
			spchAV_tmp[trial_index] <- result$spchAV
			spurAV_tmp[trial_index] <- result$spchAV - result$freeAV
		}
		spurAV[power_index] <- mean(spurAV_tmp)
		spchAV[power_index] <- mean(spchAV_tmp)
		freeAV[power_index] <- mean(freeAV_tmp)
		maxSPAV[power_index] <- sd(spurAV_tmp) + spurAV[power_index]
	}
	#-------- Plot resuts
	pdf(sprintf("SpurSim%dmin%s.pdf", integMin, randFn[randomProcess]))
	plot(spurPower, sqrt(spchAV/3), log='xy', ylim=c(5e-5, 5e-2), pch=20, col='black', xlab='Input Spurious Power [scaled by Tsys]', ylab='Detected Spurious Power [scaled by Tsys]', main= sprintf("%d-min integ", integMin)); lines(spurPower, sqrt(spchAV/3), col='black')
	points(spurPower, sqrt(freeAV/3), pch=20, col='darkgreen'); lines(spurPower, sqrt(freeAV/3), col='darkgreen')
	polygon( c(spurPower, spurPower[length(spurPower):1]), c(sqrt(maxSPAV/3), sqrt(spurAV[length(spurPower):1]/3)), col='lightskyblue', border=NA)
	points(spurPower,  sqrt(spurAV/3), pch=20, col='blue'); lines(spurPower, sqrt(spurAV/3), col='blue')
	abline(h=4.3e-4, col='red'); text(7e-2, 4.8e-4, "Requirement", cex=0.5)
	labels <- expression(paste(sigma, spth), paste(sigma, th), paste(sigma, sp))
	legend("topleft", legend=labels, col=c('black', 'darkgreen', 'blue'), lty=1)
	dev.off()
}

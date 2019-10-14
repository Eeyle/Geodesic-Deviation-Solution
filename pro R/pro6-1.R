
library(deSolve);


# Variable definitions

M <- 1; # mass
R <- 320.6; # initial radius
e <- 1/11; # eccentricity of the orbit
p <- 2400/11; # parameter to do with size of derivatives
s <- 0.1; # parameter to do with size of orbit

Art <- 2*M/(R^2)*((1-2*M/R)^(-1))*((1-3*M/R)^(-1/2)); # intermediate coefficients
Atr <- 2*M/(R^2)*(1-2*M/R)*((1-3*M/R)^(-1/2));
Arr <- (-3*M/(R^3)*((1-3*M/R)^(-1)))*(1-2*M/R);
Apr <- -2*R*((M/(R^3))^(1/2))*(1-2*M/R)*((1-3*M/R)^(-1/2));
Arp <-  2/R*((M/(R^3))^(1/2))*((1-3*M/R)^(-1/2));

wt <- (1 - 3*M/R)^(-1/2); # time frequency
w <- ((M/(R^3))^(1/2))*((1 - 3*M/R)^(-1/2))*((1 - 6*M/R)^(1/2)); # radial frequency
wp <- ((M/(R^3))^(1/2))*((1 - 3*M/R)^(-1/2)); # angle phi frequency
nr0 <- -206; # cos coefficient term for r
nt0 <- Art/w*nr0; # sin coefficient term for t 
np0 <- Arp/w*nr0; # sin coefficient term for phi
Nr0 <- 1000; # constant term for r
Nt0 <- -3/2*M/(R^2)*((1 - 3*M/R)^(-3/2))*Nr0; # constant term for t
Np0 <- -3/2*((M/R)^(1/2))/(R^2)*(1 - 2*M/R)*((1 - 3*M/R)^(-3/2))*Nr0;

# weird time factors
tfSol = 5500
tfNum = 1


# Our deviation solution. Found by analytically substituting a generic second-order Taylor series into 
# each of the equations of motion predicted by the geodesic equation in the Schwarzschild metric.

t <- seq(from = 0, to = 10, by = 0.0001);

y0 <- wt*(tfSol*t) - s*(nt0*sin(w*(tfSol*t)) +  Nt0*(tfSol*t)); # values of time as a function of proper time. 10000 is a strange scale factor 
y1 <- R + s*(nr0*cos(w*(tfSol*t)) - Nr0); # values of r
y2 <- pi/2; # value of theta is kept constant
y3 <- wp*(tfSol*t) - s*(np0*sin(w*(tfSol*t)) + Np0*(tfSol*t)); # values of phi

ycart1 <- y1 * cos(y3); # convert r and phi to 2d cartesian coordinates for plotting
ycart2 <- y1 * sin(y3);


# Numerical solution. Numerically solve the equations of motion.

dqdchi <- function(t, y, parms) {

	e <- 0.2; # global variables dont work
	p <- 12;
	
	#e <- 1/3;
	#p <- 266.66;
	a <- M # lord knows why I used a

	with(as.list(y), {
		#dR <- ( a * p^2 * (p-2-2*e)^(0.5) * (p-2+2*e)^(0.5) ) / ( (p-2-2*e*cos(tfSol*t)) * (1+e*cos(tfSol*t))^2 * (p-6-2*e*cos(tfSol*t))^(0.5) );
		dR <- ( a * (p^2) * (p-2-2*e)^(0.5) * (p-2+2*e)^(0.5) ) / ( (p-2-2*e*cos(tfSol*t)*(1+e*cos(tfSol*t)))^2 * (p-6-2*e*cos(tfSol*t))^(0.5) );
		dP <- ( p^(1/2) ) / ( (p-6-2*e*cos(tfSol*t))^(1/2) );
		list(c(dR, dP))
	})
}

#t1 <- seq(from = 0, to = 12.5, by = (1/2)^13);
yinit <- c(y1[1], y3[1]);
out <- ode(y = yinit, 
		times = tfNum*t, # lord knows why we need 12.5 
		func = dqdchi, 
		parm = NULL, 
		method = c("ode45"),
		rtol = 0.0001);

# Get the four variables from the ode solver output
q0 <- out[,1];
q1 <- out[,2];
q2 <- p * M * 1/(1 + e*cos(q0)); 
q3 <- out[,3];

qcart1 <- q2 * cos(q3); # and convert r and phi to 2d cartesian coordinates for plotting
qcart2 <- q2 * sin(q3);



# Calculate the error

# Combine the time q0 and y0 into one long vector, since they may be different, and sort it. 
combinedTime <- sort(c(q0/tfNum, y0/tfSol));

# Then interpolate q1, y1, q3, and y3 over those particular times.
q1inter <- approx(t, q2, xout=combinedTime);
q3inter <- approx(t, q3, xout=combinedTime);
y1inter <- approx(t, y1, xout=combinedTime);
y3inter <- approx(t, y3, xout=combinedTime);

errorr <- (y1inter$y - q1inter$y) / q1inter$y;
errorp <- (y3inter$y - q3inter$y) / q3inter$y;

#q1inter <- smooth.spline(q0, q1, n=length(y1), method="fmm");
#q3inter <- smooth.spline(q0, q3, n=length(y1), method="fmm");

#errorr <- (y1 - q1inter$y) / q1inter$y;
#errorp <- (y3 - q3inter$y) / q3inter$y;


# Plot

par(mfrow=c(1, 2));

plot(qcart1, qcart2, type="l", col="black", main="Path around black hole", xlab="", ylab="");
lines(ycart1, ycart2, col="green");
legend(-240, 240, legend=c("Numeric", "Analytic"), col=c("black", "green"), lty=1);

plot(q1inter$x, errorr, type="l", col="black", main="Relative error", ylim=c(-0.2, 0.4), xlim=c(0, 5),
	ylab="relative error (%)", xlab="time (s)", axes=FALSE);
lines(q3inter$x, errorp, type="l", col="blue");
axis(1, at=0:5, pos=0);
axis(2, at=seq(-1, 1, 0.1)); 
legend(3, 0.4, legend=c("Radius", "Angle"), col=c("black", "blue"), lty=1);






library(deSolve);

paramsV <- c(c(M=1, R=395, e=1/3, p=266.66, s=0.1, nr0=-950, tfSol=10000, tfNum=1.25, e2=0.2, p2=12, lx=-360, ly=400, lex=4, ley=0.4, yf=1), 
		c(M=1, R=320.6, e=1/11, p=2400/11, s=0.1, nr0=-206, tfSol=5500, tfNum=1, e2=0.2, p2=12, lx=-240, ly=240, lex=3, ley=0.4, yf=1), 
		c(M=1, R=1119.5, e=2/100, p=1020, s=0.1, nr0=-195, tfSol=20000, tfNum=1.25, e2=0.02, p2=1020, lx=-1000, ly=1000, lex=3, ley=0.4, yf=2));
params <- array(paramsV, c(11, 3), list(c("M", "R", "e", "p", "s", "nr0" "tfSol", "tfNum", "e2", "p2", "lx", "ly", "lex", "ley", "yf")));

par(mfrow=c(3, 2));

for (i in 1:3) {
	


# Variable definitions

M <- params["M", i]; # mass
R <- params["R", i]; # initial radius
e <- params["e", i]; # eccentricity of the orbit
p <- params["p", i]; # parameter to do with size of derivatives
s <- params["s", i]; # parameter to do with size of orbit

Art <- 2*M/(R^2)*((1-2*M/R)^(-1))*((1-3*M/R)^(-1/2)); # intermediate coefficients
Atr <- 2*M/(R^2)*(1-2*M/R)*((1-3*M/R)^(-1/2));
Arr <- (-3*M/(R^3)*((1-3*M/R)^(-1)))*(1-2*M/R);
Apr <- -2*R*((M/(R^3))^(1/2))*(1-2*M/R)*((1-3*M/R)^(-1/2));
Arp <-  2/R*((M/(R^3))^(1/2))*((1-3*M/R)^(-1/2));

wt <- (1 - 3*M/R)^(-1/2); # time frequency
w <- ((M/(R^3))^(1/2))*((1 - 3*M/R)^(-1/2))*((1 - 6*M/R)^(1/2)); # radial frequency
wp <- ((M/(R^3))^(1/2))*((1 - 3*M/R)^(-1/2)); # angle phi frequency
nr0 <- params["nr0", i]; # cos coefficient term for r
nt0 <- Art/w*nr0; # sin coefficient term for t 
np0 <- Arp/w*nr0; # sin coefficient term for phi
Nr0 <- 1000; # constant term for r
Nt0 <- -3/2*M/(R^2)*((1 - 3*M/R)^(-3/2))*Nr0; # constant term for t
Np0 <- -3/2*((M/R)^(1/2))/(R^2)*(1 - 2*M/R)*((1 - 3*M/R)^(-3/2))*Nr0;

# weird time factors
tfSol = params["tfSol", i]
tfNum = params["tfNum", i]


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

	e <- params["e2", i]; # global variables dont work
	p <- params["p2", i];
	
	#e <- 1/3;
	#p <- 266.66;
	a <- M # lord knows why I used a

	with(as.list(y), {
		dR <- ( a * p^2 * (p-2-2*e)^(0.5) * (p-2+2*e)^(0.5) ) / ( (p-2-2*e*cos(tfSol*t)) * (1+e*cos(tfSol*t))^2 * (p-6-2*e*cos(tfSol*t))^(0.5) );
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
		rtol = 0.001);

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
y3inter <- approx(t, params["yf", i]*y3, xout=combinedTime);

errorr <- (y1inter$y - q1inter$y) / q1inter$y;
errorp <- (y3inter$y - q3inter$y) / q3inter$y;

#q1inter <- smooth.spline(q0, q1, n=length(y1), method="fmm");
#q3inter <- smooth.spline(q0, q3, n=length(y1), method="fmm");

#errorr <- (y1 - q1inter$y) / q1inter$y;
#errorp <- (y3 - q3inter$y) / q3inter$y;


# Plot

plot(qcart1, qcart2, type="l", col="black", main="Path around black hole", xlab="", ylab="");
lines(ycart1, ycart2, col="green");
legend(params["lx", i], params["ly", i], legend=c("Numeric", "Analytic"), col=c("black", "green"), lty=1);

plot(q1inter$x, errorr, type="l", col="black", main="Relative error", ylim=c(-0.2, 0.4), xlim=c(0, 5),
	ylab="relative error (%)", xlab="time (s)", axes=FALSE);
lines(q3inter$x, errorp, type="l", col="blue");
axis(1, at=2*0:6, pos=0);
axis(2, at=seq(-1, 1, 0.1)); 
legend(params["lex", i], params["ley", i], legend=c("Radius", "Angle"), col=c("black", "blue"), lty=1);
	
}


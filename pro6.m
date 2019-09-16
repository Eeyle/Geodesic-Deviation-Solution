%% ---- R200e1/3
% Run these one section at a time with the Run Section button. Ensure
% pro3.m and pro5.m are included in the current folder and are on the path.

% -- Variables

M = 1;
R = 395;
e = 1/3;
p = 266.0+2/3;
s = 0.1;

Art = 2*M/(R^2)*((1-2*M/R)^(-1))*((1-3*M/R)^(-1/2));
Atr = 2*M/(R^2)*(1-2*M/R)*((1-3*M/R)^(-1/2));
Arr = (-3*M/(R^3)*((1-3*M/R)^(-1)))*(1-2*M/R);
Apr = -2*R*((M/(R^3))^(1/2))*(1-2*M/R)*((1-3*M/R)^(-1/2));
Arp =  2/R*((M/(R^3))^(1/2))*((1-3*M/R)^(-1/2));

wt = (1 - 3*M/R)^(-1/2);
w = ((M/(R^3))^(1/2))*((1 - 3*M/R)^(-1/2))*((1 - 6*M/R)^(1/2));
wp = ((M/(R^3))^(1/2))*((1 - 3*M/R)^(-1/2));
nr0 = -950;
nt0 = Art/w*nr0;
np0 = Arp/w*nr0;
Nr0 = 1000;
Nt0 = -3/2*M/(R^2)*((1 - 3*M/R)^(-3/2))*Nr0;
Np0 = -3/2*((M/R)^(1/2))/(R^2)*(1 - 2*M/R)*((1 - 3*M/R)^(-3/2))*Nr0;


% -- Our solution

t = 0:10:100000; 

y = zeros(length(t), 4);
y(:,1) = wt*t - s*(nt0*sin(w*t) +  Nt0*t);
y(:,2) = R + s*(nr0*cos(w*t) - Nr0);
y(:,3) = pi/2;
y(:,4) = wp*t - s*(np0*sin(w*t) + Np0*t);

ycart = zeros(length(t),2);
ycart(:,1) = y(:,2).*cos(y(:,4)); % convert from polar to cartesian
ycart(:,2) = y(:,2).*sin(y(:,4)); 


% -- Numerical solution

options = odeset('Refine',512,'RelTol',0.0001);
[X1,q1] = ode45(@(X1,q1,e,p,a)pro3(X1,q1,0.2,12,M), [0 12.5], [0], options);
[X3,q3] = ode45(@(X3,q3,e,p,a)pro5(X3,q3,0.2,12,M), [0 12.5], [0], options);
q = zeros(length(q1), 3);
q(:,1) = q1;
%q(:,2) = (1./(1+e*cos(X1)))*(p);
zabadabadoo = 1./(1+e*cos(X1));
zabadabadee = p*M*zabadabadoo;
q(:,2) = zabadabadee;
q(:,3) = spline(X3,q3,X1);

qcart = zeros(length(q(:,1)),2);
qcart(:,1) = q(:,2).*cos(q(:,3)); % convert from polar to cartesian
qcart(:,2) = q(:,2).*sin(q(:,3));
  

% -- Plotting

viscircles([0 0], 3, 'Color', 'y'); hold on; % Event horizon * 1.5
plot(0,0,'k.'); hold on; % Black hole
axis([-400 400 -400 400]);  

l1 = animatedline('Color','b');       % Animated lines to draw the paths
l2 = animatedline('Color','g');       % Animated lines to draw the paths

for k = 1:max(length(qcart(:,1)), length(ycart(:,1)))

    if k <= length(qcart(:,1))
        addpoints(l1, qcart(k,1), qcart(k,2));
    end
    if k <= length(ycart(:,1))
       addpoints(l2, ycart(k,1), ycart(k,2));
    end

    pause(0.00005);
    drawnow;
end




%% ---- R200e1/11 with error


% -- Variables

M = 1;
R = 320.6;
e = 1/11;
p = 2400/11;
s = 0.1;

Art = 2*M/(R^2)*((1-2*M/R)^(-1))*((1-3*M/R)^(-1/2));
Atr = 2*M/(R^2)*(1-2*M/R)*((1-3*M/R)^(-1/2));
Arr = (-3*M/(R^3)*((1-3*M/R)^(-1)))*(1-2*M/R);
Apr = -2*R*((M/(R^3))^(1/2))*(1-2*M/R)*((1-3*M/R)^(-1/2));
Arp =  2/R*((M/(R^3))^(1/2))*((1-3*M/R)^(-1/2));

wt = (1 - 3*M/R)^(-1/2);
w = ((M/(R^3))^(1/2))*((1 - 3*M/R)^(-1/2))*((1 - 6*M/R)^(1/2));
wp = ((M/(R^3))^(1/2))*((1 - 3*M/R)^(-1/2));
nr0 = -206;
nt0 = Art/w*nr0;
np0 = Arp/w*nr0;
Nr0 = 1000;
Nt0 = -3/2*M/(R^2)*((1 - 3*M/R)^(-3/2))*Nr0;
Np0 = -3/2*((M/R)^(1/2))/(R^2)*(1 - 2*M/R)*((1 - 3*M/R)^(-3/2))*Nr0;


% -- Our solution

t = 0:8:55000; 

y = zeros(length(t), 4);
y(:,1) = wt*t - s*(nt0*sin(w*t) +  Nt0*t);
y(:,2) = R + s*(nr0*cos(w*t) - Nr0);
y(:,3) = pi/2;
y(:,4) = wp*t - s*(np0*sin(w*t) + Np0*t);

ycart = zeros(length(t),2);
ycart(:,1) = y(:,2).*cos(y(:,4)); % convert from polar to cartesian
ycart(:,2) = y(:,2).*sin(y(:,4)); 


% -- Numerical solution

options = odeset('Refine',512,'RelTol',0.00001);
[X1,q1] = ode45(@(X1,q1,e,p,a)pro3(X1,q1,0.2,12,M), [0 10], [0], options);
[X3,q3] = ode45(@(X3,q3,e,p,a)pro5(X3,q3,0.2,12,M), [0 10], [0], options);
q = zeros(length(q1), 3);
q(:,1) = q1;
%q(:,2) = (1./(1+e*cos(X1)))*(p*M);
zabadabadoo = 1./(1+e*cos(X1));
zabadabadee = p*M*zabadabadoo;
q(:,2) = zabadabadee;
q(:,3) = spline(X3,q3,X1);

qcart = zeros(length(q(:,1)),2);
qcart(:,1) = q(:,2).*cos(q(:,3)); % convert from polar to cartesian
qcart(:,2) = q(:,2).*sin(q(:,3));
  

% -- Finding the error

qatt = zeros(length(t),3);
iHateMatlab = spline(X1, q(:,2), t/(200000/10)); 
qatt(:,2) = iHateMatlab;
qatt(:,3) = spline(X3, q3, t/(200000/10));

errorr = (y(:,2)-qatt(:,2))./qatt(:,2);
errorp = (y(:,4)-qatt(:,3))./qatt(:,3);

% -- Plotting

figure('units','normalized','outerposition',[0 0 1 1]);
%set(gca,'Position',[0.1 0.1 0.75 0.85]);

leftPlot = subplot(1,2,1);
viscircles([0 0], 3, 'Color', 'y'); hold on; % Event horizon * 1.5
plot(0,0,'k.'); hold on; % Black hole
axis([-250 250 -250 250]);  
title("The Orbit");

rightPlot = subplot(1,2,2);
plot(1:length(t), zeros(length(t),1), 'k'); hold on; % black line
axis([0 length(t) -0.2 0.4]);
title("Error in R");
xticklabels({});
xlabel("Time");
ylabel("Relative Error");

l1 = animatedline('Color','b','Parent',leftPlot);       % Animated lines to draw the paths
l2 = animatedline('Color','r','Parent',leftPlot);
l3 = animatedline('Color','k','Parent',rightPlot);

for k = 1:max([length(qcart(:,1)), length(ycart(:,1)), length(errorr)])
    subplot(1,2,1);
    if k+1000 <= length(qcart(:,1))
       addpoints(l1, qcart(k+1000,1), qcart(k+1000,2));
    end
    if k <= length(ycart(:,1))
       addpoints(l2, ycart(k,1), ycart(k,2));
    end
    
    subplot(1,2,2);
    if (k-500 <= length(errorr) && k > 500)
        addpoints(l3, k-500, errorr(k-500))
    end
    
    %pause(0.0001);
    drawnow; hold on;
end

%{
 for l = k:(k+84)
    subplot(1,2,2);
    addpoints(l3, l-500, errorr(l-500));
    drawnow; hold on;
 end
%}



%% ---- R1000e1/100


% -- Variables

M = 1;
R = 1119.5;
e = 2/100;
p = 1020;
s = 0.1;

Art = 2*M/(R^2)*((1-2*M/R)^(-1))*((1-3*M/R)^(-1/2));
Atr = 2*M/(R^2)*(1-2*M/R)*((1-3*M/R)^(-1/2));
Arr = (-3*M/(R^3)*((1-3*M/R)^(-1)))*(1-2*M/R);
Apr = -2*R*((M/(R^3))^(1/2))*(1-2*M/R)*((1-3*M/R)^(-1/2));
Arp =  2/R*((M/(R^3))^(1/2))*((1-3*M/R)^(-1/2));

wt = (1 - 3*M/R)^(-1/2);
w = ((M/(R^3))^(1/2))*((1 - 3*M/R)^(-1/2))*((1 - 6*M/R)^(1/2));
wp = ((M/(R^3))^(1/2))*((1 - 3*M/R)^(-1/2));
nr0 = -195;
nt0 = Art/w*nr0;
np0 = Arp/w*nr0;
Nr0 = 1000;
Nt0 = -3/2*M/(R^2)*((1 - 3*M/R)^(-3/2))*Nr0;
Np0 = -3/2*((M/R)^(1/2))/(R^2)*(1 - 2*M/R)*((1 - 3*M/R)^(-3/2))*Nr0;


% -- Our solution

t = 0:40:200000; 

y = zeros(length(t), 4);
y(:,1) = wt*t - s*(nt0*sin(w*t) +  Nt0*t);
y(:,2) = R + s*(nr0*cos(w*t) - Nr0);
y(:,3) = pi/2;
y(:,4) = wp*t - s*(np0*sin(w*t) + Np0*t);

ycart = zeros(length(t),2);
ycart(:,1) = y(:,2).*cos(y(:,4)); % convert from polar to cartesian
ycart(:,2) = y(:,2).*sin(y(:,4)); 


% -- Numerical solution

options = odeset('Refine',128,'RelTol',0.0001);
[X1,q1] = ode45(@(X1,q1,e,p,a)pro3(X1,q1,0.02,1020,M), [0 10], [0], options);
[X3,q3] = ode45(@(X3,q3,e,p,a)pro5(X3,q3,0.02,1020,M), [0 10], [0], options);
q = zeros(length(q1), 3);
q(:,1) = q1;
%q(:,2) = (1./(1+e*cos(X1)))*(p);
zabadabadoo = 1./(1+e*cos(X1));
zabadabadee = p*M*zabadabadoo;
q(:,2) = zabadabadee;
q(:,3) = spline(X3,q3,X1);

qcart = zeros(length(q(:,1)),2);
qcart(:,1) = q(:,2).*cos(q(:,3)); % convert from polar to cartesian
qcart(:,2) = q(:,2).*sin(q(:,3));


% -- Plotting

viscircles([0 0], 3, 'Color', 'y'); hold on; % Event horizon * 1.5
plot(0,0,'k.'); hold on; % Black hole
axis([-1100 1100 -1100 1100]);  

l1 = animatedline('Color','b');       % Animated lines to draw the paths
l2 = animatedline('Color','g');       % Animated lines to draw the paths

for k = 1:max(length(qcart(:,1)), length(ycart(:,1)))
    if k <= length(qcart(:,1))
        addpoints(l1, qcart(k,1), qcart(k,2));
    end
    if k <= length(ycart(:,1))
       addpoints(l2, ycart(k,1), ycart(k,2));
    end

    pause(0.001);
    drawnow;
end


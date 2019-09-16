function dqdchi = pro3(X, q, e, p, a)

dqdchi(1) = a*p^2*(p-2-2*e)^(1/2)*(p-2+2*e)^(1/2)/((p-2-2*e*cos(X))*(1+e*cos(X))^2*(p-6-2*e*cos(X))^(1/2)); % dt/dX

end
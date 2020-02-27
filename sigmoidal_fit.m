function out=sigmoidal_fit(x,y,bl,amp,phi,lambda,x0,y0)
out = bl+amp./(1+exp(-1*lambda*((x-x0)*sin(phi)+(y-y0)*cos(phi))));
end
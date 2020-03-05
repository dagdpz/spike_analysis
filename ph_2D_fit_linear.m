function out=ph_2D_fit_linear(x,y,bl,slope,phi)
out = bl+slope.*(x*sin(phi)+y*cos(phi));
end
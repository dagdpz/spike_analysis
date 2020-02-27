function out=linear_fit(x,y,bl,slope,phi)
out = bl+slope.*(x*sin(phi)+y*cos(phi));
end
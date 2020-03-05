function out=ph_2D_fit_gaussian_1_RF(x,y,bl,phi,xmax,ymax,zmax,sx,sy,minratio)
        y_sd_12=sy + sx*((minratio-sy/sx)*(sign(minratio-sy/sx)+1)/2);
        out = bl + zmax*exp(-((x-xmax)*cos(phi*-1)-(y-ymax)*sin(phi*-1)).^2/(2*sx^2) -((x-xmax)*sin(phi*-1)+(y-ymax)*cos(phi*-1)).^2/(2*y_sd_12^2));
end
function out=ph_2D_fit_gaussian_2nd_RF(x,y,bl,phi,xmax,ymax,zmax,sx,sy,minratio,xmax1,ymax1,phi1,sx1,sy1)

vector_between_two_peaks=xmax+ymax*1i-xmax1-ymax1*1i;
y_sd_1=sy1 + sx1*((minratio-sy1/sx1)*(sign(minratio-sy1/sx1)+1)/2);
y_sd_2=sy + sx*((minratio-sy/sx)*(sign(minratio-sy/sx)+1)/2);
gaussian1_sd_in_direction=2*sqrt((sx1*cos(angle(vector_between_two_peaks)   -phi1))^2 + (y_sd_1*sin(angle(vector_between_two_peaks)   -phi1))^2);
gaussian2_sd_in_direction=2*sqrt((sx*cos(angle(vector_between_two_peaks*-1)-phi))^2 + (y_sd_2*sin(angle(vector_between_two_peaks*-1)-phi))^2);
radial_distance_to_correct=abs(vector_between_two_peaks)-gaussian1_sd_in_direction-gaussian2_sd_in_direction;
if radial_distance_to_correct>0
    radial_distance_to_correct=0;
end
y_shift=sin(angle(vector_between_two_peaks))*abs(radial_distance_to_correct);
x_shift=cos(angle(vector_between_two_peaks))*abs(radial_distance_to_correct);
xmax=xmax+x_shift;
ymax=ymax+y_shift;
out = bl + zmax*exp(-((x-xmax)*cos(phi*-1)-(y-ymax)*sin(phi*-1)).^2/(2*sx^2) -((x-xmax)*sin(phi*-1)+(y-ymax)*cos(phi*-1)).^2/(2*y_sd_2^2));
end

function out=twoDgaussian_with_2_RFs_fit(x,y,bl,phi1,xmax1,ymax1,zmax1,sx1,sy1,zmax2,dist,ratio,minratio)


% ,phi2,xmax2,ymax2,zmax2,sx2,sy2

phi2=phi1;
xmax2=xmax1+dist*cos(phi1-pi/2)/2;
ymax2=ymax1+dist*sin(phi1-pi/2)/2;
xmax1=xmax1-dist*cos(phi1-pi/2)/2;
ymax1=ymax1-dist*sin(phi1-pi/2)/2;
sx2=sx1*ratio;
sy2=sy1*ratio;
vector_between_two_peaks=xmax2+ymax2*1i-xmax1-ymax1*1i;


y_sd_1=sy1 + sx1*((minratio-sy1/sx1)*(sign(minratio-sy1/sx1)+1)/2);
y_sd_2=sy2 +sx2*((minratio-sy2/sx2)*(sign(minratio-sy2/sx2)+1)/2);
gaussian1_sd_in_direction=2*sqrt((sx1*cos(angle(vector_between_two_peaks)   -phi1))^2 + (y_sd_1*sin(angle(vector_between_two_peaks)   -phi1))^2);
gaussian2_sd_in_direction=2*sqrt((sx2*cos(angle(vector_between_two_peaks*-1)-phi2))^2 + (y_sd_2*sin(angle(vector_between_two_peaks*-1)-phi2))^2);
radial_distance_to_correct=abs(vector_between_two_peaks)-gaussian1_sd_in_direction-gaussian2_sd_in_direction;
if radial_distance_to_correct>0
    radial_distance_to_correct=0;
end
y_shift=sin(angle(vector_between_two_peaks))*abs(radial_distance_to_correct)/2;
x_shift=cos(angle(vector_between_two_peaks))*abs(radial_distance_to_correct)/2;


xmax1=xmax1-x_shift;
ymax1=ymax1-y_shift;
xmax2=xmax2+x_shift;
ymax2=ymax2+y_shift;

out = bl + zmax1*...
    exp(-((x-xmax1)*cos(phi1*-1)-(y-ymax1)*sin(phi1*-1)).^2/(2*(sx1)^2)...
    -((x-xmax1)*sin(phi1*-1)+(y-ymax1)*cos(phi1*-1)).^2/(2*(y_sd_1)^2))+...
    zmax2*sign(zmax1)*sign(zmax2)*-1*...
    exp(-((x-xmax2)*cos(phi2*-1)-(y-ymax2)*sin(phi2*-1)).^2/(2*(sx2)^2)...
    -((x-xmax2)*sin(phi2*-1)+(y-ymax2)*cos(phi2*-1)).^2/(2*(y_sd_2)^2));
end

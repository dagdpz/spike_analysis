function out=sigmoidal_fit(x,y,bl,amp,phi,lambda,x0,y0)
        %out = a(1)+a(2)./(1+exp(-1*a(4)*([x(:,1)-a(5)  x(:,2)-a(6)]*[sin(a(3));cos(a(3))])));
        out = bl+amp./(1+exp(-1*lambda*((x-x0)*sin(phi)+(y-y0)*cos(phi))));
        %out = bl+amp./(1+exp(-1*lambda*([x(:,1)-x0  x(:,2)-y0]*[sin(phi);cos(phi)])));
        
       % out = bl+amp./(1+exp(-1*lambda*([x-x0]*[sin(phi);cos(phi)]))) +  y-y0;
%         
%         
%             Out.(fitfun).bl=coe(1);
%             Out.(fitfun).amp=coe(2);
%             Out.(fitfun).phi=coe(3);
%             Out.(fitfun).lambda=coe(4);
%             Out.(fitfun).x0=coe(5);
%             Out.(fitfun).y0=coe(6);
end
function Out=fit_multiples(xin,yin,zin,fitsettings)
% if  numel(zin)<4 || numel(unique(zin))<=1
%
%     %                     [RF_tmp.Zout RF_tmp.sx1 RF_tmp.sy1 RF_tmp.phi1 RF_tmp.xmax1 RF_tmp.ymax1 RF_tmp.zmax1...
%     %                         RF_tmp.sx2 RF_tmp.sy2 RF_tmp.phi2 RF_tmp.xmax2 RF_tmp.ymax2 RF_tmp.zmax2]=deal(NaN);
%     return;
% end

xout=fitsettings.xout;
yout=fitsettings.yout;
Xout=combvec(xout,yout)';
Z=double(zin(:));

X=[xin(~isnan(Z)),yin(~isnan(Z))];
Z=Z(~isnan(Z));%-double(baseline(~isnan(Z)));

%% temporary part to fit on per position averages instead of single trials
[unique_positions, ~, pos_idx]=unique(X,'rows');
Z_mean=NaN(size(unique_positions,1),1);
for p=1:max(pos_idx)
    Z_mean(p)=mean(Z(pos_idx==p));
    Z_weights(p)=1/var(Z(pos_idx==p));
end
Z_weights(isinf(Z_weights))=mean(Z_weights(~isinf(Z_weights))); %% is this correct?? can weight be bigger than 1? it can....

X=unique_positions;
x=X(:,1);
y=X(:,2);
dx=min(diff(unique(x)))/2;
dy=min(diff(unique(y)))/2;

opts = optimset('lsqcurvefit');
opts.Display='off';

fitfuns=fitsettings.fittypes; %[FO, G, O] = FIT(X, Y, ...)
for f=1:numel(fitfuns)
    fitfun=fitfuns{f};
    %clear coe
z=Z_mean;
    switch fitfun
        case 'linear'
            %  %baseline       %slope         rotation
            LB=[min(Z)         0               -inf         ];
            X0=[0              1               0            ];
            UB=[max(Z)         inf             inf           ];
            
            fitT=fittype('linear_fit( x, y, bl, slope, phi )','dependent',{'z'},'independent',{'x','y'},'coefficients',{'bl', 'slope', 'phi'});
        case 'sigmoidal'
            %  %baseline       %amp           rotation     lambda (min so that 3/4 are covering one target)                     x0                   y0
            LB=[min(Z)         0               -inf         0                                                      min(x)+dx         min(y)+dy];
            X0=[0              max(abs(Z))     0            0.5                                                     0                     0];
            UB=[max(Z)         2*max(abs(Z))   inf          -1*log(1/4)/min([dx,dy])                                max(x)-dx         max(y)-dy];
            
            fitT=fittype('sigmoidal_fit( x, y, bl, amp, phi, lambda, x0, y0 )','dependent',{'z'},'independent',{'x','y'},'coefficients',{'bl', 'amp', 'phi', 'lambda', 'x0', 'y0'});
            
        case 'gaussian1'
            range_factor=fitsettings.range_factor;
            x_range=(max(max(xin))-min(min(xin)))*range_factor; if x_range==0; x_range=1; end;
            y_range=(max(max(yin))-min(min(yin)))*range_factor; if y_range==0; y_range=1; end;
            sd_max=fitsettings.sd_max_x;
            zscaling=max(Z)-min(Z);
            start_pos_xy_3=[mean(x) mean(y)];%[nanmean(xin(~isnan(Z)).*abs(Z))/nanmean(abs(Z)) nanmean(yin(~isnan(Z)).*abs(Z))/nanmean(abs(Z))];
            
            sd_min_ratio=fitsettings.sd_x_min_ratio;
            
            
            % baseline      rotation x0                 y0                    peak             sx                           sy       sigma x(ratio to sd_max_x)  ratio xy
            LB=[min([z;0])  -inf        min(x)+dx          min(y)+dy               -1.5*zscaling    sd_min_ratio*sd_max          sd_min_ratio*sd_max              ];
            X0=[0           0        start_pos_xy_3(1)  start_pos_xy_3(2)       zscaling         (sd_min_ratio+1)/2*sd_max    (sd_min_ratio+1)/2*sd_max        ];
            UB=[max(z)      inf       max(x)-dx          max(y)-dy               1.5*zscaling     sd_max                       sd_max                           ];
            
            fitT=fittype(['twoDgaussian_1_RF_fit( x, y, bl, phi, xmax, ymax, zmax, sx, sy, ' num2str(fitsettings.sd_xy_min_ratio) ' )'],'dependent',{'z'},'independent',{'x','y'},'coefficients',{'bl', 'phi', 'xmax', 'ymax', 'zmax', 'sx', 'sy'});
            
        case 'gaussian2'
            z=Z_mean-nanmean(Z_mean);
            range_factor=fitsettings.range_factor;
            x_range=(max(max(xin))-min(min(xin)))*range_factor; if x_range==0; x_range=1; end;
            y_range=(max(max(yin))-min(min(yin)))*range_factor; if y_range==0; y_range=1; end;
            sd_max_x=fitsettings.sd_max_x;
            sd_max_y=fitsettings.sd_max_y;
            
            zscaling=max(Z)-min(Z);%max(abs(Z));
            start_pos_xy_1=[mean(x(x>0)) mean(y)];%[18*sign(nanmean(xin(~isnan(Z)).*Z)/nanmean(Z)) 0]; % 18 here is arbitrary though
            start_pos_xy_2=start_pos_xy_1*-1;
            %start_pos_xy_3=[nanmean(xin(~isnan(Z)).*abs(Z))/nanmean(abs(Z)) nanmean(yin(~isnan(Z)).*abs(Z))/nanmean(abs(Z))];
            
            sd_min_ratio=fitsettings.sd_x_min_ratio;
            sd_max=fitsettings.sd_max_x;
            
            % baseline rotation x0                      y0                  peak                sx                           sy                                  igmax(ratio to sd_max_x)  ratio xy
            LB=[-inf     min(x)+dx             min(y)+dy           -1.5*zscaling       sd_min_ratio*sd_max          sd_min_ratio*sd_max         ...
                -inf     min(x)+dx             min(y)+dy           0                   sd_min_ratio*sd_max          sd_min_ratio*sd_max         ];
            X0=[pi/2  start_pos_xy_1(1)     start_pos_xy_1(2)   -zscaling                   (sd_min_ratio+1)/2*sd_max    (sd_min_ratio+1)/2*sd_max   ...
                pi/2  start_pos_xy_2(1)     start_pos_xy_2(2)   zscaling                   (sd_min_ratio+1)/2*sd_max    (sd_min_ratio+1)/2*sd_max   ];
            UB=[inf    max(x)-dx             max(y)-dy           0                   sd_max                       sd_max                      ...
               inf    max(x)-dx             max(y)-dy           1.5*zscaling        sd_max                       sd_max                      ];
            
            fitT=fittype(['twoDgaussian_with_2_opposing_RFs_fit( x, y, phi1, xmax1, ymax1, zmax1, sx1, sy1, phi2, xmax2, ymax2, zmax2, sx2, sy2, ' num2str(fitsettings.sd_xy_min_ratio) ' )'],'dependent',{'z'},'independent',{'x','y'},'coefficients',{'phi1', 'xmax1', 'ymax1', 'zmax1', 'sx1', 'sy1', 'phi2', 'xmax2', 'ymax2', 'zmax2', 'sx2', 'sy2'});
            
    end
    same_bounds=LB==UB;
    LB(same_bounds)=LB(same_bounds)-0.01;
    UB(same_bounds)=UB(same_bounds)+0.01;
    fitopts=fitoptions('method','NonlinearLeastSquares','Lower',LB,'Upper',UB,'Start',X0);
    if sum(~isnan(z)) > numel(X0)
        [fitobj, Goodness] =fit([x,y],z,fitT,fitopts);
    else
        [fitobj, Goodness] =fit([x,y],zeros(size(X0)),fitT,fitopts);
    end
    FN=fieldnames(fitobj);
    for fn=1:numel(FN)
        Out.(fitfun).(FN{fn})=fitobj.(FN{fn});
    end
    Out.(fitfun).R2=Goodness.rsquare;
    R2_adjusted(f)=Goodness.adjrsquare;
    Out.(fitfun).R2_adjusted=R2_adjusted(f);
    [~,~,~,fitfun_valid(f)]=stepwisefit(fitobj(x,y),z,'display','off');
    %fitfun_valid(f)=true;
    
    Out.(fitfun).fitfun_valid=fitfun_valid(f);
    Zout=fitobj(Xout(:,1),Xout(:,2));
    Out.(fitfun).Zout=reshape(Zout,numel(xout),numel(yout));
    
    switch fitfun %% so here we make sure the outputs are actually meaningful
        case 'gaussian1'
            minratio=fitsettings.sd_xy_min_ratio;
            sx=Out.(fitfun).sx;
            sy=Out.(fitfun).sy;
            sy=sy +sx*((minratio-sy/sx)*(sign(minratio-sy/sx)+1)/2);
            Out.(fitfun).sy=sy;
        case 'gaussian2'
            
            minratio=fitsettings.sd_xy_min_ratio;
            xmax1=Out.(fitfun).xmax1;
            ymax1=Out.(fitfun).ymax1;
            xmax2=Out.(fitfun).xmax2;
            ymax2=Out.(fitfun).ymax2;
            sx1=Out.(fitfun).sx1;
            sy1=Out.(fitfun).sy1;
            sx2=Out.(fitfun).sx2;
            sy2=Out.(fitfun).sy2;
            phi1=Out.(fitfun).phi1;
            phi2=Out.(fitfun).phi2;
            
            vector_between_two_peaks=xmax2+ymax2*1i-xmax1-ymax1*1i;
            y_sd_1=sy1 + sx1*((minratio-sy1/sx1)*(sign(minratio-sy1/sx1)+1)/2);
            y_sd_2=sy2 + sx2*((minratio-sy2/sx2)*(sign(minratio-sy2/sx2)+1)/2);
            gaussian1_sd_in_direction=2*sqrt((sx1*cos(angle(vector_between_two_peaks)   -phi1))^2 + (y_sd_1*sin(angle(vector_between_two_peaks)   -phi1))^2);
            gaussian2_sd_in_direction=2*sqrt((sx2*cos(angle(vector_between_two_peaks*-1)-phi2))^2 + (y_sd_2*sin(angle(vector_between_two_peaks*-1)-phi2))^2);
            radial_distance_to_correct=abs(vector_between_two_peaks)-gaussian1_sd_in_direction-gaussian2_sd_in_direction;
            if radial_distance_to_correct>0
                radial_distance_to_correct=0;
            end
            y_shift=sin(angle(vector_between_two_peaks))*abs(radial_distance_to_correct);
            x_shift=cos(angle(vector_between_two_peaks))*abs(radial_distance_to_correct);
            
            Out.(fitfun).sy1=y_sd_1;
            Out.(fitfun).sy2=y_sd_2;
            Out.(fitfun).xmax2=xmax2+x_shift;
            Out.(fitfun).ymax2=ymax2+y_shift;
            
    end
    
end
Out.bestfit='none';
Out.secondbestfit='none';
Zout=NaN(size(X));
R2_a=0;
if any(fitfun_valid)
%     best_R2_index=R2_adjusted==max(R2_adjusted(fitfun_valid));
%     secondbest_R2_index=R2_adjusted==max(R2_adjusted(~best_R2_index));
%     R2_a=max(R2_adjusted(fitfun_valid));

best_R2_index=R2_adjusted==max(R2_adjusted);
    secondbest_R2_index=R2_adjusted==max(R2_adjusted(~best_R2_index));
    R2_a=max(R2_adjusted);
    
    Out.bestfit=fitfuns{best_R2_index};
    Out.secondbestfit=fitfuns{secondbest_R2_index};
    %R2_a=Out.(Out.bestfit).R2;
    Zout=Out.(Out.bestfit).Zout;
end
Out.R2_adjusted=R2_a;
Out.Zout=Zout;
Out.none=struct;
end

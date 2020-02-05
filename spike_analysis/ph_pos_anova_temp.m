function anova_struct=ph_pos_anova_temp(keys,Pos_trial)

%FR,States,choices,Position_idx,Fixation_idx,Positions,Fixation,hands

States                              ={Pos_trial.state}';
hands                               =vertcat(Pos_trial.hand)+1;
choices                             =vertcat(Pos_trial.choice)+1;
FR                                  =vertcat(Pos_trial.FR);
[Positions,~,Position_idx]          =unique([Pos_trial.position]');
[Fixation,~,Fixation_idx]           =unique([Pos_trial.fixation]');


[~, ~, StateIndex]          =unique(States);
u_cho                       =unique(choices)';
u_hnd                       =unique(hands)';
INCHnamepart={'in','ch'};
LHRHnamepart={'NH','LH','RH'};

%% epochs_for_multicomparison is something different actually
states_for_multicomparison=keys.epoch_for_multicomparison(ismember(keys.epoch_for_multicomparison,keys.EPOCHS(:,1)));

multicomp_states=states_for_multicomparison(ismember(states_for_multicomparison,States'));
labelsp={'LS','-','RS'};
labelIN={'IN','-','CH'};
true_labels={'false','true'};
%true_labels={'false','true'};

%% Independently for Instructed and Choice !!
for ch=u_cho
    for hn=u_hnd
        idx=choices==ch & hands==hn;
        if numel(Fixation)==1
            varnames={'State','position'};
            [anova_out.p,anova_out.table,anova_out.stats,anova_out.terms] =...
                anovan(FR(idx),[StateIndex(idx),Position_idx(idx)],'model','full','varnames',varnames,'display',keys.display_anova_tables);
        else
            varnames={'State','position','fixation'};
            [anova_out.p,anova_out.table,anova_out.stats,anova_out.terms] =...
                anovan(FR(idx),[StateIndex(idx),Position_idx(idx),Fixation_idx(idx)],'model','full','varnames',varnames,'display',keys.display_anova_tables);
            anova_struct.([INCHnamepart{ch} '_' LHRHnamepart{hn} '_fixation_main_effect'])=anova_out.p(3)<0.05; %main effect on epoch!
            
        end
        
        anova_struct.([INCHnamepart{ch} '_' LHRHnamepart{hn} '_epoch_main_effect'])=anova_out.p(1)<0.05; %main effect on epoch!
        anova_struct.([INCHnamepart{ch} '_' LHRHnamepart{hn} '_position_main_effect'])=anova_out.p(2)<0.05; %main effect on epoch!
        
        
        for s=multicomp_states(:)'
            idx=ismember(States,s) & choices==ch & hands==hn;
            if any(~isnan(FR(idx)))
                if numel(Fixation)==1
                    
                    idx_L= ismember(Position_idx,find(real(Positions)<0));
                    idx_R= ismember(Position_idx,find(real(Positions)>0));
                    [anova_outs.p] = anovan(FR(idx),{Position_idx(idx)},'model','full','varnames',{'Positions'},'display',keys.display_anova_tables);
                    h=anova_outs.p(1)<0.05;
                    labelindexsp=h*sign(nanmean(FR(idx & idx_R)) - nanmean(FR(idx & idx_L)))+2; labelindexsp(isnan(labelindexsp))=2;
                    anova_struct.([INCHnamepart{ch} '_' LHRHnamepart{hn} '_' s{:} '_position'])=labelsp{labelindexsp};
                else
                    % excluding positions with only one fixation point and
                    % vice versa % --- DOESNT WORK
                    %                     unique_combinations=unique([Fixation_idx(idx),Position_idx(idx)],'rows');
                    %                     potential_indexes=1:max([Fixation_idx(idx);Position_idx(idx)]);
                    %                     frequency_pos=hist(unique_combinations(:,2),potential_indexes);
                    %                     frequency_fix=hist(unique_combinations(:,1),potential_indexes);
                    %                     idx_p=ismember(Position_idx,potential_indexes(frequency_pos>2));
                    %                     idx_f=ismember(Fixation_idx,potential_indexes(frequency_fix>1));
                    %                     idx=idx & idx_p & idx_f;
                    
                    [anova_outs.p] = anovan(FR(idx),[Position_idx(idx),Fixation_idx(idx)],'model','full','varnames',{'Positions','Fixations'},'display',keys.display_anova_tables);
                    %[anova_outs.p] = anovan(FR(idx),[Fixation_idx(idx)],'model','full','varnames',{'Fixations'},'display',keys.display_anova_tables);
                    if isnan(anova_outs.p(2))
                        %isnan(anova_outs.p(1)) ||
                        [temp_p] = anovan(FR(idx),[Position_idx(idx),Fixation_idx(idx)],'model',[1 0;1 1],'varnames',{'Positions','Fixations'},'display',keys.display_anova_tables);
                        %                     [anova_outs.p(1)] = anovan(FR(idx),Position_idx(idx),'model','full','varnames',{'Positions'},'display',keys.display_anova_tables);
                        %                       [anova_outs.p(2)] = anovan(FR(idx),Fixation_idx(idx),'model','full','varnames',{'Fixations'},'display',keys.display_anova_tables);
                        %                        anova_outs.p(3) = NaN;
                        anova_outs.p(1)=temp_p(1);
                        anova_outs.p(3)=temp_p(2);
                    end
                    h1=anova_outs.p(1)<0.05;
                    h2=anova_outs.p(2)<0.05;
                    h3=anova_outs.p(3)<0.05;
                    
                    anova_struct.([INCHnamepart{ch} '_' LHRHnamepart{hn} '_' s{:} '_position'])=true_labels{h1+1};
                    anova_struct.([INCHnamepart{ch} '_' LHRHnamepart{hn} '_' s{:} '_fixation'])=true_labels{h2+1};
                    anova_struct.([INCHnamepart{ch} '_' LHRHnamepart{hn} '_' s{:} '_PxF'])=true_labels{h3+1};
                    
                    anova_struct.([INCHnamepart{ch} '_' LHRHnamepart{hn} '_' s{:} '_position_PV'])=anova_outs.p(1);
                    anova_struct.([INCHnamepart{ch} '_' LHRHnamepart{hn} '_' s{:} '_fixation_PV'])=anova_outs.p(2);
                    anova_struct.([INCHnamepart{ch} '_' LHRHnamepart{hn} '_' s{:} '_PxF_PV'])=anova_outs.p(3);
                    
                    if strcmp(s{:},'Fhol')
                       fix1=nanmean(FR(Fixation_idx==1&idx));
                       fix2=nanmean(FR(Fixation_idx==2&idx));
                       fix3=nanmean(FR(Fixation_idx==3&idx));
                       if fix3>fix2 && fix2>fix1
                            anova_struct.([INCHnamepart{ch} '_' LHRHnamepart{hn} '_' s{:} '_monotonous'])='RS';
                       elseif fix3<fix2 && fix2<fix1
                            anova_struct.([INCHnamepart{ch} '_' LHRHnamepart{hn} '_' s{:} '_monotonous'])='LS';
                       else
                            anova_struct.([INCHnamepart{ch} '_' LHRHnamepart{hn} '_' s{:} '_monotonous'])='-';
                       end
                    end
                end
                
                %                 %% choice part for position with strongest response
                %                 if ch==1 && ismember(2,u_cho) %%instructed
                %                     for p=1:max(Position_idx)
                %                         AAA(p)=nanmean(FR(Position_idx == p & idx));
                %                     end
                %                     [~,Position_with_maximum_response_in_instructed(hn,ismember(multicomp_states,s))]=max(AAA);
                %                 elseif ch==2  && ismember(1,u_cho) %% choice  assuming that ch== 1 happened before
                %                     maximum_response_position=Position_with_maximum_response_in_instructed(hn,ismember(multicomp_states,s));
                %
                %                     idx_IN=ismember(States,s) & choices==1 & hands==hn & Position_idx==maximum_response_position;
                %                     idx_CH=ismember(States,s) & choices==2 & hands==hn & Position_idx==maximum_response_position;
                %
                %                     [h,~]=ttest2(FR(idx_IN),FR(idx_CH));
                %                     ES=nanmean(FR(idx_CH))-nanmean(FR(idx_IN));
                %                     labelindex=h*sign(ES)+2; labelindex(isnan(labelindex))=2;
                %                     anova_struct.([INCHnamepart{ch} '_' LHRHnamepart{hn} '_' s{:} '_RF_choice1'])=labelIN{labelindex};
                %
                %                     % combining both targets for choice
                %                     Opposite_position=(-1)*real(Positions(maximum_response_position))+1i*imag(Positions(maximum_response_position));
                %                     [~, opposite_position_index]=min(Positions-Opposite_position);
                %                     idx_IN=ismember(States,s) & choices==1 & hands==hn & Position_idx==maximum_response_position;
                %                     idx_CH=ismember(States,s) & choices==2 & hands==hn & (Position_idx==maximum_response_position | Position_idx==opposite_position_index);
                %                     [h,~]=ttest2(FR(idx_IN),FR(idx_CH));
                %                     ES=nanmean(FR(idx_CH))-nanmean(FR(idx_IN));
                %                     labelindex=h*sign(ES)+2; labelindex(isnan(labelindex))=2;
                %                     anova_struct.([INCHnamepart{ch} '_' LHRHnamepart{hn} '_' s{:} '_RF_choice2'])=labelIN{labelindex};
                %
                %                 end
                
                %% choice part for position with strongest response
                
                idx_L= ismember(Position_idx,find(real(Positions)<0));
                idx_R= ismember(Position_idx,find(real(Positions)>0));
                idxLR= [idx_L,false(size(idx_L)),idx_R];
                if ch==1 && ismember(2,u_cho) %%instructed
                    Position_with_maximum_response_in_instructed(hn,ismember(multicomp_states,s))=sign(nanmean(FR(idx_R & idx)) - nanmean(FR(idx_L & idx)));
                elseif ch==2  && ismember(1,u_cho) %% choice  assuming that ch== 1 happened before
                    maximum_response_position=Position_with_maximum_response_in_instructed(hn,ismember(multicomp_states,s))+2;
                    maximum_response_position(isnan(maximum_response_position))=2;
                    idx_IN=ismember(States,s) & choices==1 & hands==hn & idxLR(:,maximum_response_position);
                    idx_CH=ismember(States,s) & choices==2 & hands==hn & idxLR(:,maximum_response_position);
                    
                    [h,~]=ttest2(FR(idx_IN),FR(idx_CH));
                    ES=nanmean(FR(idx_CH))-nanmean(FR(idx_IN));
                    labelindex=h*sign(ES)+2; labelindex(isnan(labelindex))=2;
                    anova_struct.([INCHnamepart{ch} '_' LHRHnamepart{hn} '_' s{:} '_RF_choice1'])=labelIN{labelindex};
                    
                    % combining both targets for choice
                    Opposite_position=(-1)*real(Positions(maximum_response_position))+1i*imag(Positions(maximum_response_position));
                    [~, opposite_position_index]=min(Positions-Opposite_position);
                    idx_IN=ismember(States,s) & choices==1 & hands==hn & idxLR(:,maximum_response_position);
                    idx_CH=ismember(States,s) & choices==2 & hands==hn;
                    [h,~]=ttest2(FR(idx_IN),FR(idx_CH));
                    ES=nanmean(FR(idx_CH))-nanmean(FR(idx_IN));
                    labelindex=h*sign(ES)+2; labelindex(isnan(labelindex))=2;
                    anova_struct.([INCHnamepart{ch} '_' LHRHnamepart{hn} '_' s{:} '_RF_choice2'])=labelIN{labelindex};
                    
                end
                
            else
                anova_struct.([INCHnamepart{ch} '_' LHRHnamepart{hn} '_' s{:} '_position'])=labelsp{2};
            end
            
            
            
            
        end
    end
end
end

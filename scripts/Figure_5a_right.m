clear all;
close all;
clc;
%% Figures 5a right, from Unimodal and bimodal access to WM (2019) Wolff et al.
cd('D:\UBA_WM') %path to main dir.
addpath(genpath(cd))

params.n_folds=8; % number of training and testing sets
params.reps=50; % number of reptitions
params.toi=[.1 .4]; % time-window of interest
params.span=10; % width of each segment within sliding-window (in ms)
params.hz=500; % sample rate of data
%
do_decoding=0; %1=run the decodng (takes long time!) or 0= load in previous output
%%
if do_decoding
    for sub=1:30
        fprintf(['Doing ' num2str(sub) '\n'])
        
        load(['AUD_exp_' num2str(sub) '.mat'])
        
        cue_cond=Results(:,3);
        
        incl_a=setdiff(1:size(Results,1), ft_imp_aud.bad_trials);
        incl_v=setdiff(1:size(Results,1), ft_imp_vis.bad_trials);
        incl=intersect(incl_a, incl_v);
        
        span=params.span/(1000/params.hz);
        
        %%%%%%%%%%%%
        % visual impulse
        data_v=ft_imp_vis.trial(incl,:,ft_imp_vis.time>params.toi(1)&ft_imp_vis.time<=params.toi(2));
        data_v=bsxfun(@minus,data_v,mean(data_v,3)); % mean center over whole window
        data_v=movmean(data_v,span,3,'Endpoints','discard');
        data_v=data_v(:,:,1:span:end);
        data_v=reshape(data_v,[size(data_v,1),size(data_v,2)*size(data_v,3)]);
        
        % auditory impulse
        data_a=ft_imp_aud.trial(incl,:,ft_imp_aud.time>params.toi(1)&ft_imp_aud.time<=params.toi(2));
        data_a=bsxfun(@minus,data_a,mean(data_a,3)); % mean center over whole window
        data_a=movmean(data_a,span,3,'Endpoints','discard');
        data_a=data_a(:,:,1:span:end);
        data_a=reshape(data_a,[size(data_a,1),size(data_a,2)*size(data_a,3)]);
        
        % cued item; train on visual impulse, test on auditory impulse
        conds=round(log2((Results(incl,4))./270).*2)+1; % cued item frequency, convert to condition labels
        temp_b=nan(params.reps,1);
        for rep=1:params.reps
            distance_b = mahal_func_ordinal_kfold_b_sep(data_a,data_v,conds,params.n_folds);
            temp_b(rep,1)=mean(distance_b,1);
        end
        cued_imp_aud_v_dec(sub,1)=mean(temp_b,1);
        
        % cued item; train on auditory impulse, test on visual impulse
        conds=round(log2((Results(incl,4))./270).*2)+1; % cued item frequency, convert to condition labels
        temp_b=nan(params.reps,1);
        for rep=1:params.reps
            distance_b = mahal_func_ordinal_kfold_b_sep(data_v,data_a,conds,params.n_folds);
            temp_b(rep,1)=mean(distance_b,1);
        end
        cued_imp_vis_a_dec(sub,1)=mean(temp_b,1);
    end
    cued_imps_dec_b=(cued_imp_aud_v_dec+cued_imp_vis_a_dec)./2; % average
else
    load('Figure_5a_right_results')
end
%% significance testing
n_perms=100000;
p_imps_cued=GroupPermTest(cued_imps_dec_b,100000,2);
%%
ci_imps_cued=bootci(100000,@mean,cued_imps_dec_b);
%% plot boxplots and error bars
pos=[1];
figure
title('Auditory WM')
hold all
b1=boxplot([cued_imps_dec_b],...
    'positions',pos,'Widths',0.1,'Symbol','ko','Labels',{'aud. imp. & vis. imp.'});
set(findobj(gcf,'LineStyle','--'),'LineStyle','-')
set(b1(:,1),'color','b');
plot(pos(1),mean(cued_imps_dec_b,1),'o','MarkerFaceColor','b','MarkerEdgeColor','none','MarkerSize',10)
plot([pos(1) pos(1)],ci_imps_cued','Color','b','LineWidth',3)
plot([0 4],[0 0 ],'Color','k','LineWidth',.5,'LineStyle','--')
ylim([-.005 .005])
ylabel('Cross-generalization')
set(gca,'TickDir','out')

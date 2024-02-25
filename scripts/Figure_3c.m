clear all;
close all;
clc;
%% Figures 3c, from Unimodal and bimodal acces to WM (2019) Wolff et al.
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
        
        span=params.span/(1000/params.hz);
        
        %%%%%%%%%%%%
        % auditory impulse
        incl=setdiff(1:size(Results,1),ft_imp_aud.bad_trials); % good trials to be included
        data=ft_imp_aud.trial(incl,:,ft_imp_aud.time>params.toi(1)&ft_imp_aud.time<=params.toi(2));
        data=bsxfun(@minus,data,mean(data,3)); % mean center over whole window
        data=movmean(data,span,3,'Endpoints','discard');
        data=data(:,:,1:span:end);
        data=reshape(data,[size(data,1),size(data,2)*size(data,3)]);
        
        % cued item
        conds=round(log2((Results(incl,4))./270).*2)+1; % cued item frequency, convert to condition labels
        temp_b=nan(params.reps,1);
        for rep=1:params.reps
            distance_b = mahal_func_ordinal_kfold_b(data,conds,params.n_folds);
            temp_b(rep,1)=mean(distance_b,1);
        end
        cued_imp_aud_dec(sub,1)=mean(temp_b,1);
        
        % uncued item
        conds=round(log2((Results(incl,5))./270).*2)+1; % cued item frequency, convert to condition labels
        temp_b=nan(params.reps,1);
        for rep=1:params.reps
            distance_b = mahal_func_ordinal_kfold_b(data,conds,params.n_folds);
            temp_b(rep,1)=mean(distance_b,1);
        end
        uncued_imp_aud_dec(sub,1)=mean(temp_b,1);
        
        %%%%%%%%%%%%
        % visual impulse
        incl=setdiff(1:size(Results,1),ft_imp_vis.bad_trials); % good trials to be included
        data=ft_imp_vis.trial(incl,:,ft_imp_vis.time>params.toi(1)&ft_imp_vis.time<=params.toi(2));
        data=bsxfun(@minus,data,mean(data,3));
        data=movmean(data,span,3,'Endpoints','discard');
        data=data(:,:,1:span:end);
        data=reshape(data,[size(data,1),size(data,2)*size(data,3)]);
        
        % cued item
        conds=round(log2((Results(incl,4))./270).*2)+1; % cued item frequency, convert to condition labels
        temp_b=nan(params.reps,1);
        for rep=1:params.reps
            distance_b = mahal_func_ordinal_kfold_b(data,conds,params.n_folds);
            temp_b(rep,1)=mean(distance_b,1);
        end
        cued_imp_vis_dec(sub,1)=mean(temp_b,1);
        
        % uncued item
        conds=round(log2((Results(incl,5))./270).*2)+1; % cued item frequency, convert to condition labels
        temp_b=nan(params.reps,1);
        for rep=1:params.reps
            distance_b = mahal_func_ordinal_kfold_b(data,conds,params.n_folds);
            temp_b(rep,1)=mean(distance_b,1);
        end
        uncued_imp_vis_dec(sub,1)=mean(temp_b,1);
    end
else
    load('Figure_3c_results')
end
%% significance testing
n_perms=100000;
p_aud_cued=GroupPermTest(cued_imp_aud_dec,n_perms,1);
p_aud_uncued=GroupPermTest(uncued_imp_aud_dec,n_perms,1);
p_vis_cued=GroupPermTest(cued_imp_vis_dec,n_perms,1);
p_vis_uncued=GroupPermTest(uncued_imp_vis_dec,n_perms,1);
%% make table for JASP, for Bayesian analyses
AUD_WM_imps_table=table(cued_imp_aud_dec,uncued_imp_aud_dec,...
    cued_imp_vis_dec,uncued_imp_vis_dec,...
    'VariableNames',{'imp_aud_cued','imp_aud_uncued','imp_vis_cued','imp_vis_uncued'});
writetable(AUD_WM_imps_table,fullfile([pwd '\results\'],'AUD_WM_imps_table.txt'),'Delimiter',' ')
%% make C.I. for plots
ci_impv_cued=bootci(n_perms,@mean,cued_imp_vis_dec);
ci_impv_uncued=bootci(n_perms,@mean,uncued_imp_vis_dec);
ci_impa_cued=bootci(n_perms,@mean,cued_imp_aud_dec);
ci_impa_uncued=bootci(n_perms,@mean,uncued_imp_aud_dec);
%% Figure 3c
pos=[.9 1.15 1.55 1.8];
figure
title('Auditory WM')
hold all
b1=boxplot([cued_imp_aud_dec,uncued_imp_aud_dec,cued_imp_vis_dec,uncued_imp_vis_dec],...
    'positions',pos,'Widths',0.1,'Symbol','ko','Labels',{'aud. cued','aud. uncued','vis. cued','vis. uncued'});
set(findobj(gcf,'LineStyle','--'),'LineStyle','-')
set(b1(:,1),'color','b');set(b1(:,2),'color','k');set(b1(:,3),'color','b');set(b1(:,4),'color','k');
plot(pos(1),mean(cued_imp_aud_dec,1),'o','MarkerFaceColor','b','MarkerEdgeColor','none','MarkerSize',10)
plot(pos(2),mean(uncued_imp_aud_dec,1),'o','MarkerFaceColor','k','MarkerEdgeColor','none','MarkerSize',10)
plot(pos(3),mean(cued_imp_vis_dec,1),'o','MarkerFaceColor','b','MarkerEdgeColor','none','MarkerSize',10)
plot(pos(4),mean(uncued_imp_vis_dec,1),'o','MarkerFaceColor','k','MarkerEdgeColor','none','MarkerSize',10)
plot([pos(1) pos(1)],ci_impa_cued','Color','b','LineWidth',3)
plot([pos(2) pos(2)],ci_impa_uncued','Color','k','LineWidth',3)
plot([pos(3) pos(3)],ci_impv_cued','Color','b','LineWidth',3)
plot([pos(4) pos(4)],ci_impv_uncued','Color','k','LineWidth',3)
plot([0 4],[0 0 ],'Color','k','LineWidth',.5,'LineStyle','--')
ylim([-.002 .0025])
ylabel('Decoding accuracy')
set(gca,'TickDir','out')
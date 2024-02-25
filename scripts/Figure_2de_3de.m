clear all;
close all;
clc;
%% Figures 2de, 3de, from Unimodal and bimodal acces to WM (2019) Wolff et al.
cd('D:\UBA_WM') %path to main dir.
addpath(genpath(cd))

params.n_folds=8; % number of training and testing sets
params.reps=50; % number of repititions
params.steps=8; % steps sliding-window moves in (in ms)
params.span=10; % width of each segment within sliding-window (in ms)
params.w_length=100; % length of sliding window (in ms)
params.s_factor=16; % sd of guassian smoothing kernel for decoding time-course (in ms)
params.hz=500; % sample rate of data
%
do_decoding=0; %1=run the decodng (takes long time!) or 0= load in previous output
%%
if do_decoding
    for sub=1:28
        fprintf(['Doing ' num2str(sub) '\n'])
        load(['VIS_exp_' num2str(sub) '.mat'])
        
        % Memory item 1
        params.toi=[-.053 .9]; % time-window of interest
        incl=setdiff(1:size(Results,1),ft_mem1.bad_trials); % good trials to be included
        data=ft_mem1.trial(incl,:,:);
        
        params.time_dat=ft_mem1.time;
        params.theta=Results(incl,1).*2; %
        
        [cos_d,distances,time_mem]=decoding_dynamics_time_course_vis(data,params);
        
        mem1_dec(sub,:)=cos_d;
        mem1_dists(sub,:,:)=distances;
        
        mem1_dec_temp=cos_d;
        mem1_dists_temp=distances;
        
        % Memory item 2
        incl=setdiff(1:size(Results,1),ft_mem2.bad_trials); % good trials to be included
        data=ft_mem2.trial(incl,:,:);
        params.theta=Results(incl,2).*2; %
        
        [cos_d,distances]=decoding_dynamics_time_course_vis(data,params);
        
        mem2_dec(sub,:)=cos_d;
        mem2_dists(sub,:,:)=distances;
        
        mem2_dec_temp=cos_d;
        mem2_dists_temp=distances;
        
        %%%%%%%%%%%%
        % auditory impulse
        params.toi=[-0.053 .5];
        params.time_dat=ft_imp_aud.time;
        incl=setdiff(1:size(Results,1),ft_imp_aud.bad_trials); % good trials to be included
        data=ft_imp_aud.trial(incl,:,:);
        
        % cued item
        params.theta=Results(incl,4).*2;
        [cos_d,distances,time_imp]=decoding_dynamics_time_course_vis(data,params);
        cued_imp_aud_dec(sub,:)=cos_d;
        cued_imp_aud_dists(sub,:,:)=distances;
        
        cued_imp_aud_dec_temp=cos_d;
        cued_imp_aud_dists_temp=distances;
        
        % uncued item
        params.theta=Results(incl,5).*2;
        [cos_d,distances]=decoding_dynamics_time_course_vis(data,params);
        uncued_imp_aud_dec(sub,:)=cos_d;
        uncued_imp_aud_dists(sub,:,:)=distances;
        
        uncued_imp_aud_dec_temp=cos_d;
        uncued_imp_aud_dists_temp=distances;
        
        %%%%%%%%%%%%
        % visual impulse
        params.toi=[-0.053 .5];
        params.time_dat=ft_imp_vis.time;
        incl=setdiff(1:size(Results,1),ft_imp_vis.bad_trials); % good trials to be included
        data=ft_imp_vis.trial(incl,:,:);
        
        % cued item
        params.theta=Results(incl,4).*2;
        [cos_d,distances]=decoding_dynamics_time_course_vis(data,params);
        cued_imp_vis_dec(sub,:)=cos_d;
        cued_imp_vis_dists(sub,:,:)=distances;
        
        cued_imp_vis_dec_temp=cos_d;
        cued_imp_vis_dists_temp=distances;
        
        % uncued item
        params.theta=Results(incl,5).*2;
        [cos_d,distances]=decoding_dynamics_time_course_vis(data,params);
        uncued_imp_vis_dec(sub,:)=cos_d;
        uncued_imp_vis_dists(sub,:,:)=distances;
        
        uncued_imp_vis_dec_temp=cos_d;
        uncued_imp_vis_dists_temp=distances;
        
    end
else
    load('Figure_2de_3de_results')
end
ang_space=-67.5:22.5:90;
%% significance testing
n_perms=100000;
% Memory items
[datobs, datrnd] = cluster_test_helper(mem1_dec',n_perms); % randomize data at each time-point
[~, p_mem1] = cluster_test(datobs,datrnd,1,0.05,0.05); % do cluster-corrected test (one-sided)

[datobs, datrnd] = cluster_test_helper(mem2_dec',n_perms); % randomize data at each time-point
[~, p_mem2] = cluster_test(datobs,datrnd,1,0.05,0.05); % do cluster-corrected test (one-sided)

% auditory impulse
[datobs, datrnd] = cluster_test_helper(cued_imp_aud_dec',n_perms); % randomize data at each time-point
[~, p_cued_imp_aud] = cluster_test(datobs,datrnd,1,0.05,0.05); % do cluster-corrected test (one-sided)

[datobs, datrnd] = cluster_test_helper(uncued_imp_aud_dec',n_perms); % randomize data at each time-point
[~, p_uncued_imp_aud] = cluster_test(datobs,datrnd,1,0.05,0.05); % do cluster-corrected test (one-sided)

% visual impulse
[datobs, datrnd] = cluster_test_helper(cued_imp_vis_dec',n_perms); % randomize data at each time-point
[~, p_cued_imp_vis] = cluster_test(datobs,datrnd,1,0.05,0.05); % do cluster-corrected test (one-sided)

[datobs, datrnd] = cluster_test_helper(uncued_imp_vis_dec',n_perms); % randomize data at each time-point
[~, p_uncued_imp_vis] = cluster_test(datobs,datrnd,1,0.05,0.05); % do cluster-corrected test (one-sided)
%% Figure 2d, left
clims=[-0.0075 0.015];
figure;
imagesc(time_mem,ang_space,squeeze(mean(mem1_dists,1)),clims)
xlim([-0.05 .9])
colorbar
pbaspect([2,0.5,0.5]);set(gca,'TickDir','out');ax = gca;ax.YTick = ang_space(2:2:end);ax.XTick=linspace(-.1,.9,11);
ylabel('Orientation difference (degrees)')
xlabel('Time (s) - relative to item 1 onset')
colormap 'hot'
title('item 1')
%% Figure 2d, right
clims=[-0.0075 0.015];
figure;
imagesc(time_mem,ang_space,squeeze(mean(mem2_dists,1)),clims)
xlim([-0.05 .9])
colorbar
pbaspect([2,0.5,0.5]);set(gca,'TickDir','out');ax = gca;ax.YTick = ang_space(2:2:end);ax.XTick=linspace(-.1,.9,11);
ylabel('Orientation difference (degrees)')
xlabel('Time (s) - relative to item 2 onset')
colormap 'hot'
title('item 2')
%%
ci_mem1=bootci(n_perms,@mean,(mem1_dec)); % get confidence interval for plotting
ci_mem2=bootci(n_perms,@mean,(mem2_dec)); % get confidence interval for plotting
%% Figure 2e, left
pclustu_mem1 = unique(p_mem1);
npclust_mem1 = nnz(pclustu_mem1 < 0.05);
figure
title('Visual WM')
hold all
plot(time_mem,mean(mem1_dec(:,:),1),'Color',[1 .6 0],'LineWidth',4)
line('XData', [-0.05 .9], 'YData', [0 0], 'LineStyle', '--','LineWidth', 1, 'Color','k');
fill([0,0,.2,.2],[-0.001,-0.0008,-0.0008,-0.001],[0 0 0],'EdgeColor','none')
plot(time_mem,ci_mem1,'Color',[1 .6 0 1],'LineWidth',.5, 'LineStyle', ':')
for ipclust = 1:npclust_mem1
    currind  = p_mem1 == pclustu_mem1(ipclust);
    fill([min(time_mem(currind)),min(time_mem(currind)),max(time_mem(currind)),max(time_mem(currind))],[0.0075,0.0074,0.0074,0.0075],[1 .6 0],'EdgeColor','none')
    h=fill([min(time_mem(currind)),min(time_mem(currind)),max(time_mem(currind)),max(time_mem(currind))],[0,0.008,0.008,0],[1 .6 0],'EdgeColor','none');
    set(h,'facealpha',0.1);
end
xlim([-.05 .9]);ylim([-0.00075 0.0075]); % limits od x and y axes
set(gca,'TickDir','out')
xlabel('Time (s) - relative to memory item 1');ylabel('Decoding accuracy')
legend('item 1','Location','northeast')
ax=gca;ax.XTick=linspace(-.1,.9,11);pbaspect([1,1,1]) % ratio between axes
title('item 1')
%% Figure 2e, right
pclustu_mem2 = unique(p_mem2);
npclust_mem2 = nnz(pclustu_mem2 < 0.05);
figure
title('Visual WM')
hold all
plot(time_mem,mean(mem2_dec(:,:),1),'Color',[.7 .33 0],'LineWidth',4)
line('XData', [-0.05 .9], 'YData', [0 0], 'LineStyle', '--','LineWidth', 1, 'Color','k');
fill([0,0,.2,.2],[-0.001,-0.0008,-0.0008,-0.001],[0 0 0],'EdgeColor','none')
plot(time_mem,ci_mem2,'Color',[.7 .33 0 1],'LineWidth',.5, 'LineStyle', ':')
for ipclust = 1:npclust_mem2
    currind  = p_mem2 == pclustu_mem2(ipclust);
    fill([min(time_mem(currind)),min(time_mem(currind)),max(time_mem(currind)),max(time_mem(currind))],[0.0075,0.0074,0.0074,0.0075],[.7 .33 0],'EdgeColor','none')
    h=fill([min(time_mem(currind)),min(time_mem(currind)),max(time_mem(currind)),max(time_mem(currind))],[0,0.008,0.008,0],[.7 .33 0],'EdgeColor','none');
    set(h,'facealpha',0.1);
end
xlim([-.05 .9]);ylim([-0.00075 0.0075]); % limits od x and y axes
set(gca,'TickDir','out')
xlabel('Time (s) - relative to memory item 2');ylabel('Decoding accuracy')
legend('item 2','Location','northeast')
ax=gca;ax.XTick=linspace(-.1,.9,11);pbaspect([1,1,1]) % ratio between axes
title('item 2')
%% Figure 3d left
clims=[-0.0015 0.003];
figure
subplot(2,1,1)
imagesc(time_imp,ang_space,squeeze(mean(cued_imp_aud_dists,1)),clims)
xlim([-0.05 .5])
colorbar
pbaspect([2,0.5,0.5]);set(gca,'TickDir','out');ax = gca;ax.YTick = ang_space(2:2:end);ax.XTick=linspace(-.1,.5,7);
ylabel('Orientation difference (degrees)')
xlabel('Time (s) - relative to auditory impulse')
colormap 'hot'
title('cued item')

subplot(2,1,2)
imagesc(time_imp,ang_space,squeeze(mean(uncued_imp_aud_dists,1)),clims)
xlim([-0.05 .5])
colorbar
pbaspect([2,0.5,0.5]);set(gca,'TickDir','out');ax = gca;ax.YTick = ang_space(2:2:end);ax.XTick=linspace(-.1,.5,7);
ylabel('Orientation difference (degrees)')
xlabel('Time (s) - relative to auditory impulse')
colormap 'hot'
title('uncued item')
%% Figure 3a right
clims=[-0.0015 0.003];
figure;
subplot(2,1,1)
imagesc(time_imp,ang_space,squeeze(mean(cued_imp_vis_dists(:,:,:),1)),clims)
xlim([-0.05 .5])
colorbar
pbaspect([2,0.5,0.5]);set(gca,'TickDir','out');ax = gca;ax.YTick = ang_space(2:2:end);ax.XTick=linspace(-.1,.5,7);
ylabel('Orientation difference (degrees)')
xlabel('Time (s) - relative to visual impulse')
colormap 'hot'
title('cued item')

subplot(2,1,2)
imagesc(time_imp,ang_space,squeeze(mean(uncued_imp_vis_dists,1)),clims)
xlim([-0.05 .5])
colorbar
pbaspect([2,0.5,0.5]);set(gca,'TickDir','out');ax = gca;ax.YTick = ang_space(2:2:end);ax.XTick=linspace(-.1,.5,7);
ylabel('Orientation difference (degrees)')
xlabel('Time (s) - relative to visual impulse')
colormap 'hot'
title('uncued item')
%%
ci_cued_imp_aud=bootci(n_perms,@mean,(cued_imp_aud_dec));
ci_uncued_imp_aud=bootci(n_perms,@mean,(uncued_imp_aud_dec));

ci_cued_imp_vis=bootci(n_perms,@mean,(cued_imp_vis_dec));
ci_uncued_imp_vis=bootci(n_perms,@mean,(uncued_imp_vis_dec));
%% Figure 3e left
pclustu_impa_cued = unique(p_cued_imp_aud);
npclust_impa_cued = nnz(pclustu_impa_cued < 0.05);
pclustu_impa_uncued = unique(p_uncued_imp_aud);
npclust_impa_uncued = nnz(pclustu_impa_uncued < 0.05);
figure
hold all
plot(time_imp,mean(cued_imp_aud_dec,1),'Color',[0 0 1],'LineWidth',4)
plot(time_imp,mean(uncued_imp_aud_dec,1),'Color',[0 0 0],'LineWidth',4)
line('XData', [-0.05 .5], 'YData', [0 0], 'LineStyle', '--','LineWidth', 1, 'Color','k');
fill([0,0,.1,.1],[-0.001,-0.0008,-0.0008,-0.001],[0 0 0],'EdgeColor','none')
plot(time_imp,ci_cued_imp_aud,'Color',[0 0 1 1],'LineWidth',.5, 'LineStyle', ':')
plot(time_imp,ci_uncued_imp_aud,'Color',[0 0 0 1],'LineWidth',.5, 'LineStyle', ':')
for ipclust = 1:npclust_impa_cued
    currind  = p_cued_imp_aud == pclustu_impa_cued(ipclust);
    fill([min(time_imp(currind)),min(time_imp(currind)),max(time_imp(currind)),max(time_imp(currind))],[0.002,0.00195,0.00195,0.002],[0 0 1],'EdgeColor','none')
    h=fill([min(time_imp(currind)),min(time_imp(currind)),max(time_imp(currind)),max(time_imp(currind))],[0,0.002,0.002,0],[0 0 1],'EdgeColor','none');
    set(h,'facealpha',0.1);
end
for ipclust = 1:npclust_impa_uncued
    currind  = p_uncued_imp_aud == pclustu_impa_uncued(ipclust);
    fill([min(time_imp(currind)),min(time_imp(currind)),max(time_imp(currind)),max(time_imp(currind))],[0.0019,0.00185,0.00185,0.0019],[0 0 0],'EdgeColor','none')
    h=fill([min(time_imp(currind)),min(time_imp(currind)),max(time_imp(currind)),max(time_imp(currind))],[0,0.00185,0.00185,0],[0 0 0],'EdgeColor','none');
    set(h,'facealpha',0.1);
end
xlim([-.05 .5]);ylim([-0.00075 0.002]); % limits od x and y axes
set(gca,'TickDir','out')
xlabel('Time (s) - relative to auditory impulse');ylabel('Decoding accuracy')
legend('cued','uncued','Location','northwest')
ax=gca;ax.XTick=linspace(-.1,.9,11);pbaspect([1,1,1]) % ratio between axes
%% Figure 3e right
pclustu_impv_cued = unique(p_cued_imp_vis);
npclust_impv_cued = nnz(pclustu_impv_cued < 0.05);
pclustu_impv_uncued = unique(p_uncued_imp_vis);
npclust_impv_uncued = nnz(pclustu_impv_uncued < 0.05);
figure
hold all
plot(time_imp,mean(cued_imp_vis_dec,1),'Color',[0 0 1],'LineWidth',4)
plot(time_imp,mean(uncued_imp_vis_dec,1),'Color',[0 0 0],'LineWidth',4)
line('XData', [-0.05 .5], 'YData', [0 0], 'LineStyle', '--','LineWidth', 1, 'Color','k');
fill([0,0,.1,.1],[-0.001,-0.0008,-0.0008,-0.001],[0 0 0],'EdgeColor','none')
plot(time_imp,ci_cued_imp_vis,'Color',[0 0 1 1],'LineWidth',.5, 'LineStyle', ':')
plot(time_imp,ci_uncued_imp_vis,'Color',[0 0 0 1],'LineWidth',.5, 'LineStyle', ':')
for ipclust = 1:npclust_impv_cued
    currind  = p_cued_imp_vis == pclustu_impv_cued(ipclust);
    fill([min(time_imp(currind)),min(time_imp(currind)),max(time_imp(currind)),max(time_imp(currind))],[0.002,0.00195,0.00195,0.002],[0 0 1],'EdgeColor','none')
    h=fill([min(time_imp(currind)),min(time_imp(currind)),max(time_imp(currind)),max(time_imp(currind))],[0,0.002,0.002,0],[0 0 1],'EdgeColor','none');
    set(h,'facealpha',0.1);
end
for ipclust = 1:npclust_impv_uncued
    currind  = p_uncued_imp_vis == pclustu_impv_uncued(ipclust);
    fill([min(time_imp(currind)),min(time_imp(currind)),max(time_imp(currind)),max(time_imp(currind))],[0.0019,0.00185,0.00185,0.0019],[0 0 0],'EdgeColor','none')
    h=fill([min(time_imp(currind)),min(time_imp(currind)),max(time_imp(currind)),max(time_imp(currind))],[0,0.00185,0.00185,0],[0 0 0],'EdgeColor','none');
    set(h,'facealpha',0.1);
end
xlim([-.05 .5]);ylim([-0.00075 0.002]); % limits od x and y axes
set(gca,'TickDir','out')
xlabel('Time (s) - relative to visual impulse');ylabel('Decoding accuracy')
legend('cued','uncued','Location','northwest')
ax=gca;ax.XTick=linspace(-.1,.9,11);pbaspect([1,1,1]) % ratio between axes

%%
function [cos_d,distances,time_dec]=decoding_dynamics_time_course_vis(data,params)

w_length=params.w_length/(1000/params.hz);
span=params.span/(1000/params.hz);
steps=params.steps/(1000/params.hz);
s_factor=(params.s_factor/(1000/params.hz))/steps;
time_dec=params.time_dat(params.time_dat>params.toi(1)&params.time_dat<=params.toi(2));
time_dec=time_dec(1:steps:end);

dat_dec=nan(size(data,1),size(data,2)*(w_length/span),length(time_dec));
for t=1:length(time_dec)
    [~,ind]=min(abs(params.time_dat-time_dec(t)));
    temp_dat=bsxfun(@minus,data(:,:,ind-(w_length-1):ind),...
        mean(data(:,:,ind-(w_length-1):ind),3));
    temp_dat = movmean(temp_dat,span,3,'Endpoints','discard');
    temp_dat=temp_dat(:,:,1:span:end);
    dat_dec(:,:,t)=reshape(temp_dat,[size(temp_dat,1),size(temp_dat,2)*size(temp_dat,3)]);
end
temp_c=nan(params.reps,size(dat_dec,3));
temp_dists=nan(params.reps,length(unique(params.theta)),size(dat_dec,3));
theta=params.theta;
n_folds=params.n_folds;
for rep=1:params.reps
    [distance_c,distances] = mahal_func_theta_kfold_b(dat_dec,theta,n_folds);
    temp_c(rep,:)=mean(distance_c,1);
    temp_dists(rep,:,:)=squeeze(mean(distances,2));
end
cos_d=smoothdata(squeeze(mean(temp_c,1)),2,'gaussian',s_factor*5);
distances=squeeze(mean(temp_dists,1));
end



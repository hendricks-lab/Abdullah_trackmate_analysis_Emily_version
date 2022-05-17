%% Analyze and Plot Rg and MSD per cell
%% Abdullah R. Chaudhary

close all; clear all; clc; 
set(0,'DefaultFigureWindowStyle','docked')
addpath('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/New_General_codes/');
addpath('/Volumes/Emily_htt/New_General_codes/plotSpread/');
addpath('/Volumes/Emily_htt/New_General_codes/raacampbell-notBoxPlot-7d90c27/code/');
addpath('/Volumes/Emily_htt/New_General_codes/raacampbell-notBoxPlot-7d90c27/code/+NBP/');
addpath('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/violin/');
addpath('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/breakyaxis/');
addpath('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/Statistical_testing/');
addpath('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/');
addpath('/Volumes/Emily_htt/AC_codes_epmodified_20200603/')

plus_processive_runs=[];

colour_WT=[0.5 0.5 0.5];
colour_S421D=[0 0.4 0.4];
colour_WThtt=[0.3 0.3 0.3];
colour_S421Dhtt=[0 0.4 0.68];

for k_choose = 1:4

if k_choose == 1    % WT r5
    col1=colour_WT;
    pth1='/Volumes/Emily_htt/Huntingtin_Project/EE/WT/Rab5/WT_r5_rg/';
    pth2='/Volumes/Emily_htt/Huntingtin_Project/EE/WT/Rab5/WT_r5_msd/';
    pth3='/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/Dirbias/WT_r5/';
    addpath(pth1);
    addpath(pth2);
    addpath(pth3);
    fl1='Rg_per_cell.mat';
    fl2='MSD_per_cell.mat';
    fl3='plus_processive_runs.mat';
    an=0;
    DT=0.12;
    kf=0;
    proc_runs=[];
    proc_times=[];
    rng=[];
    tlp=-0.5;
    msd_all=[];
    mean_rg_all=[];
    slp_all=[];
elseif k_choose == 2    % S421D r5
    col1=colour_S421D;
    pth1='/Volumes/Emily_htt/Huntingtin_Project/EE/S421D/Rab5/S421D_r5_rg/';
    pth2='/Volumes/Emily_htt/Huntingtin_Project/EE/S421D/Rab5/S421D_r5_msd/';
    pth3='/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/Dirbias/S421D_r5/';
    addpath(pth1);
    addpath(pth2);
    addpath(pth3);
    fl1='Rg_per_cell.mat';
    fl2='MSD_per_cell.mat';
    fl3='plus_processive_runs.mat';
    an=0;
    DT=0.12;
    kf=0;
    proc_runs=[];
    proc_times=[];
    rng=[];
    tlp=-0.5;
    msd_all=[];
    mean_rg_all=[];
    slp_all=[];

elseif k_choose == 3   % WT +WT htt r5 
    col1=colour_WThtt;
    pth1='/Volumes/Emily_htt/Huntingtin_Project/EE/WT_trans_mchWThtt/Rab5/WT_trans_WThtt_rg/';
    pth2='/Volumes/Emily_htt/Huntingtin_Project/EE/WT_trans_mchWThtt/Rab5/WT_trans_WThtt_msd/';
    pth3='/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/Dirbias/WT+WThtt_r5/';
    addpath(pth1);
    addpath(pth2);
    addpath(pth3);
    fl1='Rg_per_cell.mat';
    fl2='MSD_per_cell.mat';
    fl3='plus_processive_runs.mat';
    an=0;
    DT=0.12;
    kf=0;
    proc_runs=[];
    proc_times=[];
    rng=[];
    tlp=1;
    msd_all=[];
    mean_rg_all=[];
    slp_all=[];
    
elseif k_choose == 4   % WT+S421Dhtt r5
    col1=colour_S421Dhtt;
    pth1='/Volumes/Emily_htt/Huntingtin_Project/EE/WT_trans_S421Dhtt/Rab5/WT_trans_S421Dhtt_rg/';
    pth2='/Volumes/Emily_htt/Huntingtin_Project/EE/WT_trans_S421Dhtt/Rab5/WT_trans_S421Dhtt_msd/';
    pth3='/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/Dirbias/WT+S421Dhtt_r5/';
    addpath(pth1);
    addpath(pth2);
    addpath(pth3);
    fl1='Rg_per_cell.mat';
    fl2='MSD_per_cell.mat';
    fl3='plus_processive_runs.mat';
    an=0;
    DT=0.12;
    kf=0;
    proc_runs=[];
    proc_times=[];
    rng=[];
    tlp=1.5;
    msd_all=[];
    mean_rg_all=[];
    slp_all=[];
end

load(fullfile(pth1,fl1));
load(fullfile(pth2,fl2));
load(fl3); %Emily added

plus_processive_runs{k_choose}=outward_runs_per_cell;

[timek2,msd2] = average_MSD_per_cell_emmodified4kchoose(res1,col1,k_choose);

% MSD Analysis:
for ik=1:numel(res1)
    timek=res1{ik}.timek;
    msd=res1{ik}.msd;
    logt=res1{ik}.logtk;
    logmsd=res1{ik}.log_msd;
    alph=res1{ik}.slp;
    if isnan(msd) %Emily added because the msd had some cells without any data
        k_choose %1 for WT 2 for S421D, helps identify the empty cell
        ik
        warning='Empty cell in msd analysis'
    else
    % Figure(1)
    figure(k_choose*1e1), hold on, 
    p1=plot(timek,msd,'Color',col1);
    p1.Color(4)=0.25;
    xlabel('Time interval (sec)'); ylabel('MSD (\mum^2)'); 
    publication_fig(0,0,1);
    
    % Figure(1)
    figure(k_choose*1e3), hold on, 
    p1=plot(logt,logmsd,'Color',col1);
    p1.Color(4)=0.25;
    xlabel('Log[Time interval (sec)]'); ylabel('Log[MSD (\mum^2)]'); 
    publication_fig(0,0,1);
    
    msd_all=[msd_all,msd];
    slp_all=[slp_all,alph];
    
    rg_per_cell=rg_all{ik};
    if isempty(rg_per_cell) %Emily added because the Rg had some cells without any data
        k_choose %1 for WT 2 for S421D, helps identify the empty cell
        ik
        warning='Empty cell in Rg analysis'
        %clear rg_per_cell
    else
    [mean_rg,pci]=mle(rg_per_cell,'distribution','exp'); 
    mean_rg_all=[mean_rg_all,mean_rg];
    
    %Emily plotting the max likelihood estimation of the rg and to figure out what's wrong
%     rg_x=0:0.03:2.5;
%     rg_hist_per_cell=hist(rg_per_cell,rg_x);
%     rg_hist_mean=hist(mean_rg_all,rg_x);
    
%     figure(k_choose*8), hold on, 
%     bar(rg_x, rg_hist_mean./sum(rg_hist_mean), 'BarWidth',1,'FaceColor',[0, 0.75, 0.75],'edgecolor',[0, 0.75, 0.75],'Facealpha',1);
%     hold on,
%     stairs(rg_x, rg_hist_per_cell./sum(rg_hist_per_cell), 'LineWidth',2,'Color','k');
%     hold on, 
%     xlabel('Radius of Gyration (\mum)'); ylabel('Number of Trajectories');
%     xlim([0 2.1]);
%     publication_fig(0,0,1);

%     figure(k_choose*9), hold on,
%     scatter(1:numel(rg_per_cell),rg_per_cell,'r')
%     plot(1:numel(mean_rg_all),mean_rg_all,'g');
    clear rg_per_cell
    end
    end
%     [mean_rg,pci]=mle(rg_per_cell,'distribution','exp'); %original
%     %Abdullah's code but for some reason my data had some NaNs in it so
%     %we had to add the if: else statement. 
%     mean_rg_all=[mean_rg_all,mean_rg];
%     clear rg_per_cell
end

msd_data{k_choose}=msd_all';
rg_data{k_choose}=mean_rg_all';
timek_f{k_choose}=timek2;
MSD_f{k_choose}=msd2;
alph_f{k_choose}=slp_all; %slopes from the fits to the MSD curve

clear msd_all mean_rg_all rg_per_cell timek2 msd2 slp_all outward_runs_per_cell;
%rmpath(pth1)
%rmpath(pth2)
end

figure('Name','Alpha_per_cell_box','NumberTitle','off'), hold on,
h1=notBoxPlot(alph_f{1},1,'style','line');
set(h1.data,'MarkerFaceColor',colour_WT,'MarkerEdgeColor',colour_WT);
h2=notBoxPlot(alph_f{2},2,'style','line');
set(h2.data,'MarkerFaceColor',colour_S421D,'MarkerEdgeColor',colour_S421D);
h3=notBoxPlot(alph_f{3},3,'style','line');
set(h3.data,'MarkerFaceColor',colour_WThtt,'MarkerEdgeColor',colour_WThtt);
h4=notBoxPlot(alph_f{4},4,'style','line');
set(h4.data,'MarkerFaceColor',colour_S421Dhtt,'MarkerEdgeColor',colour_S421Dhtt);
%set(gca,'Ygrid','on');
xlabel(''); ylabel('\alpha');
xlim([0.5 4.5]);
set(gca,'xticklabel',{[]})
publication_fig(0,0,1)

figure('Name','Rg_per_cell_box','NumberTitle','off'), hold on,
h1=notBoxPlot(rg_data{1},1,'style','line');
set(h1.data,'MarkerFaceColor',colour_WT,'MarkerEdgeColor',colour_WT);
h2=notBoxPlot(rg_data{2},2,'style','line');
set(h2.data,'MarkerFaceColor',colour_S421D, 'MarkerEdgeColor',colour_S421D);
h3=notBoxPlot(rg_data{3},3,'style','line');
set(h3.data,'MarkerFaceColor',colour_WThtt,'MarkerEdgeColor',colour_WThtt);
h4=notBoxPlot(rg_data{4},4,'style','line');
set(h4.data,'MarkerFaceColor',colour_S421Dhtt,'MarkerEdgeColor',colour_S421Dhtt);
%set(gca,'Ygrid','on');
xlabel(''); ylabel('Radius of Gyration (\mum)');
xlim([0.5 4.5]);
set(gca,'xticklabel',{[]})
publication_fig(0,0,1)


figure('Name','Alpha_per_cell_violin','NumberTitle','off'), hold on,
violin(alph_f, 'xlabel',{'WT','S421D','WT+WThtt','WT+S421Dhtt'},'facecolor',[colour_WT;colour_S421D;colour_WThtt;colour_S421Dhtt],'edgecolor','k','mc','k','medc','k--');
publication_fig(0,0,1)
axisHandle = gca; 
    set(findall(gca, 'Type', 'Line'),'LineWidth',2);
    set(gca,'LineWidth',2);
    set(gca,'FontSize',24);
    set(gca, 'FontName', 'Arial');
    set(gca,'XColor',[0 0 0],'YColor',[0 0 0]);
    set(gca,'Box','on');
ylabel('Radius of Gyration (\mum)');
ylim([0.5 2.25]);

figure('Name','Alpha_per_cell_violin_and_box','NumberTitle','off'), hold on,
violin(alph_f, 'xlabel',{'WT','S421D','WT+WThtt','WT+S421Dhtt'},'facecolor',[colour_WT;colour_S421D;colour_WThtt;colour_S421Dhtt],'facealpha',0.3,'edgecolor','k','mc','k','medc','k--');
h1=notBoxPlot_nodotorline(alph_f{1},1,'style','line');
set(h1.data,'MarkerFaceColor','none','MarkerEdgeColor',colour_WT,'MarkerSize', 5);
h2=notBoxPlot_nodotorline(alph_f{2},2,'style','line');
set(h2.data,'MarkerFaceColor','none','MarkerEdgeColor',colour_S421D,'MarkerSize', 5);
h3=notBoxPlot_nodotorline(alph_f{3},3,'style','line');
set(h3.data,'MarkerFaceColor','none','MarkerEdgeColor',colour_WThtt,'MarkerSize', 5);
h4=notBoxPlot_nodotorline(alph_f{4},4,'style','line');
set(h4.data,'MarkerFaceColor','none','MarkerEdgeColor',colour_S421Dhtt,'MarkerSize', 5);
publication_fig(0,0,1)
axisHandle = gca; 
    set(findall(gca, 'Type', 'Line'),'LineWidth',2);
    set(gca,'LineWidth',2);
    set(gca,'FontSize',24);
    set(gca, 'FontName', 'Arial');
    set(gca,'XColor',[0 0 0],'YColor',[0 0 0]);
    set(gca,'Box','on');
ylabel('\alpha');
ylim([0.5 2.25]);
hFig=findall(0,'type','figure');
hLeg=findobj(hFig(1,1),'type','legend');
set(hLeg,'visible','off');

figure('Name','Rg_per_cell_violin','NumberTitle','off'), hold on,
violin(rg_data, 'xlabel',{'WT','S421D','WT+WThtt','WT+S421Dhtt'},'facecolor',[colour_WT;colour_S421D;colour_WThtt;colour_S421Dhtt],'edgecolor','k','bw',0.1,'mc','k','medc','k--');
xlim([0.5 4.5]);
publication_fig(0,0,1)
axisHandle = gca; 
    set(findall(gca, 'Type', 'Line'),'LineWidth',2);
    set(gca,'LineWidth',2);
    set(gca,'FontSize',24);
    set(gca, 'FontName', 'Arial');
    set(gca,'XColor',[0 0 0],'YColor',[0 0 0]);
    set(gca,'Box','on');
ylabel('Radius of Gyration (\mum)');
ylim([-0.25 2]);

figure('Name','Rg_per_cell_violin_and_box','NumberTitle','off'), hold on,
violin(rg_data, 'xlabel',{'WT','S421D','WT+WThtt','WT+S421Dhtt'},'facecolor',[colour_WT;colour_S421D;colour_WThtt;colour_S421Dhtt],'facealpha',0.3,'edgecolor','k','mc','k','medc','k--');
h1=notBoxPlot_nodotorline(rg_data{1},1,'style','line');
set(h1.data,'MarkerFaceColor','none','MarkerEdgeColor',colour_WT, 'MarkerSize', 5);
h2=notBoxPlot_nodotorline(rg_data{2},2,'style','line');
set(h2.data,'MarkerFaceColor','none','MarkerEdgeColor',colour_S421D,'MarkerSize', 5);
h3=notBoxPlot_nodotorline(rg_data{3},3,'style','line');
set(h3.data,'MarkerFaceColor','none','MarkerEdgeColor',[0.3 0.3 0.3],'MarkerSize', 5);
h4=notBoxPlot_nodotorline(rg_data{4},4,'style','line');
set(h4.data,'MarkerFaceColor','none','MarkerEdgeColor',colour_S421Dhtt,'MarkerSize', 5);
xlim([0.5 4.5]);
ylim([-0.25 2]);
publication_fig(0,0,1)
axisHandle = gca; 
    set(findall(gca, 'Type', 'Line'),'LineWidth',2);
    set(gca,'LineWidth',2);
    set(gca,'FontSize',24);
    set(gca, 'FontName', 'Arial');
    set(gca,'XColor',[0 0 0],'YColor',[0 0 0]);
    set(gca,'Box','on');
ylabel('Radius of Gyration (\mum)');
hFig=findall(0,'type','figure');
hLeg=findobj(hFig(1,1),'type','legend');
set(hLeg,'visible','off')

[WT_wthtt,p_WT_wthtt]=ttest2(rg_data{1}, rg_data{3})
[WT_S421D,p_WT_S421D]=ttest2(rg_data{1}, rg_data{2})
[S421D_WTS421Dhtt,p_S421D_WTS421Dhtt]=ttest2(rg_data{2}, rg_data{4})
[WThtt_WTS421Dhtt,p_WThtt_WTS421Dhtt]=ttest2(rg_data{3}, rg_data{4})
[WT_S421Dhtt,p_WT_S421Dhtt]=ttest2(rg_data{1}, rg_data{4})

[WT_wthtt_alpha,p_WT_wthtt_alpha]=ttest2(alph_f{1}, alph_f{3})
[WT_S421D_alpha,p_WT_S421D_alpha]=ttest2(alph_f{1}, alph_f{2})
[S421D_WTS421Dhtt_alpha,p_S421D_WTS421Dhtt_alpha]=ttest2(alph_f{2}, alph_f{4})
[WThtt_WTS421Dhtt_alpha,p_WThtt_WTS421Dhtt_alpha]=ttest2(alph_f{3}, alph_f{4})
[WT_S421Dhtt_alpha,p_WT_S421Dhtt_alpha]=ttest2(alph_f{1}, alph_f{4})
% 
% cd('/Volumes/Emily_htt/Huntingtin_Project/EE/Results/20200710_eeonly_fastexp/')
% saveas(figure(50),'20200715_Rg_per_cell_fastexp.png')
% 
% % MSD:
% figure, 
% cdfplot(msd_data{1});
% xlabel('MSD (\mum^2)'); 
% ylabel('CDF'); 
% publication_fig(0,0,1)
% 
% % MSD:
% figure(12), 
% cdfplot(msd_data{2});
% xlabel('MSD (\mum^2)'); 
% ylabel('CDF'); 
% publication_fig(0,0,1)
% 
% % MSD:
% figure(14), 
% cdfplot(msd_data{3});
% xlabel('MSD (\mum^2)'); 
% ylabel('CDF'); 
% publication_fig(0,0,1)
% 
% % MSD:
% figure(16), 
% cdfplot(msd_data{4});
% xlabel('MSD (\mum^2)'); 
% ylabel('CDF'); 
% publication_fig(0,0,1)

bootstrapped_all=[];

nsamp=min([numel(rg_all{1}),numel(rg_all{2}), numel(rg_all{3}), numel(rg_all{4})])

for dataset=1:4
bstrp_WT_lyso=Loic_bootstrap_code_04092019_em(rg_all{dataset},nsamp,1000,5);
bootstrapped_dataset{dataset}=bstrp_WT_lyso';
end


difference_WT_S421D = bootstrapped_dataset{1}.bstrap_means - bootstrapped_dataset{2}.bstrap_means;
difference_fake = bootstrapped_dataset{1}.bstrap_means - bootstrapped_dataset{1}.bstrap_means; %testing the bootstrapping algorithm
difference_WT_WThtt = bootstrapped_dataset{1}.bstrap_means - bootstrapped_dataset{3}.bstrap_means;
difference_WThtt_S421Dhtt = bootstrapped_dataset{3}.bstrap_means - bootstrapped_dataset{4}.bstrap_means;
difference_S421D_S421Dhtt = bootstrapped_dataset{2}.bstrap_means - bootstrapped_dataset{4}.bstrap_means;

figure('Name','Rg_bootstrapping_histogram','NumberTitle','off'), hold on,
histogram(difference_WT_S421D,'facecolor',colour_WT, 'edgecolor','none','facealpha',0.5), hold on,
histogram(difference_WT_WThtt, 'facecolor',colour_S421D, 'edgecolor','none','facealpha',0.2), hold on,
histogram(difference_WThtt_S421Dhtt, 'facecolor',colour_WThtt, 'edgecolor','none','facealpha',0.2), hold on,
histogram(difference_S421D_S421Dhtt, 'facecolor',colour_S421Dhtt, 'edgecolor','none','facealpha',0.2), hold on,
xlabel('Radius of Gyration (\mum)');
ylabel('Number of Trajectories');
publication_fig(0,0,1);

% figure('Name','Rg_bootstrapping_histogram_fake','NumberTitle','off'), hold on,
% histogram(difference_fake,'facecolor',[0 0 1], 'edgecolor','none','facealpha',0.5), hold on,
% xlabel('Radius of Gyration (\mum)');
% ylabel('Number of Trajectories');
% publication_fig(0,0,1);

figure('Name','Rg_bootstrapping_histogram_WT_S421D','NumberTitle','off'), hold on,
histogram(difference_WT_S421D,'facecolor',colour_WT, 'edgecolor','none','facealpha',0.5), hold on,
xlabel('Radius of Gyration (\mum)');
ylabel('Number of Trajectories');
publication_fig(0,0,1);

figure('Name','Rg_bootstrapping_histogram_WT_WThtt','NumberTitle','off'), hold on,
histogram(difference_WT_WThtt, 'facecolor',colour_WThtt, 'edgecolor','none','facealpha',0.2), hold on,
xlabel('Radius of Gyration (\mum)');
ylabel('Number of Trajectories');
publication_fig(0,0,1);

figure('Name','Rg_bootstrapping_histogram_WThtt_S421Dhtt','NumberTitle','off'), hold on,
histogram(difference_WThtt_S421Dhtt, 'facecolor',[0 0.75 0.75], 'edgecolor','none','facealpha',0.2), hold on,
xlabel('Radius of Gyration (\mum)');
ylabel('Number of Trajectories');
publication_fig(0,0,1);

figure('Name','Rg_bootstrapping_histogram_S421D_S421Dhtt','NumberTitle','off'), hold on,
histogram(difference_S421D_S421Dhtt, 'facecolor',colour_S421Dhtt, 'edgecolor','none','facealpha',0.5), hold on,
xlabel('Radius of Gyration (\mum)');
ylabel('Number of Trajectories');
publication_fig(0,0,1);

% nsamp=min([numel(alph_f{1}),numel(alph_f{2}), numel(alph_f{3}), numel(alph_f{4})])
% 
% for dataset=1:4
% bstrp_WT_lyso_alpha=Loic_bootstrap_code_04092019_em(alph_f{dataset},nsamp,1000,5);
% bootstrapped_dataset_alpha{dataset}=bstrp_WT_lyso_alpha';
% end
% 
% 
% difference_WT_S421D_alpha = bootstrapped_dataset_alpha{1}.bstrap_means - bootstrapped_dataset_alpha{2}.bstrap_means;
% difference_fake_alpha = bootstrapped_dataset_alpha{1}.bstrap_means - bootstrapped_dataset_alpha{1}.bstrap_means; %testing the bootstrapping
% difference_WT_WThtt_alpha = bootstrapped_dataset_alpha{1}.bstrap_means - bootstrapped_dataset_alpha{3}.bstrap_means;
% difference_WThtt_S421Dhtt_alpha = bootstrapped_dataset_alpha{3}.bstrap_means - bootstrapped_dataset_alpha{4}.bstrap_means;
% difference_S421D_S421Dhtt_alpha = bootstrapped_dataset_alpha{2}.bstrap_means - bootstrapped_dataset_alpha{4}.bstrap_means;
% 
% figure('Name','Bootstrapping_histogram_alpha','NumberTitle','off'), hold on,
% histogram(difference_WT_S421D_alpha,'facecolor',colour_WT, 'edgecolor','none','facealpha',0.5), hold on,
% histogram(difference_WT_WThtt_alpha, 'facecolor',[0 0 0], 'edgecolor','none','facealpha',0.2), hold on,
% histogram(difference_WThtt_S421Dhtt_alpha, 'facecolor',[0 0.75 0.75], 'edgecolor','none','facealpha',0.2), hold on,
% histogram(difference_S421D_S421Dhtt_alpha, 'facecolor',colour_S421Dhtt, 'edgecolor','none','facealpha',0.2), hold on,
% xlabel('\alpha');
% ylabel('Number of Trajectories');
% publication_fig(0,0,1);
% 
% figure('Name','Bootstrapping_fake_alpha','NumberTitle','off'), hold on,
% histogram(difference_fake_alpha,'facecolor',[0 0 1], 'edgecolor','none','facealpha',0.5), hold on,
% xlabel('\alpha');
% ylabel('Number of Trajectories');
% xlim([-1 1]);
% ylim([0 2000]);
% publication_fig(0,0,1);
% 
% 
% figure('Name','Bootstrapping_histogram_WT_S421D_alpha','NumberTitle','off'), hold on,
% histogram(difference_WT_S421D_alpha,'facecolor',colour_WT, 'edgecolor','none','facealpha',0.5), hold on,
% xlabel('\alpha');
% ylabel('Number of Trajectories');
% publication_fig(0,0,1);
% 
% figure('Name','Bootstrapping_histogram_WT_WThtt_alpha','NumberTitle','off'), hold on,
% histogram(difference_WT_WThtt_alpha, 'facecolor',[0 0 0], 'edgecolor','none','facealpha',0.2), hold on,
% xlabel('\alpha');
% ylabel('Number of Trajectories');
% publication_fig(0,0,1);
% 
% figure('Name','Bootstrapping_histogram_WThtt_S421Dhtt_alpha','NumberTitle','off'), hold on,
% histogram(difference_WThtt_S421Dhtt_alpha, 'facecolor',[0 0.75 0.75], 'edgecolor','none','facealpha',0.2), hold on,
% xlabel('\alpha');
% ylabel('Number of Trajectories');
% publication_fig(0,0,1);
% 
% figure('Name','Bootstrapping_histogram_S421D_S421Dhtt_alpha','NumberTitle','off'), hold on,
% histogram(difference_S421D_S421Dhtt_alpha, 'facecolor',colour_S421Dhtt, 'edgecolor','none','facealpha',0.5), hold on,
% xlabel('\alpha');
% ylabel('Number of Trajectories');
% publication_fig(0,0,1);
% 
% 
% 
figure('Name','Directional bias run length box and violin','NumberTitle','off'), hold on,
violin(plus_processive_runs, 'xlabel',{'WT','S421D','WT+WThtt','WT+S421Dhtt'},'facecolor',[colour_WT;colour_S421D;colour_WThtt;colour_S421Dhtt],'facealpha',0.3,'edgecolor','k','mc','k','medc','k--');
h1=notBoxPlot_nodotorline(plus_processive_runs{1},1,'style','line');
set(h1.data,'MarkerFaceColor','none','MarkerEdgeColor',colour_WT, 'MarkerSize', 5);
h2=notBoxPlot_nodotorline(plus_processive_runs{2},2,'style','line');
set(h2.data,'MarkerFaceColor','none','MarkerEdgeColor',colour_S421D,'MarkerSize', 5);
h3=notBoxPlot_nodotorline(plus_processive_runs{3},3,'style','line');
set(h3.data,'MarkerFaceColor','none','MarkerEdgeColor',[0.3 0.3 0.3],'MarkerSize', 5);
h4=notBoxPlot_nodotorline(plus_processive_runs{4},4,'style','line');
set(h4.data,'MarkerFaceColor','none','MarkerEdgeColor',colour_S421Dhtt,'MarkerSize', 5);

publication_fig(0,0,1)
axisHandle = gca; 
    set(findall(gca, 'Type', 'Line'),'LineWidth',2);
    set(gca,'LineWidth',2);
    set(gca,'FontSize',24);
    set(gca, 'FontName', 'Arial');
    set(gca,'XColor',[0 0 0],'YColor',[0 0 0]);
    set(gca,'Box','on');
ylabel('Fraction of processive runs moving outward');
 ylim([0.2 0.75]);
hFig=findall(0,'type','figure');
hLeg=findobj(hFig(1,1),'type','legend');
set(hLeg,'visible','off')

[WT_wthtt_dirbias,p_WT_wthtt_dirbias]=ttest2(plus_processive_runs{1}, plus_processive_runs{3})
[WT_S421D_dirbias,p_WT_S421D_dirbias]=ttest2(plus_processive_runs{1}, plus_processive_runs{2})
[S421D_WTS421Dhtt_dirbias,p_S421D_WTS421Dhtt_dirbias]=ttest2(plus_processive_runs{2}, plus_processive_runs{4})
[WThtt_WTS421Dhtt_dirbias,p_WThtt_WTS421Dhtt_dirbias]=ttest2(plus_processive_runs{3}, plus_processive_runs{4})
[WT_S421Dhtt_dirbias,p_WT_S421Dhtt_dirbias]=ttest2(plus_processive_runs{1}, plus_processive_runs{4})

% fileprefix='20210727_htttrans_r5_transhtt_per_cell_a_'
% % 
% tempdir = '/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/Figures_motility/';   % Your destination folder
% FolderName = tempdir;   % Your destination folder
% FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
% for iFig = 1:length(FigList)
%   FigHandle = FigList(iFig);
%   FigName   = get(FigHandle, 'Name');
%   savefig(FigHandle, fullfile(FolderName, [append(fileprefix,FigName), '.fig'])); 
%   saveas(FigHandle, fullfile(FolderName, [append(fileprefix,FigName), '.png']));
% end

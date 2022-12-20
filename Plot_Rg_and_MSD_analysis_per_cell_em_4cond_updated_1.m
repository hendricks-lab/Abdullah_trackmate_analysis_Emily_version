%% Analyze and Plot Rg and MSD per cell
%% Abdullah R. Chaudhary

close all; clear all; clc; 
set(0,'DefaultFigureWindowStyle','docked')
addpath('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/New_General_codes/');
addpath('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/New_General_codes/plotSpread/');
addpath('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/New_General_codes/raacampbell-notBoxPlot-7d90c27/code/');
addpath('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis//New_General_codes/raacampbell-notBoxPlot-7d90c27/code/+NBP/');
addpath('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/violin/');
addpath('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/Statistical_testing/');
addpath('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/');
addpath('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/AC_codes_epmodified_20200603/')

plus_processive_runs=[];
% %Rab5 colours
% colour_WT=[0.5 0.5 0.5];
% colour_S421D=[0 0.4 0.4];
% colour_WThtt=[0.3 0.3 0.3];
% colour_S421Dhtt=[0 0.4 0.68];

% %Lyso colours
% colour_WT=[0.729 0.831 0.957];
% colour_S421D=[0.757 0.867 0.776];
% colour_WThtt=[0.153 0.227 0.373];
% colour_S421Dhtt=[0.071 0.212 0.141];

%Rab5 Lat A colours
% colour_WT=[0.5 0.5 0.5];
% colour_S421D=[0 0.4 0.4];
% colour_WThtt=[0.3 0.3 0.3]; % Actually WT+latA
% colour_S421Dhtt=[0 0.2 0.2]; % Actually S421D+latA

%Quantum dot LatA colours
colour_WT=[0.8 0.8 0.8];
colour_S421D=[0 0.75 0.75];
colour_WThtt=[0.3 0.3 0.3];  % Actually WT+latA
colour_S421Dhtt=[0 0.2 0.2];% Actually S421D+latA

%R5 vs qd
% colour_WT=[0.5 0.5 0.5];
% colour_S421D=[0 0.4 0.4];
% colour_WThtt=[0.8 0.8 0.8]; %actually WT quantum dots
% colour_S421Dhtt=[0 0.75 0.75]; %actually S421D quantum dots

%Rab5 NZ colours
% colour_WT=[0.5 0.5 0.5];
% colour_S421D=[0.74 0 0.74];

%Lyso NZ colours
% colour_WT=[0.5 0.5 0.5];
% colour_S421D=[0.9 0.2 0.5];

for k_choose = 1:4

if k_choose == 1    % WT
    col1=colour_WT;
%     %Rab5 WT
%     pth1='/Volumes/Emily_htt/Huntingtin_Project/EE/WT/Rab5/WT_r5_rg/';
%     pth2='/Volumes/Emily_htt/Huntingtin_Project/EE/WT/Rab5/WT_r5_msd/';
%     pth3='/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/Dirbias/WT_r5_updated/';
%     %Lyso WT
%     pth1='/Volumes/Emily_htt/Huntingtin_Project/EE/WT/Lysotracker/lyso_rg/';
%     pth2='/Volumes/Emily_htt/Huntingtin_Project/EE/WT/Lysotracker/lyso_msd/';
%     pth3='/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/Dirbias/WT_lyso_updated/';
    %Qdots WT
    pth1='/Volumes/Emily_htt/Huntingtin_Project/EE/WT/EE_Qdot/WT_all_ee_rg_tm2/';
    pth2='/Volumes/Emily_htt/Huntingtin_Project/EE/WT/EE_Qdot/WT_all_ee_msd_tm2/';
    pth3='/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/Dirbias/WT_qd/';
    %Nocodazole control Rab5(WT)
%     pth1='/Volumes/Emily_htt_2/S421D_htt_2/Rab5/20220205_r5+NZ/20220205_r5+NZ_Tiffs/20220205_r5_control_rg/';
%     pth2='/Volumes/Emily_htt_2/S421D_htt_2/Rab5/20220205_r5+NZ/20220205_r5+NZ_Tiffs/20220205_r5_control_msd/';
%     pth3='/Volumes/Emily_htt_2/S421D_htt_2/Rab5/20220205_r5+NZ/20220205_r5+NZ_Tiffs/20220205_r5_control_pos/';
    % Nocodazole control Lyso (WT)
%     pth1='/Volumes/Emily_htt_2/S421D_htt_2/Lyso/Nocodazole/20220122_lyso_NZ/20220122_Tiffs/20220122_control_rg/';
%     pth2='/Volumes/Emily_htt_2/S421D_htt_2/Lyso/Nocodazole/20220122_lyso_NZ/20220122_Tiffs/20220122_control_msd/';
%     pth3='/Volumes/Emily_htt_2/S421D_htt_2/Lyso/Nocodazole/20220122_lyso_NZ/20220122_Tiffs/20220122_control_pos/';
    addpath(pth1);
    addpath(pth2);
    addpath(pth3);
    fl1='Rg_per_cell.mat';
    fl2='MSD_per_cell.mat';
    fl3='frac_proc_out_per_cell.mat';
    fl4='frac_proc_in_per_cell.mat';
    fl5='frac_diff_per_cell.mat';
    fl6='frac_stat_per_cell.mat';
    fl7='frac_diff_out_per_cell.mat';
    fl8='frac_diff_in_per_cell.mat';
    fl9='frac_proc_per_cell.mat';
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
elseif k_choose == 2    % S421D
    col1=colour_S421D;
    %Rab5 S421D
%     pth1='/Volumes/Emily_htt/Huntingtin_Project/EE/S421D/Rab5/S421D_r5_rg/';
%     pth2='/Volumes/Emily_htt/Huntingtin_Project/EE/S421D/Rab5/S421D_r5_msd/';
%     pth3='/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/Dirbias/S421D_r5_updated/';
    %Lyso S421D
%     pth1='/Volumes/Emily_htt/Huntingtin_Project/EE/S421D/Lysotracker/S421D_lyso_rg/';
%     pth2='/Volumes/Emily_htt/Huntingtin_Project/EE/S421D/Lysotracker/S421D_lyso_msd/';
%     pth3='/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/Dirbias/S421D_lyso_updated/';
%     %Qdots S421D
    pth1='/Volumes/Emily_htt/Huntingtin_Project/EE/S421D/Qdot_EE/S421D_all_ee_rg_tm2/';
    pth2='/Volumes/Emily_htt/Huntingtin_Project/EE/S421D/Qdot_EE/S421D_all_ee_msd_tm2/';
    pth3='/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/Dirbias/S421D_qd/';
    %Nocodazole treated Rab5(WT)
%     pth1='/Volumes/Emily_htt_2/S421D_htt_2/Rab5/20220205_r5+NZ/20220205_r5+NZ_Tiffs/20220205_r5_NZ_rg/';
%     pth2='/Volumes/Emily_htt_2/S421D_htt_2/Rab5/20220205_r5+NZ/20220205_r5+NZ_Tiffs/20220205_r5_NZ_msd/';
%     pth3='/Volumes/Emily_htt_2/S421D_htt_2/Rab5/20220205_r5+NZ/20220205_r5+NZ_Tiffs/20220205_r5_NZ_pos/';
    % Nocodazole treated Lyso (WT)
%     pth1='/Volumes/Emily_htt_2/S421D_htt_2/Lyso/Nocodazole/20220122_lyso_NZ/20220122_Tiffs/20220122_NZ_rg/';
%     pth2='/Volumes/Emily_htt_2/S421D_htt_2/Lyso/Nocodazole/20220122_lyso_NZ/20220122_Tiffs/20220122_NZ_msd/';
%     pth3='/Volumes/Emily_htt_2/S421D_htt_2/Lyso/Nocodazole/20220122_lyso_NZ/20220122_Tiffs/20220122_NZ_pos/';
    addpath(pth1);
    addpath(pth2);
    addpath(pth3);
    fl1='Rg_per_cell.mat';
    fl2='MSD_per_cell.mat';
    fl3='frac_proc_out_per_cell.mat';
    fl4='frac_proc_in_per_cell.mat';
    fl5='frac_diff_per_cell.mat';
    fl6='frac_stat_per_cell.mat';
    fl7='frac_diff_out_per_cell.mat';
    fl8='frac_diff_in_per_cell.mat';
    fl9='frac_proc_per_cell.mat';
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

elseif k_choose == 3   % WT +WT htt
    col1=colour_WThtt;
    %Rab5 WThtt
%     pth1='/Volumes/Emily_htt/Huntingtin_Project/EE/WT_trans_mchWThtt/Rab5/WT_trans_WThtt_rg/';
%     pth2='/Volumes/Emily_htt/Huntingtin_Project/EE/WT_trans_mchWThtt/Rab5/WT_trans_WThtt_msd/';
%     pth3='/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/Dirbias/WT+WThtt_r5_updated/';
%     %Lyso WThtt
%     pth1='/Volumes/Emily_htt/Huntingtin_Project/EE/WT_trans_mchWThtt/Lyso/WT_trans_WThtt_lyso_rg/';
%     pth2='/Volumes/Emily_htt/Huntingtin_Project/EE/WT_trans_mchWThtt/Lyso/WT_trans_WThtt_lyso_msd/';
%     pth3='/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/Dirbias/WT+WThtt_lyso_updated/';
    %Rab5 LatA WT
%     pth1='/Volumes/Emily_htt/Huntingtin_Project/EE/WT/Rab5_Latrunculin/WT_r5_latA_rg/';
%     pth2='/Volumes/Emily_htt/Huntingtin_Project/EE/WT/Rab5_Latrunculin/WT_r5_latA_msd/';
%     pth3='/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/Dirbias/WT_r5_latA/';
    %Qdots WT
%     pth1='/Volumes/Emily_htt/Huntingtin_Project/EE/WT/EE_Qdot/WT_all_ee_rg_tm2/';
%     pth2='/Volumes/Emily_htt/Huntingtin_Project/EE/WT/EE_Qdot/WT_all_ee_msd_tm2/';
%     pth3='/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/Dirbias/WT_qd/';
%     Qdots WT+LatA
    pth1='/Volumes/Emily_htt/Huntingtin_Project/EE/WT/Latrunculin_Qdot/WT_latA_rg/';
    pth2='/Volumes/Emily_htt/Huntingtin_Project/EE/WT/Latrunculin_Qdot/WT_latA_msd/';
    pth3='/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/Dirbias/WT_qd_latA/';
    addpath(pth1);
    addpath(pth2);
    addpath(pth3);
    fl1='Rg_per_cell.mat';
    fl2='MSD_per_cell.mat';
    fl3='frac_proc_out_per_cell.mat';
    fl4='frac_proc_in_per_cell.mat';
    fl5='frac_diff_per_cell.mat';
    fl6='frac_stat_per_cell.mat';
    fl7='frac_diff_out_per_cell.mat';
    fl8='frac_diff_in_per_cell.mat';
    fl9='frac_proc_per_cell.mat';
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
    
elseif k_choose == 4   % WT+S421Dhtt
    col1=colour_S421Dhtt;
%     %Rab5 S421Dhtt
%     pth1='/Volumes/Emily_htt/Huntingtin_Project/EE/WT_trans_S421Dhtt/Rab5/WT_trans_S421Dhtt_rg/';
%     pth2='/Volumes/Emily_htt/Huntingtin_Project/EE/WT_trans_S421Dhtt/Rab5/WT_trans_S421Dhtt_msd/';
%     pth3='/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/Dirbias/WT+S421Dhtt_r5_updated/';
%     %Lyso S421Dhtt
%     pth1='/Volumes/Emily_htt/Huntingtin_Project/EE/WT_trans_S421Dhtt/Lyso/WT_trans_S421Dhtt_lyso_rg/';
%     pth2='/Volumes/Emily_htt/Huntingtin_Project/EE/WT_trans_S421Dhtt/Lyso/WT_trans_S421Dhtt_lyso_msd/';
%     pth3='/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/Dirbias/WT+S421Dhtt_lyso_updated/';
    %Rab5 lat A S421D
%     pth1='/Volumes/Emily_htt/Huntingtin_Project/EE/S421D/Rab5_Latrunculin/S421D_r5_latA_rg/';
%     pth2='/Volumes/Emily_htt/Huntingtin_Project/EE/S421D/Rab5_Latrunculin/S421D_r5_latA_msd/';
%     pth3='/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/Dirbias/S421D_r5_latA/';
    %Qdots S421D
%     pth1='/Volumes/Emily_htt/Huntingtin_Project/EE/S421D/Qdot_EE/S421D_all_ee_rg_tm2/';
%     pth2='/Volumes/Emily_htt/Huntingtin_Project/EE/S421D/Qdot_EE/S421D_all_ee_msd_tm2/';
%     pth3='/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/Dirbias/S421D_qd/';
%     %Qdots S421D+latA
    pth1='/Volumes/Emily_htt/Huntingtin_Project/EE/S421D/Qdot_Latrunculin/S421D_latA_rg/';
    pth2='/Volumes/Emily_htt/Huntingtin_Project/EE/S421D/Qdot_Latrunculin/S421D_latA_msd/';
    pth3='/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/Dirbias/S421D_qd_latA/';
%     
    addpath(pth1);
    addpath(pth2);
    addpath(pth3);
    fl1='Rg_per_cell.mat';
    fl2='MSD_per_cell.mat';
    fl3='frac_proc_out_per_cell.mat';
    fl4='frac_proc_in_per_cell.mat';
    fl5='frac_diff_per_cell.mat';
    fl6='frac_stat_per_cell.mat';
    fl7='frac_diff_out_per_cell.mat';
    fl8='frac_diff_in_per_cell.mat';
    fl9='frac_proc_per_cell.mat';
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
load(fullfile(pth3,fl3)); %Emily added
load(fullfile(pth3,fl4)); %Emily added 20220927
load(fullfile(pth3,fl5)); %Emily added 20220927
load(fullfile(pth3,fl6)); %Emily added 20220927
load(fullfile(pth3,fl7)); %Emily added 20221020
load(fullfile(pth3,fl8)); %Emily added 20221020
load(fullfile(pth3,fl9));

%Emily adding in a longer minimum run length version 20221024
proc_frac_t{k_choose}=t_proc_per_cell;
plus_proc_frac_t{k_choose}=t_proc_out_per_cell;
minus_proc_frac_t{k_choose}=t_proc_in_per_cell;
diff_frac_t{k_choose}=t_diff_per_cell;
plus_diff_frac_t{k_choose}=t_diff_out_per_cell; %Emily added 20221020
minus_diff_frac_t{k_choose}=t_diff_in_per_cell; %Emily added 20221020
stat_frac_t{k_choose}=t_stat_per_cell;

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

%Emily added for testing if min_lr was set correctly 20221020
figure('Name','Diff_runs_in_out_box','NumberTitle','off'), hold on,
% violin(plus_diff_frac_t, 'xlabel',{'WT','S421D','WT+WThtt','WT+S421Dhtt', 'WT', 'S421D','WT+WThtt','WT+S421Dhtt'},'facecolor',[colour_WT;colour_S421D;colour_WThtt;colour_S421Dhtt;colour_WT;colour_S421D;colour_WThtt;colour_S421Dhtt],'facealpha',0.3,'edgecolor','k','mc','k','medc','k--');
violin(plus_diff_frac_t, 'xlabel',{'WT','S421D','WT+WThtt','WT+S421Dhtt'},'facecolor',[colour_WT;colour_S421D;colour_WThtt;colour_S421Dhtt],'facealpha',0.3,'edgecolor','k','mc','k','medc','k--');
h5=notBoxPlot_nodotorline(plus_diff_frac_t{1},1,'style','line');
set(h5.data,'MarkerFaceColor','none','MarkerEdgeColor',colour_WT, 'MarkerSize', 5);
h6=notBoxPlot_nodotorline(plus_diff_frac_t{2},2,'style','line');
set(h6.data,'MarkerFaceColor','none','MarkerEdgeColor',colour_S421D,'MarkerSize', 5);
h7=notBoxPlot_nodotorline(plus_diff_frac_t{3},3,'style','line');
set(h7.data,'MarkerFaceColor','none','MarkerEdgeColor',[0.3 0.3 0.3],'MarkerSize', 5);
h8=notBoxPlot_nodotorline(plus_diff_frac_t{4},4,'style','line');
set(h8.data,'MarkerFaceColor','none','MarkerEdgeColor',colour_S421Dhtt,'MarkerSize', 5);
ylabel('Fraction of time of diffusive runs moving outwards');
publication_fig(0,0,1);
hFig=findall(0,'type','figure');
hLeg=findobj(hFig(1,1),'type','legend');
set(hLeg,'visible','off');

figure('Name','Directional bias run length box','NumberTitle','off'), hold on,
violin(plus_proc_frac_t, 'xlabel',{'WT','S421D','WT+WThtt','WT+S421Dhtt'},'facecolor',[colour_WT;colour_S421D;colour_WThtt;colour_S421Dhtt],'facealpha',0.3,'edgecolor','k','mc','k','medc','k--');
% violin(plus_frac_proc, 'xlabel',{'WT','S421D','WT+latA','WT+latA'},'facecolor',[colour_WT;colour_S421D;colour_WThtt;colour_S421Dhtt],'facealpha',0.3,'edgecolor','k','mc','k','medc','k--');
h1=notBoxPlot_nodotorline(plus_proc_frac_t{1},1,'style','line');
set(h1.data,'MarkerFaceColor','none','MarkerEdgeColor',colour_WT, 'MarkerSize', 5);
h2=notBoxPlot_nodotorline(plus_proc_frac_t{2},2,'style','line');
set(h2.data,'MarkerFaceColor','none','MarkerEdgeColor',colour_S421D,'MarkerSize', 5);
h3=notBoxPlot_nodotorline(plus_proc_frac_t{3},3,'style','line');
set(h3.data,'MarkerFaceColor','none','MarkerEdgeColor',[0.3 0.3 0.3],'MarkerSize', 5);
h4=notBoxPlot_nodotorline(plus_proc_frac_t{4},4,'style','line');
set(h4.data,'MarkerFaceColor','none','MarkerEdgeColor',colour_S421Dhtt,'MarkerSize', 5);
ylabel('Fraction of time of processive runs moving outwards');
publication_fig(0,0,1)
hFig=findall(0,'type','figure');
hLeg=findobj(hFig(1,1),'type','legend');
set(hLeg,'visible','off');

plus_frac_proc_WT_S421D_only{1}=plus_proc_frac_t{1};
plus_frac_proc_WT_S421D_only{2}=plus_proc_frac_t{2};

figure('Name','Directional bias run length box only WT and S421D','NumberTitle','off'), hold on,
violin(plus_frac_proc_WT_S421D_only, 'xlabel',{'WT','S421D'},'facecolor',[colour_WT;colour_S421D],'facealpha',0.3,'edgecolor','k','mc','k','medc','k--');
h1=notBoxPlot_nodotorline(plus_proc_frac_t{1},1,'style','line');
set(h1.data,'MarkerFaceColor','none','MarkerEdgeColor',colour_WT, 'MarkerSize', 5);
h2=notBoxPlot_nodotorline(plus_proc_frac_t{2},2,'style','line');
set(h2.data,'MarkerFaceColor','none','MarkerEdgeColor',colour_S421D,'MarkerSize', 5);
ylabel('Fraction of time of processive runs moving outwards');
publication_fig(0,0,1)
hFig=findall(0,'type','figure');
hLeg=findobj(hFig(1,1),'type','legend');
set(hLeg,'visible','off');

plus_runs_WT=mean(plus_proc_frac_t{1});
plus_runs_S421D=mean(plus_proc_frac_t{2});
plus_runs_WThtt=mean(plus_proc_frac_t{3});
plus_runs_S421Dhtt=mean(plus_proc_frac_t{4});
minus_runs_WT=mean(minus_proc_frac_t{1});
minus_runs_S421D=mean(minus_proc_frac_t{2});
minus_runs_WThtt=mean(minus_proc_frac_t{3});
minus_runs_S421Dhtt=mean(minus_proc_frac_t{4});
proc_runs_WT=mean(proc_frac_t{1});
proc_runs_S421D=mean(proc_frac_t{2});
proc_runs_WThtt=mean(proc_frac_t{3});
proc_runs_S421Dhtt=mean(proc_frac_t{4});
diff_frac_WT=mean(diff_frac_t{1});
diff_frac_S421D=mean(diff_frac_t{2});
diff_frac_WThtt=mean(diff_frac_t{3});
diff_frac_S421Dhtt=mean(diff_frac_t{4});
stat_frac_WT=mean(stat_frac_t{1});
stat_frac_S421D=mean(stat_frac_t{2});
stat_frac_WThtt=mean(stat_frac_t{3});
stat_frac_S421Dhtt=mean(stat_frac_t{4});

%SEM for directional bias parameters
error_proc_WT=std(proc_frac_t{1})/(sqrt(length(proc_frac_t{1})));
error_diff_WT=std(diff_frac_t{1})/(sqrt(length(diff_frac_t{1})));
error_stat_WT=std(stat_frac_t{1})/(sqrt(length(stat_frac_t{1})));
error_proc_S421D=std(proc_frac_t{2})/(sqrt(length(proc_frac_t{2})));
error_diff_S421D=std(diff_frac_t{2})/(sqrt(length(diff_frac_t{2})));
error_stat_S421D=std(stat_frac_t{2})/(sqrt(length(stat_frac_t{2})));
error_proc_WThtt=std(proc_frac_t{3})/(sqrt(length(proc_frac_t{3})));
error_diff_WThtt=std(diff_frac_t{3})/(sqrt(length(diff_frac_t{3})));
error_stat_WThtt=std(stat_frac_t{3})/(sqrt(length(stat_frac_t{3})));
error_proc_S421Dhtt=std(proc_frac_t{4})/(sqrt(length(proc_frac_t{4})));
error_diff_S421Dhtt=std(diff_frac_t{4})/(sqrt(length(diff_frac_t{4})));
error_stat_S421Dhtt=std(stat_frac_t{4})/(sqrt(length(stat_frac_t{4})));
%creating a variable with all the errors for all the conditions
all_err_bars=[error_proc_WT error_proc_S421D error_proc_WThtt error_proc_S421Dhtt; error_diff_WT error_diff_S421D error_diff_WThtt error_diff_S421Dhtt; error_stat_WT error_stat_S421D error_stat_WThtt error_stat_S421Dhtt];

values=[proc_runs_WT proc_runs_S421D proc_runs_WThtt proc_runs_S421Dhtt; diff_frac_WT diff_frac_S421D diff_frac_WThtt diff_frac_S421Dhtt;stat_frac_WT stat_frac_S421D stat_frac_WThtt stat_frac_S421Dhtt];
figure('Name','Bar plot proc stat diff','NumberTitle','off'), hold on,
b2=bar(values,'grouped');
b2(1).FaceColor=colour_WT;
b2(2).FaceColor=colour_S421D;
b2(3).FaceColor=colour_WThtt;
b2(4).FaceColor=colour_S421Dhtt;
b2(1).FaceAlpha=0.5;
b2(2).FaceAlpha=0.5;
b2(3).FaceAlpha=0.5;
b2(4).FaceAlpha=0.5;
hold on
[ngroups,nbars] = size(values);
groupwidth = min(0.8, nbars/(nbars + 1.5));
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
for i = 1:nbars
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, values(:,i), all_err_bars(:,i), 'k', 'linestyle', 'none','LineWidth',2);
end
xticks([1 2 3]);
xticklabels({'Processive', 'Diffusive','Stationary'});
ylabel('Fraction of runs');
publication_fig(0,0,1)

%only WT S421D bar plot proc/diff/stat
all_err_bars_2=all_err_bars(1:3,1:2);
values_2=[proc_runs_WT proc_runs_S421D; diff_frac_WT diff_frac_S421D; stat_frac_WT stat_frac_S421D;];
figure('Name','Bar plot proc stat diff only WT and S421D','NumberTitle','off'), hold on,
b2=bar(values_2, 'grouped');
b2(1).FaceColor=colour_WT;
b2(2).FaceColor=colour_S421D;
b2(1).FaceAlpha=0.5;
b2(2).FaceAlpha=0.5;
[ngroups,nbars] = size(values_2);
groupwidth = min(0.8, nbars/(nbars + 1.5));
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
for i = 1:nbars
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, values_2(:,i), all_err_bars_2(:,i), 'k', 'linestyle', 'none','LineWidth',2);
end
xticks([1 2 3]);
xticklabels({'Processive', 'Diffusive','Stationary'});
ylabel('Fraction of runs');
publication_fig(0,0,1)

%Bootstrapping Rg and Alpha values

% bootstrapped_all=[];
% 
% nsamp=min([numel(rg_all{1}),numel(rg_all{2}), numel(rg_all{3}), numel(rg_all{4})])
% 
% for dataset=1:4
% bstrp_WT_lyso=Loic_bootstrap_code_04092019_em(rg_all{dataset},nsamp,1000,5);
% bootstrapped_dataset{dataset}=bstrp_WT_lyso';
% end
% 
% 
% difference_WT_S421D = bootstrapped_dataset{1}.bstrap_means - bootstrapped_dataset{2}.bstrap_means;
% difference_fake = bootstrapped_dataset{1}.bstrap_means - bootstrapped_dataset{1}.bstrap_means; %testing the bootstrapping algorithm
% difference_WT_WThtt = bootstrapped_dataset{1}.bstrap_means - bootstrapped_dataset{3}.bstrap_means;
% difference_WThtt_S421Dhtt = bootstrapped_dataset{3}.bstrap_means - bootstrapped_dataset{4}.bstrap_means;
% difference_S421D_S421Dhtt = bootstrapped_dataset{2}.bstrap_means - bootstrapped_dataset{4}.bstrap_means;
% 
% figure('Name','Rg_bootstrapping_histogram','NumberTitle','off'), hold on,
% histogram(difference_WT_S421D,'facecolor',colour_WT, 'edgecolor','none','facealpha',0.5), hold on,
% histogram(difference_WT_WThtt, 'facecolor',colour_S421D, 'edgecolor','none','facealpha',0.2), hold on,
% histogram(difference_WThtt_S421Dhtt, 'facecolor',colour_WThtt, 'edgecolor','none','facealpha',0.2), hold on,
% histogram(difference_S421D_S421Dhtt, 'facecolor',colour_S421Dhtt, 'edgecolor','none','facealpha',0.2), hold on,
% xlabel('Radius of Gyration (\mum)');
% ylabel('Number of Trajectories');
% publication_fig(0,0,1);
% 
% % figure('Name','Rg_bootstrapping_histogram_fake','NumberTitle','off'), hold on,
% % histogram(difference_fake,'facecolor',[0 0 1], 'edgecolor','none','facealpha',0.5), hold on,
% % xlabel('Radius of Gyration (\mum)');
% % ylabel('Number of Trajectories');
% % publication_fig(0,0,1);
% 
% figure('Name','Rg_bootstrapping_histogram_WT_S421D','NumberTitle','off'), hold on,
% histogram(difference_WT_S421D,'facecolor',colour_WT, 'edgecolor','none','facealpha',0.5), hold on,
% xlabel('Radius of Gyration (\mum)');
% ylabel('Number of Trajectories');
% publication_fig(0,0,1);
% 
% figure('Name','Rg_bootstrapping_histogram_WT_WThtt','NumberTitle','off'), hold on,
% histogram(difference_WT_WThtt, 'facecolor',colour_WThtt, 'edgecolor','none','facealpha',0.2), hold on,
% xlabel('Radius of Gyration (\mum)');
% ylabel('Number of Trajectories');
% publication_fig(0,0,1);
% 
% figure('Name','Rg_bootstrapping_histogram_WThtt_S421Dhtt','NumberTitle','off'), hold on,
% histogram(difference_WThtt_S421Dhtt, 'facecolor',[0 0.75 0.75], 'edgecolor','none','facealpha',0.2), hold on,
% xlabel('Radius of Gyration (\mum)');
% ylabel('Number of Trajectories');
% publication_fig(0,0,1);
% 
% figure('Name','Rg_bootstrapping_histogram_S421D_S421Dhtt','NumberTitle','off'), hold on,
% histogram(difference_S421D_S421Dhtt, 'facecolor',colour_S421Dhtt, 'edgecolor','none','facealpha',0.5), hold on,
% xlabel('Radius of Gyration (\mum)');
% ylabel('Number of Trajectories');
% publication_fig(0,0,1);

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

%Statistical testing (T tests)
[WT_wthtt_rg,p_WT_wthtt_rg]=ttest2(rg_data{1}, rg_data{3});
[WT_S421D_rg,p_WT_S421D_rg]=ttest2(rg_data{1}, rg_data{2});
[S421D_WTS421Dhtt_rg,p_S421D_WTS421Dhtt_rg]=ttest2(rg_data{2}, rg_data{4});
[WThtt_WTS421Dhtt_rg,p_WThtt_WTS421Dhtt_rg]=ttest2(rg_data{3}, rg_data{4});
[WT_S421Dhtt_rg,p_WT_S421Dhtt_rg]=ttest2(rg_data{1}, rg_data{4});

rg_stats=[p_WT_S421D_rg;p_WT_wthtt_rg; p_WThtt_WTS421Dhtt_rg; p_S421D_WTS421Dhtt_rg; p_WT_S421Dhtt_rg];

[WT_wthtt_alpha,p_WT_wthtt_alpha]=ttest2(alph_f{1}, alph_f{3});
[WT_S421D_alpha,p_WT_S421D_alpha]=ttest2(alph_f{1}, alph_f{2});
[S421D_WTS421Dhtt_alpha,p_S421D_WTS421Dhtt_alpha]=ttest2(alph_f{2}, alph_f{4});
[WThtt_WTS421Dhtt_alpha,p_WThtt_WTS421Dhtt_alpha]=ttest2(alph_f{3}, alph_f{4});
[WT_S421Dhtt_alpha,p_WT_S421Dhtt_alpha]=ttest2(alph_f{1}, alph_f{4});

alpha_stats=[p_WT_S421D_alpha;p_WT_wthtt_alpha; p_WThtt_WTS421Dhtt_alpha; p_S421D_WTS421Dhtt_alpha;p_WT_S421Dhtt_alpha];
%Statistical testing for directional bias

[WT_wthtt_dirbias,p_WT_wthtt_dirbias]=ttest2(plus_proc_frac_t{1}, plus_proc_frac_t{3});
[WT_S421D_dirbias,p_WT_S421D_dirbias]=ttest2(plus_proc_frac_t{1}, plus_proc_frac_t{2});
[S421D_WTS421Dhtt_dirbias,p_S421D_WTS421Dhtt_dirbias]=ttest2(plus_proc_frac_t{2}, plus_proc_frac_t{4});
[WThtt_WTS421Dhtt_dirbias,p_WThtt_WTS421Dhtt_dirbias]=ttest2(plus_proc_frac_t{3}, plus_proc_frac_t{4});
[WT_S421Dhtt_dirbias,p_WT_S421Dhtt_dirbias]=ttest2(plus_proc_frac_t{1}, plus_proc_frac_t{4});
 dirbias_stats=[p_WT_S421D_dirbias;p_WT_wthtt_dirbias; p_WThtt_WTS421Dhtt_dirbias; p_S421D_WTS421Dhtt_dirbias;p_WT_S421Dhtt_dirbias];

[WT_S421D_proc,p_WT_S421D_proc]=ttest2(proc_frac_t{1}, proc_frac_t{2});
[WT_S421D_diff,p_WT_S421D_diff]=ttest2(diff_frac_t{1}, diff_frac_t{2});
[WT_S421D_stat,p_WT_S421D_stat]=ttest2(stat_frac_t{1}, stat_frac_t{2});

[WT_S421Dhtt_proc,p_WT_S421Dhtt_proc]=ttest2(proc_frac_t{1}, proc_frac_t{4});
[WT_S421Dhtt_diff,p_WT_S421Dhtt_diff]=ttest2(diff_frac_t{1}, diff_frac_t{4});
[WT_S421Dhtt_stat,p_WT_S421Dhtt_stat]=ttest2(stat_frac_t{1}, stat_frac_t{4});

[WT_WThtt_proc,p_WT_WThtt_proc]=ttest2(proc_frac_t{1}, proc_frac_t{3});
[WT_WThtt_diff,p_WT_WThtt_diff]=ttest2(diff_frac_t{1}, diff_frac_t{3});
[WT_WThtt_stat,p_WT_WThtt_stat]=ttest2(stat_frac_t{1}, stat_frac_t{3});

[S421D_S421Dhtt_proc,p_S421D_S421Dhtt_proc]=ttest2(proc_frac_t{2}, proc_frac_t{4});
[S421D_S421Dhtt_diff,p_S421D_S421Dhtt_diff]=ttest2(diff_frac_t{2}, diff_frac_t{4});
[S421D_S421Dhtt_stat,p_S421D_S421Dhtt_stat]=ttest2(stat_frac_t{2}, stat_frac_t{4});

[WThtt_S421Dhtt_proc,p_WThtt_S421Dhtt_proc]=ttest2(proc_frac_t{3}, proc_frac_t{4});
[WThtt_S421Dhtt_diff,p_WThtt_S421Dhtt_diff]=ttest2(diff_frac_t{3}, diff_frac_t{4});
[WThtt_S421Dhtt_stat,p_WThtt_S421Dhtt_stat]=ttest2(stat_frac_t{3}, stat_frac_t{4});

proc_stats=[p_WT_S421D_proc;p_WT_WThtt_proc; p_WThtt_S421Dhtt_proc; p_S421D_S421Dhtt_proc;p_WT_S421Dhtt_proc];
diff_stats=[p_WT_S421D_diff;p_WT_WThtt_diff; p_WThtt_S421Dhtt_diff; p_S421D_S421Dhtt_diff;p_WT_S421Dhtt_diff];
stat_stats=[p_WT_S421D_stat;p_WT_WThtt_stat; p_WThtt_S421Dhtt_stat; p_S421D_S421Dhtt_stat;p_WT_S421Dhtt_stat];

titles={'WT vs S421D'; 'WT vs WThtt'; 'WThtt vs S421Dhtt'; 'S421D vs S421Dhtt'; 'WT vs S421Dhtt'};

%Auto-saving figures

% fileprefix='20221122_qdot_latA_min_lr_0pt4_transhtt_per_cell';
% tempdir_1 = '/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/Motility_Stats/';   % Your destination folder
% FolderName_1 = tempdir_1;   % Your destination folder
% 
% stats_table=table(titles, rg_stats, alpha_stats,dirbias_stats, proc_stats, diff_stats, stat_stats);
% writetable(stats_table,fullfile(FolderName_1, [fileprefix,'.csv']),'WriteRowNames',true);
% 
% tempdir = '/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/Figures_motility/';   % Your destination folder
% FolderName = tempdir;   % Your destination folder
% FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
% for iFig = 1:length(FigList)
%   FigHandle = FigList(iFig);
%   FigName   = get(FigHandle, 'Name');
%   savefig(FigHandle, fullfile(FolderName, [append(fileprefix,FigName), '.fig'])); 
%   saveas(FigHandle, fullfile(FolderName, [append(fileprefix,FigName), '.png']));
% end

%% Bootstrap MSD and posititional probability Analysis: 
%% Abdullah R. Chaudhary

close all; clear all; clc; 
set(0,'DefaultFigureWindowStyle','docked')
addpath('/Volumes/Emily_htt/New_General_codes/');
addpath('/Volumes/Emily_htt/AC_codes_epmodified_20200603/');
addpath('/Volumes/Emily_htt/New_General_codes/plotSpread/');
addpath('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/violin/');
addpath('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/');
addpath('/Volumes/Emily_htt/New_General_codes/raacampbell-notBoxPlot-7d90c27/code/');
addpath('/Volumes/Emily_htt/New_General_codes/raacampbell-notBoxPlot-7d90c27/code/+NBP/');

%% Savitsky Golay smoothing
span=10;    % 10 is good, 25 is good as well
pwr=1;      % 1 is good, 2 is good


plus_proc_times=[];
average_fraction_of_outward_runs_per_cell=[];

%% Reset these parameters everytime you run the code !!!!!!
%%%%%%%%%%%%%%%%%%%%%%%%%%
alpha_MSD=[];
ab1=[];
kg=0;
variance_plus_traj=[];
variance_minus_traj=[];

%% Savitsky Golay smoothing
span=10;    % 10 is good, 25 is good as well
pwr=1;      % 1 is good, 2 is good
%%%%%%%%%%%%%%%%%%%%%%%%%%

tic
for k_choose = 2 %need to run each k_choose separately in this code in order to get the results to save properly
    
if k_choose == 1    % WT
    col1=[0.5 0.5 0.5];
    cd('/Volumes/Emily_htt/Huntingtin_Project/EE/WT/EE_Qdot/Processed_Data/WT_all_ee_mats/');
    save_dir='/Volumes/Emily_htt/Huntingtin_Project/EE/WT/EE_Qdot/Processed_Data/WT_all_ee_pos/';
    save_dir_Rg='/Volumes/Emily_htt/Huntingtin_Project/EE/WT/EE_Qdot/Processed_Data/WT_all_ee_rg/';
    save_dir_MSD='/Volumes/Emily_htt/Huntingtin_Project/EE/WT/EE_Qdot/Processed_Data/WT_all_ee_msd/';
    save_dir_dirbias='/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/Dirbias/WT_qd/';
    fl='*.mat';
    an=0;
    DT=0.12;
    kf=0;
    proc_runs=[];
    proc_times=[];
    rng=[];
    tlp=-0.5;
    fraction_plus_proc_time=[]; %Emily added
    proc_run_em=[]; %Emily added
    proc_time_em=[]; %Emily added
    
elseif k_choose == 2  % S421D
    col1='b';
    cd('/Volumes/Emily_htt/Huntingtin_Project/EE/S421D/Qdot_EE/S421D_all_ee_mats/');
    save_dir='/Volumes/Emily_htt/Huntingtin_Project/EE/S421D/Qdot_EE/S421D_all_ee_pos/';
    save_dir_Rg='/Volumes/Emily_htt/Huntingtin_Project/EE/S421D/Qdot_EE/S421D_all_ee_rg/';
    save_dir_MSD='/Volumes/Emily_htt/Huntingtin_Project/EE/S421D/Qdot_EE/S421D_all_ee_msd/';
    save_dir_dirbias='/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/Dirbias/S421D_qd/';
    fl='*.mat';
    an=0;
    DT=0.12;
    kf=0;
    proc_runs=[];
    proc_times=[];
    rng=[];
    tlp=0.5;
    fraction_plus_proc_time=[]; %Emily added
    proc_run_em=[]; %Emily added
    proc_time_em=[]; %Emily added

end

fls=dir(fullfile(dat_dir,fl));
pos=[];
min_lr=0.4; %minimum run length is 0.4um, Emily added
time=DT*1379; %Emily added
    
for k=1:numel(fls) %k is the current file
    %pos=[];
    kg=kg+1;%still not sure why we need both kg and k, kg is the iteration number of the loop
    load(fullfile(dat_dir,fls(k).name)); %load filename
    ab = ab(~cellfun(@isempty,{ab.position})); % outputs a logical array of if the position is empty (0 for empty) or not.
    ab1{kg}=ab; %indexing by the iteration number of the loop to get the position array for that file
    rg_1=[];
    for j=1:numel(ab) %from 1 to the number of elements in the position file, which would be the number of trajectories for a file
        if numel(ab(j).position)>15 %only analyze trajectories that have more than 15 timepoints of positional data
        kf=kf+1; %again not sure why we need this and j, now I sortof do because j is reset to zero at the end of the loop for the file while kf is not.
        xk{kf}=smooth(ab(j).xk,span,'sgolay',pwr); %smoothing the position function so we don't get too many sharp switches in direction which are artificial
        yk{kf}=smooth(ab(j).yk,span,'sgolay',pwr); %smoothing the position function so we don't get too many sharp switches in direction which are artificial
        %line 80 is calculating the 2D MSD, which is used for calculating
        %the per cell alpha value I'm pretty sure. Check the function if
        %you want to know how it is calculating the MSD.
        [deltat, msdpts, sem, log_deltat, log_msdpts, alpha_1, DiffCoef] = MSD_2d_emmessed(xk{kf},yk{kf}, DT, k_choose);
        
        % Emily added this next section to determine the directional bias of trajectories per
        % cell
        curr_pos=[xk{kf}, yk{kf}];%Emily added for current position datapoints of the trajectory
        time=[1:numel(curr_pos)]'.*DT;
        res=analyze_run_length_reversals_v3_per_cell_bias(time,curr_pos,min_lr,0,k_choose);
    
        proc_run_em=[proc_run_em;res.proc_run_em];
        proc_time_em=[proc_time_em;res.proc_time_em];
    
        fraction_plus_proc_time=[fraction_plus_proc_time;(sum(res.proc_time_em(find(res.proc_run_em>0)))/sum(res.proc_time_em))]; %Emily added to find plus end runs
         
      
        if alpha_1>0 %not analyzing the tracks that are so flat they look
%         like they have negative slope I guess
        rg_0=func_Rg_Linda_v2(xk{kf},yk{kf}); 
        % (above) calculating the radius of gyration (rg) from the smoothed
        % position data of the trajectory
        rg_1=[rg_1,rg_0]; %adds the currently calculated rg to the previous one in a list. 
        dat_pos_tim{kf}=smooth(ab(j).position,10,'sgolay',1); %gives another smoothed position file for the current trajectory?
        else %if alpha is less than or equal to zero, this condition will run
%         rg_1=[]; %commented out because the value of rg_1 was resetting
%         for the same cell multiple times, providing lists of length 2 for
%         each cell which caused poor mle fitting.
        dat_pos_tim{kf}=[];%this was still required, not entirely sure why.
        end
        
        else %do not analyze trajectories with <15 timepoints of positions.
        end
      
    end %the end of the function analyzing trajectories
    pos=[pos,dat_pos_tim];
    position{k}=pos;
    
    save([save_dir, num2str(k_choose) 'Position_per_cell'],'position');
    
    % Clear 1D MSD
    K=1:40;
    Ndat=numel(xk);
    %below is a calculation of the 1D MSD
    res=msd_fun_v1_emmessedwithit(K,pos,DT,Ndat,0);% change the value to 0 to stop plotting
    pos=[]; %clearing the values stored in pos
    res2{k}=res; %saving this 1D MSD into res2, indexed by file number (cell)
    rg_all1{k}=rg_1;%the list of rgs created on line 87 is saved into rg_all1, indexed also by file number (cell)
    clear res dat_pos_tim xk yk rg_2;
    j=0;%returning the value of j to 0 so that the next trajectory will be analyzed from the first position
    
    average_fraction_of_outward_runs_per_cell= [average_fraction_of_outward_runs_per_cell;nanmean(fraction_plus_proc_time)]; %Emily added
    
end
   
 
     
    %below is saving the variables into particular files so that they may be
    %plotted by the next function.

    %below is saving the variables into particular files so that they may be
    %plotted by the next function.
    rg_all=rg_all1;
    save([save_dir_Rg, 'Rg_per_cell'],'rg_all');
    
    res1=res2;
    save([save_dir_MSD, 'MSD_per_cell'],'res1');
    clear rg_all res1 res2 rg_all1;
    
    outward_runs_per_cell=average_fraction_of_outward_runs_per_cell;
    save([save_dir_dirbias, 'plus_processive_runs'], 'outward_runs_per_cell');
    
end

%clear ab xk yk pos kf;

% figure('Name','Directional bias run length box and violin','NumberTitle','off'), hold on,
% violin(average_fraction_of_outward_runs_per_cell, 'xlabel',{'WT','S421D'},'facecolor',[0 0 0;0 0.4 0.4],'edgecolor','k','bw',0.1,'mc','k','medc','k--');
% h1=notBoxPlot_nodotorline(average_fraction_of_outward_runs_per_cell{1},1,'style','line');
% set(h1.data,'MarkerFaceColor','none','MarkerEdgeColor',[0 0 0], 'MarkerSize', 5);
% h2=notBoxPlot_nodotorline(average_fraction_of_outward_runs_per_cell{2},2,'style','line');
% set(h2.data,'MarkerFaceColor','none','MarkerEdgeColor',[0 0.4 0.4],'MarkerSize', 5);
% publication_fig(0,0,1)
% axisHandle = gca; 
%     set(findall(gca, 'Type', 'Line'),'LineWidth',2);
%     set(gca,'LineWidth',2);
%     set(gca,'FontSize',24);
%     set(gca, 'FontName', 'Arial');
%     set(gca,'XColor',[0 0 0],'YColor',[0 0 0]);
%     set(gca,'Box','on');
% ylabel('Fraction of outward processive runs');
% % ylim([-0.25 2]);
% hFig=findall(0,'type','figure');
% hLeg=findobj(hFig(1,1),'type','legend');
% set(hLeg,'visible','off')


toc

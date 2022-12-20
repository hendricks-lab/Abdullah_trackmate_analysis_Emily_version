%% Bootstrap MSD and posititional probability Analysis: 
%% Abdullah R. Chaudhary

%Emily updated Sept 26/2022 to add the fraction of runs moving outward,
%inward, diffusively, or stationary per cell.
%Emily updated Nov 21 2022 to ensure the variables are reinitialized each
%iteration of the loop.

close all; clear all; clc; 
set(0,'DefaultFigureWindowStyle','docked')
addpath('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/New_General_codes/');
addpath('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/AC_codes_epmodified_20200603/');
addpath('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/New_General_codes/plotSpread/');
addpath('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/violin/');
addpath('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/');
addpath('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/raacampbell-notBoxPlot-7d90c27/code/');
addpath('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/raacampbell-notBoxPlot-7d90c27/code/+NBP/');

%% Reset these parameters everytime you run the code !!!!!!
%%%%%%%%%%%%%%%%%%%%%%%%%%
alpha_MSD=[];
ab1=[];
kg=0;
variance_plus_traj=[];
variance_minus_traj=[];
proc_run_pos_2=[];
fl='*.mat';
an=0;
DT=0.12;
kf=0;
proc_runs=[];
proc_times=[];
rng=[];
tlp=-0.5;
fraction_proc_time=[]; %Emily added 20221121
fraction_plus_proc_time=[]; %Emily added
fraction_minus_proc_time=[]; %Emily added 20220927
fraction_diff_time=[]; %Emily added 20220927
fraction_stat_time=[]; %Emily added 20220927
fraction_plus_diff_time=[]; %Emily added 20221020
fraction_minus_diff_time=[]; %Emily added 20221020
fraction_proc=[]; %Emily added 20221121
fraction_proc_out=[]; %Emily added 20221121
fraction_proc_in=[]; %Emily added 20221121
fraction_diff=[]; %Emily added 20221121
fraction_diff_out=[]; %Emily added 20221121
fraction_diff_in=[]; %Emily added 20221121
fraction_stat=[]; %Emily added 20221121
fraction_proc_all=[]; %Emily added 20221121
fraction_proc_out_all=[]; %Emily added 20221121
fraction_proc_in_all=[]; %Emily added 20221121
fraction_diff_all=[]; %Emily added 20221121
fraction_diff_out_all=[]; %Emily added 20221121
fraction_diff_in_all=[]; %Emily added 20221121
fraction_stat_all=[]; %Emily added 20221121
plus_proc_times=[];
avg_frac_t_out_per_cell=[];
avg_frac_t_in_per_cell=[];
avg_frac_t_proc_per_cell=[];
avg_frac_t_diff_per_cell=[];
avg_frac_t_diff_out_per_cell=[]; %Emily added
avg_frac_t_diff_in_per_cell=[]; %Emily added
avg_frac_t_stat_per_cell=[]; %Emily added
proc_run=[]; %Emily added
diff_run=[]; %Emily added 20220927
stat_run=[]; %Emily added 20220927
proc_time=[]; %Emily added
stat_time=[]; %Emily added 20220927
diff_time=[]; % Emily added 20220927

%% Savitsky Golay smoothing
span=10;    % 10 is good, 25 is good as well
pwr=1;      % 1 is good, 2 is good
%%%%%%%%%%%%%%%%%%%%%%%%%%

tic
for k_choose = 1 %Don't forget to change min_lr based on your cargo! 0.4 for Rab5, 0.6 for Lyso
    
if k_choose == 1    % WT endogenous
    col1=[0.5 0.5 0.5];
    %Rab5
    dat_dir = '/Volumes/Emily_htt/Huntingtin_Project/EE/WT/Rab5/WT_r5_mats/';
    save_dir='/Volumes/Emily_htt/Huntingtin_Project/EE/WT/Rab5/WT_r5_pos/';
    save_dir_Rg='/Volumes/Emily_htt/Huntingtin_Project/EE/WT/Rab5/WT_r5_rg/';
    save_dir_MSD='/Volumes/Emily_htt/Huntingtin_Project/EE/WT/Rab5/WT_r5_msd/';
    save_dir_dirbias='/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/Dirbias/WT_r5_updated/';
%     %Lyso
%     dat_dir = '/Volumes/Emily_htt/Huntingtin_Project/EE/WT/Lysotracker/lyso_mats/';
%     save_dir='/Volumes/Emily_htt/Huntingtin_Project/EE/WT/Lysotracker/lyso_pos/';
%     save_dir_Rg='/Volumes/Emily_htt/Huntingtin_Project/EE/WT/Lysotracker/lyso_rg/';
%     save_dir_MSD='/Volumes/Emily_htt/Huntingtin_Project/EE/WT/Lysotracker/lyso_msd/';
%     save_dir_dirbias='/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/Dirbias/WT_lyso_updated/';
    %Qdot
%     dat_dir='/Volumes/Emily_htt/Huntingtin_Project/EE/WT/EE_Qdot/WT_all_ee_mats_tm2/';
%     save_dir='/Volumes/Emily_htt/Huntingtin_Project/EE/WT/EE_Qdot/WT_all_ee_pos_tm2/';
%     save_dir_Rg='/Volumes/Emily_htt/Huntingtin_Project/EE/WT/EE_Qdot/WT_all_ee_rg_tm2/';
%     save_dir_MSD='/Volumes/Emily_htt/Huntingtin_Project/EE/WT/EE_Qdot/WT_all_ee_msd_tm2/';
%     save_dir_dirbias='/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/Dirbias/WT_qd/';
%     % Nocodazole Rab5
%     dat_dir='/Volumes/Emily_htt_2/S421D_htt_2/Rab5/20220205_r5+NZ/20220205_r5+NZ_Tiffs/20220205_r5_control_mats/';
%     save_dir='/Volumes/Emily_htt_2/S421D_htt_2/Rab5/20220205_r5+NZ/20220205_r5+NZ_Tiffs/20220205_r5_control_pos/';
%     save_dir_Rg='/Volumes/Emily_htt_2/S421D_htt_2/Rab5/20220205_r5+NZ/20220205_r5+NZ_Tiffs/20220205_r5_control_rg/';
%     save_dir_MSD='/Volumes/Emily_htt_2/S421D_htt_2/Rab5/20220205_r5+NZ/20220205_r5+NZ_Tiffs/20220205_r5_control_msd/';
%     save_dir_dirbias='/Volumes/Emily_htt_2/S421D_htt_2/Rab5/20220205_r5+NZ/20220205_r5+NZ_Tiffs/20220205_r5_control_pos/';
    %Nocodazole Lyso
%     dat_dir='/Volumes/Emily_htt_2/S421D_htt_2/Lyso/Nocodazole/20220122_lyso_NZ/20220122_Tiffs/20220122_control_mats/';
%     save_dir='/Volumes/Emily_htt_2/S421D_htt_2/Lyso/Nocodazole/20220122_lyso_NZ/20220122_Tiffs/20220122_control_pos/';
%     save_dir_Rg='/Volumes/Emily_htt_2/S421D_htt_2/Lyso/Nocodazole/20220122_lyso_NZ/20220122_Tiffs/20220122_control_rg/';
%     save_dir_MSD='/Volumes/Emily_htt_2/S421D_htt_2/Lyso/Nocodazole/20220122_lyso_NZ/20220122_Tiffs/20220122_control_msd/';
%     save_dir_dirbias='/Volumes/Emily_htt_2/S421D_htt_2/Lyso/Nocodazole/20220122_lyso_NZ/20220122_Tiffs/20220122_control_pos/';
    
elseif k_choose == 2  % S421D endogenous
    col1='b';
    %Rab5
%     dat_dir = '/Volumes/Emily_htt/Huntingtin_Project/EE/S421D/Rab5/S421D_r5_mats/';
%     save_dir='/Volumes/Emily_htt/Huntingtin_Project/EE/S421D/Rab5/S421D_r5_pos/';
%     save_dir_Rg='/Volumes/Emily_htt/Huntingtin_Project/EE/S421D/Rab5/S421D_r5_rg/';
%     save_dir_MSD='/Volumes/Emily_htt/Huntingtin_Project/EE/S421D/Rab5/S421D_r5_msd/';
%     save_dir_dirbias='/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/Dirbias/S421D_r5_updated/';
    %Lyso
%     dat_dir = '/Volumes/Emily_htt/Huntingtin_Project/EE/S421D/Lysotracker/S421D_lyso_mats/';
%     save_dir='/Volumes/Emily_htt/Huntingtin_Project/EE/S421D/Lysotracker/S421D_lyso_pos/';
%     save_dir_Rg='/Volumes/Emily_htt/Huntingtin_Project/EE/S421D/Lysotracker/S421D_lyso_rg/';
%     save_dir_MSD='/Volumes/Emily_htt/Huntingtin_Project/EE/S421D/Lysotracker/S421D_lyso_msd/';
%     save_dir_dirbias='/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/Dirbias/S421D_lyso_updated/';
    %Qdot
%     dat_dir='/Volumes/Emily_htt/Huntingtin_Project/EE/S421D/Qdot_EE/S421D_all_ee_mats_tm2/';
%     save_dir='/Volumes/Emily_htt/Huntingtin_Project/EE/S421D/Qdot_EE/S421D_all_ee_pos_tm2/';
%     save_dir_Rg='/Volumes/Emily_htt/Huntingtin_Project/EE/S421D/Qdot_EE/S421D_all_ee_rg_tm2/';
%     save_dir_MSD='/Volumes/Emily_htt/Huntingtin_Project/EE/S421D/Qdot_EE/S421D_all_ee_msd_tm2/';
%     save_dir_dirbias='/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/Dirbias/S421D_qd/';
%     %Nocodazole Rab5
%     dat_dir='/Volumes/Emily_htt_2/S421D_htt_2/Rab5/20220205_r5+NZ/20220205_r5+NZ_Tiffs/20220205_r5_NZ_mats/';
%     save_dir='/Volumes/Emily_htt_2/S421D_htt_2/Rab5/20220205_r5+NZ/20220205_r5+NZ_Tiffs/20220205_r5_NZ_pos/';
%     save_dir_Rg='/Volumes/Emily_htt_2/S421D_htt_2/Rab5/20220205_r5+NZ/20220205_r5+NZ_Tiffs/20220205_r5_NZ_rg/';
%     save_dir_MSD='/Volumes/Emily_htt_2/S421D_htt_2/Rab5/20220205_r5+NZ/20220205_r5+NZ_Tiffs/20220205_r5_NZ_msd/';
%     save_dir_dirbias='/Volumes/Emily_htt_2/S421D_htt_2/Rab5/20220205_r5+NZ/20220205_r5+NZ_Tiffs/20220205_r5_NZ_pos/';
    %Nocodazole Lyso
% %     dat_dir='/Volumes/Emily_htt_2/S421D_htt_2/Lyso/Nocodazole/20220122_lyso_NZ/20220122_Tiffs/20220122_NZ_mats/';
% %     save_dir='/Volumes/Emily_htt_2/S421D_htt_2/Lyso/Nocodazole/20220122_lyso_NZ/20220122_Tiffs/20220122_NZ_pos/';
% %     save_dir_Rg='/Volumes/Emily_htt_2/S421D_htt_2/Lyso/Nocodazole/20220122_lyso_NZ/20220122_Tiffs/20220122_NZ_rg/';
% %     save_dir_MSD='/Volumes/Emily_htt_2/S421D_htt_2/Lyso/Nocodazole/20220122_lyso_NZ/20220122_Tiffs/20220122_NZ_msd/';
% %     save_dir_dirbias='/Volumes/Emily_htt_2/S421D_htt_2/Lyso/Nocodazole/20220122_lyso_NZ/20220122_Tiffs/20220122_NZ_pos/';

elseif k_choose == 3    % WT transfected with WT htt
    col1=[0.5 0.5 0.5];
    %Rab5
%     dat_dir = '/Volumes/Emily_htt/Huntingtin_Project/EE/WT_trans_mchWThtt/Rab5/WT_trans_WThtt_mats/';
%     save_dir='/Volumes/Emily_htt/Huntingtin_Project/EE/WT_trans_mchWThtt/Rab5/WT_trans_WThtt_pos/';
%     save_dir_Rg='/Volumes/Emily_htt/Huntingtin_Project/EE/WT_trans_mchWThtt/Rab5/WT_trans_WThtt_rg/';
%     save_dir_MSD='/Volumes/Emily_htt/Huntingtin_Project/EE/WT_trans_mchWThtt/Rab5/WT_trans_WThtt_msd/';
%     save_dir_dirbias='/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/Dirbias/WT+WThtt_r5_updated/';
    %Lyso
%     dat_dir = '/Volumes/Emily_htt/Huntingtin_Project/EE/WT_trans_mchWThtt/Lyso/WT_trans_WThtt_lyso_mats/';
%     save_dir='/Volumes/Emily_htt/Huntingtin_Project/EE/WT_trans_mchWThtt/Lyso/WT_trans_WThtt_lyso_pos/';
%     save_dir_Rg='/Volumes/Emily_htt/Huntingtin_Project/EE/WT_trans_mchWThtt/Lyso/WT_trans_WThtt_lyso_rg/';
%     save_dir_MSD='/Volumes/Emily_htt/Huntingtin_Project/EE/WT_trans_mchWThtt/Lyso/WT_trans_WThtt_lyso_msd/';
%     save_dir_dirbias='/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/Dirbias/WT+WThtt_lyso_updated/';
%     %Latrunculin A rab5
%     dat_dir='/Volumes/Emily_htt/Huntingtin_Project/EE/WT/Rab5_Latrunculin/WT_r5_latA_mats/';
%     save_dir='/Volumes/Emily_htt/Huntingtin_Project/EE/WT/Rab5_Latrunculin/WT_r5_latA_pos/';
%     save_dir_Rg='/Volumes/Emily_htt/Huntingtin_Project/EE/WT/Rab5_Latrunculin/WT_r5_latA_rg/';
%     save_dir_MSD='/Volumes/Emily_htt/Huntingtin_Project/EE/WT/Rab5_Latrunculin/WT_r5_latA_msd/';
%     save_dir_dirbias='/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/Dirbias/WT_r5_latA/';
    %Latrunculin A Qdot
%     dat_dir='/Volumes/Emily_htt/Huntingtin_Project/EE/WT/Latrunculin_Qdot/WT_latA_mats/';
%     save_dir='/Volumes/Emily_htt/Huntingtin_Project/EE/WT/Latrunculin_Qdot/WT_latA_pos/';
%     save_dir_Rg='/Volumes/Emily_htt/Huntingtin_Project/EE/WT/Latrunculin_Qdot/WT_latA_rg/';
%     save_dir_MSD='/Volumes/Emily_htt/Huntingtin_Project/EE/WT/Latrunculin_Qdot/WT_latA_msd/';
%     save_dir_dirbias='/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/Dirbias/WT_qd_latA/';
    
    elseif k_choose == 4   % WT transfected with S421D htt
    col1=[0.5 0.5 0.5];
% %     %Rab5
%     dat_dir = '/Volumes/Emily_htt/Huntingtin_Project/EE/WT_trans_S421Dhtt/Rab5/WT_trans_S421Dhtt_mats/';
%     save_dir='/Volumes/Emily_htt/Huntingtin_Project/EE/WT_trans_S421Dhtt/Rab5/WT_trans_S421Dhtt_pos/';
%     save_dir_Rg='/Volumes/Emily_htt/Huntingtin_Project/EE/WT_trans_S421Dhtt/Rab5/WT_trans_S421Dhtt_rg/';
%     save_dir_MSD='/Volumes/Emily_htt/Huntingtin_Project/EE/WT_trans_S421Dhtt/Rab5/WT_trans_S421Dhtt_msd/';
%     save_dir_dirbias='/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/Dirbias/WT+S421Dhtt_r5_updated/';
%      %Lyso
%     dat_dir = '/Volumes/Emily_htt/Huntingtin_Project/EE/WT_trans_S421Dhtt/Lyso/WT_trans_S421Dhtt_lyso_mats/';
%     save_dir='/Volumes/Emily_htt/Huntingtin_Project/EE/WT_trans_S421Dhtt/Lyso/WT_trans_S421Dhtt_lyso_pos/';
%     save_dir_Rg='/Volumes/Emily_htt/Huntingtin_Project/EE/WT_trans_S421Dhtt/Lyso/WT_trans_S421Dhtt_lyso_rg/';
%     save_dir_MSD='/Volumes/Emily_htt/Huntingtin_Project/EE/WT_trans_S421Dhtt/Lyso/WT_trans_S421Dhtt_lyso_msd/';
%     save_dir_dirbias='/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/Dirbias/WT+S421Dhtt_lyso_updated/';
    %Latrunculin A rab5
%     dat_dir='/Volumes/Emily_htt/Huntingtin_Project/EE/S421D/Rab5_Latrunculin/S421D_r5_latA_mats/';
%     save_dir='/Volumes/Emily_htt/Huntingtin_Project/EE/S421D/Rab5_Latrunculin/S421D_r5_latA_pos/';
%     save_dir_Rg='/Volumes/Emily_htt/Huntingtin_Project/EE/S421D/Rab5_Latrunculin/S421D_r5_latA_rg/';
%     save_dir_MSD='/Volumes/Emily_htt/Huntingtin_Project/EE/S421D/Rab5_Latrunculin/S421D_r5_latA_msd/';
%     save_dir_dirbias='/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/Dirbias/S421D_r5_latA/';
% %     Latrunculin A Qdot
%     dat_dir='/Volumes/Emily_htt/Huntingtin_Project/EE/S421D/Qdot_Latrunculin/S421D_latA_mats/';
%     save_dir='/Volumes/Emily_htt/Huntingtin_Project/EE/S421D/Qdot_Latrunculin/S421D_latA_pos/';
%     save_dir_Rg='/Volumes/Emily_htt/Huntingtin_Project/EE/S421D/Qdot_Latrunculin/S421D_latA_rg/';
%     save_dir_MSD='/Volumes/Emily_htt/Huntingtin_Project/EE/S421D/Qdot_Latrunculin/S421D_latA_msd/';
%     save_dir_dirbias='/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/Dirbias/S421D_qd_latA/';

end

fls=dir(fullfile(dat_dir,fl));
min_lr=0.4;

for k=1:numel(fls) %k is the current file
    pos=[];
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
        
        curr_pos=ab(j).position;
        time=[1:length(curr_pos)]'.*DT; %Emily modified 20221019 so time and position have same length (changed numel to length)

            res=analyze_run_length_reversals_v5_per_cell_bias(time,curr_pos,min_lr,0,k_choose);
            proc_run=[proc_run;res.proc_run];
            proc_time=[proc_time;res.proc_time];
        
            fraction_proc_time=[fraction_proc_time;(sum(res.proc_time))/(sum(res.proc_time)+sum(res.diff_time)+sum(res.stat_time))];
            fraction_plus_proc_time=[fraction_plus_proc_time;(sum(res.proc_time(find(res.proc_run>0)))/(sum(res.proc_time)))]; %Emily added to find plus end runs
            fraction_minus_proc_time=[fraction_minus_proc_time;(sum(res.proc_time(find(res.proc_run<0)))/(sum(res.proc_time)))]; %Emily added to find plus end runs
            fraction_diff_time=[fraction_diff_time;(sum(res.diff_time))/(sum(res.proc_time)+sum(res.diff_time)+sum(res.stat_time))]; %Emily added to find plus end runs
            fraction_stat_time=[fraction_stat_time;(sum(res.stat_time))/(sum(res.proc_time)+sum(res.diff_time)+sum(res.stat_time))]; %Emily added to find plus end runs
            %Emily added this section 20221020 for testing whether min_lr is
            %set correctly
            fraction_plus_diff_time=[fraction_plus_diff_time;(sum(res.diff_time(find(res.diff_run>0)))/sum(res.diff_time))]; %Emily added 20221020
            fraction_minus_diff_time=[fraction_minus_diff_time;(sum(res.diff_time(find(res.diff_run<0)))/sum(res.diff_time))]; %Emily added 20221020
        
        [deltat, msdpts, sem, log_deltat, log_msdpts, alpha_1, DiffCoef] = MSD_2d_emmessed(xk{kf},yk{kf}, DT, k_choose);
        if alpha_1>0 %not analyzing the tracks that are so flat they look
%         like they have negative slope I guess
        rg_0=func_Rg_Linda_v2(xk{kf},yk{kf}); 
        % (above) calculating the radius of gyration (rg) from the smoothed
        % position data of the trajectory
        rg_1=[rg_1,rg_0]; %adds the currently calculated rg to the previous one in a list. 
        dat_pos_tim{kf}=smooth(ab(j).position,10,'sgolay',1); %gives another smoothed position file for the current trajectory?
        else %if alpha is less than or equal to zero, this condition will run

        % this condition seems to only be fulfilled when there are many
        % NaNs 
        
        dat_pos_tim{kf}=[];
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
    rg_all1{k}=rg_1; %the list of rgs created on line 87 is saved into rg_all1, indexed also by file number (cell)
    clear res dat_pos_tim xk yk rg_2;
    j=0;%returning the value of j to 0 so that the next trajectory will be analyzed from the first position
    
    fraction_proc=[fraction_proc;fraction_proc_time]; %Emily added 20221121
    fraction_proc_out=[fraction_proc_out;fraction_plus_proc_time]; %Emily added 20221121
    fraction_proc_in=[fraction_proc_in;fraction_minus_proc_time]; %Emily added 20221121
    fraction_diff=[fraction_diff;fraction_diff_time]; %Emily added 20221121
    fraction_diff_out=[fraction_diff_out;fraction_plus_diff_time]; %Emily added 20221121
    fraction_diff_in=[fraction_diff_in;fraction_minus_diff_time]; %Emily added 20221121
    fraction_stat=[fraction_stat;fraction_stat_time]; %Emily added 20221121
    
    avg_frac_t_proc_per_cell= [avg_frac_t_proc_per_cell;nanmean(fraction_proc_time)]; %Emily added
    avg_frac_t_out_per_cell= [avg_frac_t_out_per_cell;nanmean(fraction_plus_proc_time)]; %Emily added
    avg_frac_t_in_per_cell= [avg_frac_t_in_per_cell;nanmean(fraction_minus_proc_time)]; %Emily added
    avg_frac_t_diff_per_cell= [avg_frac_t_diff_per_cell;nanmean(fraction_diff_time)]; %Emily added
    avg_frac_t_diff_out_per_cell= [avg_frac_t_diff_out_per_cell;nanmean(fraction_plus_diff_time)]; %Emily added
    avg_frac_t_diff_in_per_cell= [avg_frac_t_diff_in_per_cell;nanmean(fraction_minus_diff_time)]; %Emily added
    avg_frac_t_stat_per_cell= [avg_frac_t_stat_per_cell;nanmean(fraction_stat_time)]; %Emily added
    
    
%     fraction_proc_per_cell{k}=fraction_proc; %Emily added 20221121
%     fraction_proc_out_per_cell{k}=fraction_proc_out; %Emily added 20221121
%     fraction_proc_in_per_cell{k}=fraction_proc_in; %Emily added 20221121
%     fraction_diff_per_cell{k}=fraction_diff; %Emily added 20221121
%     fraction_diff_out_per_cell{k}=fraction_diff_out; %Emily added 20221121
%     fraction_diff_in_per_cell{k}=fraction_diff_in; %Emily added 20221121
%     fraction_stat_per_cell{k}=fraction_stat; %Emily added 20221121
    
    fraction_proc_time=[]; %Emily added 20221121
    fraction_plus_proc_time=[]; %Emily added 20221121
    fraction_minus_proc_time=[]; %Emily added 20221121
    fraction_diff_time=[]; %Emily added 20221121
    fraction_plus_diff_time=[]; %Emily added 20221121
    fraction_minus_diff_time=[]; %Emily added 20221121
    fraction_stat_time=[]; %Emily added 20221121
    
end

rg_1=[];

%Getting a list of the fractions processive, diffusive, and stationary for
%all trajectories
fraction_proc_all=[fraction_proc_all;fraction_proc]; %Emily added 20221121
fraction_proc_out_all=[fraction_proc_out_all;fraction_proc_out]; %Emily added 20221121
fraction_proc_in_all=[fraction_proc_in_all;fraction_proc_in]; %Emily added 20221121
fraction_diff_all=[fraction_diff_all;fraction_diff]; %Emily added 20221121
fraction_diff_out_all=[fraction_diff_out_all;fraction_diff_out]; %Emily added 20221121
fraction_diff_in_all=[fraction_diff_in_all;fraction_diff_in]; %Emily added 20221121
fraction_stat_all=[fraction_stat_all;fraction_stat]; %Emily added 20221121

fraction_proc=[]; %Emily added 20221121
fraction_proc_out=[]; %Emily added 20221121
fraction_proc_in=[]; %Emily added 20221121
fraction_diff=[]; %Emily added 20221121
fraction_diff_out=[]; %Emily added 20221121
fraction_diff_in=[]; %Emily added 20221121
fraction_stat=[]; %Emily added 20221121

% Defining variables to save
    rg_all=rg_all1;
    res1=res2;
    t_proc_per_cell=avg_frac_t_proc_per_cell;
    t_proc_out_per_cell=avg_frac_t_out_per_cell;
    t_proc_in_per_cell=avg_frac_t_in_per_cell;
    t_diff_per_cell=avg_frac_t_diff_per_cell;
    t_diff_out_per_cell=avg_frac_t_diff_out_per_cell;
    t_diff_in_per_cell= avg_frac_t_diff_in_per_cell;
    t_stat_per_cell=avg_frac_t_stat_per_cell;
    frac_proc_all=fraction_proc_all;
    frac_proc_out_all=fraction_proc_out_all;
    frac_proc_in_all=fraction_proc_in_all;
    frac_diff_all=fraction_diff_all;
    frac_diff_out_all=fraction_diff_out_all;
    frac_diff_in_all= fraction_diff_in_all;
    frac_stat_all=fraction_stat_all;
% Saving variables
    save([save_dir_Rg, 'Rg_per_cell'],'rg_all');
    save([save_dir_MSD, 'MSD_per_cell'],'res1');
    save([save_dir_dirbias, 'frac_proc_per_cell'], 't_proc_per_cell');
    save([save_dir_dirbias, 'frac_proc_out_per_cell'], 't_proc_out_per_cell');
    save([save_dir_dirbias, 'frac_proc_in_per_cell'], 't_proc_in_per_cell');
    save([save_dir_dirbias, 'frac_diff_per_cell'], 't_diff_per_cell');
    save([save_dir_dirbias, 'frac_diff_out_per_cell'], 't_diff_out_per_cell');
    save([save_dir_dirbias, 'frac_diff_in_per_cell'], 't_diff_in_per_cell');
    save([save_dir_dirbias, 'frac_stat_per_cell'], 't_stat_per_cell');
    save([save_dir_dirbias, 'proc_all'], 'frac_proc_all');
    save([save_dir_dirbias, 'proc_out_all'], 'frac_proc_out_all');
    save([save_dir_dirbias, 'proc_in_all'], 'frac_proc_in_all');
    save([save_dir_dirbias, 'diff_all'], 'frac_diff_all');
    save([save_dir_dirbias, 'diff_out_all'], 'frac_diff_out_all');
    save([save_dir_dirbias, 'diff_in_all'], 'frac_diff_in_all');
    save([save_dir_dirbias, 'stat_all'], 'frac_stat_all');
% Clearing all variables after running the code
    clear rg_all res1 res2 rg_all1;
    clear fraction_proc_per_cell;
    clear fraction_proc_out_per_cell;
    clear fraction_proc_in_per_cell;
    clear fraction_diff_per_cell;
    clear fraction_diff_out_per_cell;
    clear fraction_diff_in_per_cell;
    clear fraction_stat_per_cell;
    
end

clear ab xk yk pos kf;

toc

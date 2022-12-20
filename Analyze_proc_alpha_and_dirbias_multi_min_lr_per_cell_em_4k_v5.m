%% Bootstrap MSD and posititional probability Analysis: 
%% Abdullah R. Chaudhary

%Emily updated Sept 26/2022 to add the fraction of runs moving outward,
%inward, diffusively, or stationary per cell.
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

%% Savitsky Golay smoothing
span=10;    % 10 is good, 25 is good as well
pwr=1;      % 1 is good, 2 is good
%%%%%%%%%%%%%%%%%%%%%%%%%%

tic
k_choose = 1
    
if k_choose == 1    % WT endogenous
    col1=[0.5 0.5 0.5];
    %Rab5
    dat_dir = '/Volumes/Emily_htt/Huntingtin_Project/EE/WT/Rab5/WT_r5_mats/';
    save_dir_dirbias='/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/Dirbias/WT_r5_updated/';
     %Lyso
%     dat_dir = '/Volumes/Emily_htt/Huntingtin_Project/EE/WT/Lysotracker/lyso_mats/';
%     save_dir_dirbias='/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/Dirbias/WT_lyso_updated/';
    %Qdot
%     dat_dir='/Volumes/Emily_htt/Huntingtin_Project/EE/WT/EE_Qdot/WT_all_ee_mats_tm2/';
%     save_dir_dirbias='/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/Dirbias/WT_qd/';
    % Nocodazole Rab5
%     dat_dir='/Volumes/Emily_htt_2/S421D_htt_2/Rab5/20220205_r5+NZ/20220205_r5+NZ_Tiffs/20220205_r5_control_mats/';
%     save_dir_dirbias='/Volumes/Emily_htt_2/S421D_htt_2/Rab5/20220205_r5+NZ/20220205_r5+NZ_Tiffs/20220205_r5_control_pos/';
    %Nocodazole Lyso
%     dat_dir='/Volumes/Emily_htt_2/S421D_htt_2/Lyso/Nocodazole/20220122_lyso_NZ/20220122_Tiffs/20220122_control_mats/';
%     save_dir_dirbias='/Volumes/Emily_htt_2/S421D_htt_2/Lyso/Nocodazole/20220122_lyso_NZ/20220122_Tiffs/20220122_control_pos/';
   
elseif k_choose == 2  % S421D endogenous
    col1='b';
    %Rab5
    dat_dir = '/Volumes/Emily_htt/Huntingtin_Project/EE/S421D/Rab5/S421D_r5_mats/';
    save_dir_dirbias='/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/Dirbias/S421D_r5_updated/';
    %Lyso
%     dat_dir = '/Volumes/Emily_htt/Huntingtin_Project/EE/S421D/Lysotracker/S421D_lyso_mats/';
%     save_dir_dirbias='/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/Dirbias/S421D_lyso_updated/';
    %Qdot
%     dat_dir='/Volumes/Emily_htt/Huntingtin_Project/EE/S421D/Qdot_EE/S421D_all_ee_mats_tm2/';
%     save_dir_dirbias='/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/Dirbias/S421D_qd/';
    %Nocodazole Rab5
%     dat_dir='/Volumes/Emily_htt_2/S421D_htt_2/Rab5/20220205_r5+NZ/20220205_r5+NZ_Tiffs/20220205_r5_NZ_mats/';
%     save_dir_dirbias='/Volumes/Emily_htt_2/S421D_htt_2/Rab5/20220205_r5+NZ/20220205_r5+NZ_Tiffs/20220205_r5_NZ_pos/';
    %Nocodazole Lyso
%     dat_dir='/Volumes/Emily_htt_2/S421D_htt_2/Lyso/Nocodazole/20220122_lyso_NZ/20220122_Tiffs/20220122_NZ_mats/';
%     save_dir_dirbias='/Volumes/Emily_htt_2/S421D_htt_2/Lyso/Nocodazole/20220122_lyso_NZ/20220122_Tiffs/20220122_NZ_pos/';
    
elseif k_choose == 3    % WT transfected with WT htt
    col1=[0.5 0.5 0.5];
    %Rab5
    dat_dir = '/Volumes/Emily_htt/Huntingtin_Project/EE/WT_trans_mchWThtt/Rab5/WT_trans_WThtt_mats/';
    save_dir_dirbias='/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/Dirbias/WT+WThtt_r5_updated/';
    %Lyso
%     dat_dir = '/Volumes/Emily_htt/Huntingtin_Project/EE/WT_trans_mchWThtt/Lyso/WT_trans_WThtt_lyso_mats/';
%     save_dir_dirbias='/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/Dirbias/WT+WThtt_lyso_updated/';
    %Latrunculin A rab5
%     dat_dir='/Volumes/Emily_htt/Huntingtin_Project/EE/WT/Rab5_Latrunculin/WT_r5_latA_mats/';
%     save_dir_dirbias='/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/Dirbias/WT_r5_latA/';
    %Latrunculin A Qdot
%     dat_dir='/Volumes/Emily_htt/Huntingtin_Project/EE/WT/Latrunculin_Qdot/WT_latA_mats/';
%     save_dir_dirbias='/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/Dirbias/WT_qd_latA/';
    
    elseif k_choose == 4   % WT transfected with S421D htt
    col1=[0.5 0.5 0.5];
    %Rab5
    dat_dir = '/Volumes/Emily_htt/Huntingtin_Project/EE/WT_trans_S421Dhtt/Rab5/WT_trans_S421Dhtt_mats/';
    save_dir_dirbias='/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/Dirbias/WT+S421Dhtt_r5_updated/';
%      Lyso
%     dat_dir = '/Volumes/Emily_htt/Huntingtin_Project/EE/WT_trans_S421Dhtt/Lyso/WT_trans_S421Dhtt_lyso_mats/';
%     save_dir_dirbias='/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/Dirbias/WT+S421Dhtt_lyso_updated/';
    %Latrunculin A rab5
%     dat_dir='/Volumes/Emily_htt/Huntingtin_Project/EE/S421D/Rab5_Latrunculin/S421D_r5_latA_mats/';
%     save_dir_dirbias='/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/Dirbias/S421D_r5_latA/';
    %Latrunculin A Qdot
%     dat_dir='/Volumes/Emily_htt/Huntingtin_Project/EE/S421D/Qdot_Latrunculin/S421D_latA_mats/';
%     save_dir_dirbias='/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/Dirbias/S421D_qd_latA/';

end

ab1=[];
fl='*.mat';
DT=0.12;
all_proc_trajs=[];
all_diff_trajs=[];
proc_runs=[];
diff_runs=[];
stat_runs=[];
proc_alph1=[];
proc_alpha=[];
proc_alpha_mthd_2=[];  %Emily added 20221125
diff_alpha_mthd_2=[];  %Emily added 20221125
diff_alph1=[];
diff_alpha=[];
fraction_proc_time=[]; %Emily 20221102
fraction_plus_proc_time=[]; %Emily added
fraction_minus_proc_time=[]; %Emily added 20220927
fraction_diff_time=[]; %Emily added 20220927
fraction_stat_time=[]; %Emily added 20220927
fraction_plus_diff_time=[]; %Emily added 20221020
fraction_minus_diff_time=[]; %Emily added 20221020
proc_run=[]; %Emily added

fls=dir(fullfile(dat_dir,fl));

for iteration=1:25 %Emily turned into a loop 20221031
min_lr=0.1*iteration;
kf=0; %Emily moved to here so that these reset for each iteration 20221116
kg=0; %Emily moved to here so that these reset for each iteration 20221116

    for k=1:numel(fls) %k is the current file aka cell
        pos=[];
        kg=kg+1;%keeping track of the number of cells
        load(fullfile(dat_dir,fls(k).name), 'ab'); %load filename
        ab = ab(~cellfun(@isempty,{ab.position})); % outputs a logical array of if the position is empty (0 for empty) or not.
        ab1{k}=ab; %Emily replaced the indexing by kg with indexing by k instead 20221116
        rg_1=[];
    
    
        for j=1:numel(ab) %from 1 to the number of elements in the position file, which would be the number of trajectories for a file
            if numel(ab(j).position)>15 %only analyze trajectories that have more than 15 timepoints of positional data
            kf=kf+1; %keeping track of the number of trajectories for all cells 
            
            curr_pos=ab(j).position;
            time=(1:length(curr_pos))'.*DT; %Emily modified 20221019 so time and position have same length (changed numel to length)
                res=analyze_run_length_reversals_v5_per_cell_bias(time,curr_pos,min_lr,0,k_choose);
            
                fraction_proc_time=[fraction_proc_time;sum(res.proc_time)/(sum(res.proc_time)+sum(res.diff_time)+sum(res.stat_time))];
                fraction_diff_time=[fraction_diff_time;(sum(res.diff_time))/(sum(res.proc_time)+sum(res.diff_time)+sum(res.stat_time))]; %Emily added to find fraction of diffusive runs
                fraction_stat_time=[fraction_stat_time;(sum(res.stat_time))/(sum(res.proc_time)+sum(res.diff_time)+sum(res.stat_time))]; %Emily added to find fraction of stationary runs
                
                proc_runs=[proc_runs;res.proc_run]; %Emily added 20221116 to calculate number of runs
                diff_runs=[diff_runs;res.diff_run]; %Emily added 20221116 to calculate number of runs
                stat_runs=[stat_runs;res.stat_run]; %Emily added 20221116 to calculate number of runs
                
                fraction_plus_proc_time=[fraction_plus_proc_time;(sum(res.proc_time(find(res.proc_run>0)))/(sum(res.proc_time)))]; %Emily added to find plus end runs
                fraction_minus_proc_time=[fraction_minus_proc_time;(sum(res.proc_time(find(res.proc_run<0)))/(sum(res.proc_time)))]; %Emily added to find plus end runs

                fraction_plus_diff_time=[fraction_plus_diff_time;(sum(res.diff_time(find(res.diff_run>0)))/sum(res.diff_time))]; %Emily added 20221020
                fraction_minus_diff_time=[fraction_minus_diff_time;(sum(res.diff_time(find(res.diff_run<0)))/sum(res.diff_time))]; %Emily added 20221020
            
                %Emily new section 20221026 calculating alpha for proc runs
                if numel(res.proc_traj)>=1
                    for num_proc_runs=1:numel(res.proc_traj)
                        
                        proc_traj=cell2mat(res.proc_traj(num_proc_runs));
                        
                        %Getting the position values from each run into a
                        %cell array (20221124)
                        if length(proc_traj)>1
                            for proc_run_length=1:length(proc_traj)
                                proc_pos_per_run{num_proc_runs,proc_run_length}=proc_traj(proc_run_length);
                            end
                        else
                            proc_pos_per_run{num_diff_runs,1}=NaN;
                        end
                        %Emily original calculation of msd
                        if numel(proc_traj)>=4
                            proc_alph=msd_fun_v1_emmessedwithit_2(1:(numel(proc_traj)-2),proc_traj,DT,numel(proc_traj),0);
                            proc_alph1=[proc_alph1;proc_alph.slp]; %alphas for each cell
                        end
                    end
                    proc_alpha = [proc_alpha; proc_alph1]; %Adam changed 20221115, for all the trajectories
                    proc_alph1=[];
                    proc_traj=[]; %Emily added 20221125
                    proc_pos_per_traj.traj{j}=proc_pos_per_run;
                    proc_pos_per_run=[];
                end
                
                %Getting the average position value for each run at each time point for processive trajectories 20221124
                size_proc_pos_array=size(proc_pos_per_traj.traj{j});
                for proc_run_length=1:size_proc_pos_array(2)
                    avg_proc_pos_each_t{j}(proc_run_length)=mean(cell2mat(proc_pos_per_traj.traj{j}(:, proc_run_length)));
                end
                
                %tood is the best
                
                %Getting the alpha values after averaging run positions first (Method 2)
                if size_proc_pos_array(2)>=4
                    proc_alpha_mthd_2_calc=msd_fun_v1_emmessedwithit_2(1:(numel(avg_proc_pos_each_t{j})-2),avg_proc_pos_each_t{j},DT,numel(avg_proc_pos_each_t{j}),0);
                    proc_alpha_mthd_2=[proc_alpha_mthd_2;proc_alpha_mthd_2_calc.slp];
                end
                
                %Emily new section 20221124 calculating alpha for diff runs
                if numel(res.diff_traj)>=1
                    for num_diff_runs=1:numel(res.diff_traj)
                        diff_traj=cell2mat(res.diff_traj(num_diff_runs));
                        %Emily trying a different way of calculating alpha
                        %20221124
                        if length(diff_traj)>1
                            for diff_run_length=1:length(diff_traj)
                                diff_pos_per_run{num_diff_runs,diff_run_length}=diff_traj(diff_run_length);
                            end
                        else
                            diff_pos_per_run{num_diff_runs,1}=NaN;
                        end
                    %Emily's original method for calculating alpha
                        if numel(diff_traj)>=4
                            diff_alph=msd_fun_v1_emmessedwithit_2(1:(numel(diff_traj)-2),diff_traj,DT,numel(diff_traj),0);
                            diff_alph1=[diff_alph1;diff_alph.slp]; %alphas for each cell
                        end
                    end
                    diff_alpha = [diff_alpha; diff_alph1]; %Adam changed 20221115, for all the trajectories
                    diff_alph1=[];
                    diff_traj=[];
                    diff_pos_per_traj.traj{j}=diff_pos_per_run; %Emily added 20221125
                    diff_pos_per_run=[]; %Emily added 20221125
                end
                
                %Getting the average position value for each run at each time point for diffusive trajectories 20221124
                size_diff_pos_array=size(diff_pos_per_traj.traj{j});
                for diff_run_length=1:size_diff_pos_array(2)
                    avg_diff_pos_each_t{j}(diff_run_length)=mean(cell2mat(diff_pos_per_traj.traj{j}(:, diff_run_length)));
                end
                
                %Calculating the alpha for diffusive runs for all trajectories
                if size_diff_pos_array(2)>=4
                    diff_alpha_mthd_2_calc=msd_fun_v1_emmessedwithit_2(1:(numel(avg_diff_pos_each_t{j})-2),avg_diff_pos_each_t{j},DT,numel(avg_diff_pos_each_t{j}),0);
                    diff_alpha_mthd_2=[proc_alpha_mthd_2;proc_alpha_mthd_2_calc.slp];
                end
                
%                 all_proc_trajs{j}=res.proc_traj; %Emily removed 20221125
%                 all_diff_trajs{j}=res.diff_traj; %Emily removed 20221125
                
            else %do not analyze trajectories with <15 timepoints of positions.
            end    
            j=0;%returning the value of j to 0 so that the next trajectory will be analyzed from the first position
        
        end %end of trajectory analysis
        proc_alpha_mthd_2_per_cell{k}=proc_alpha_mthd_2;  %Emily added 20221125
        diff_alpha_mthd_2_per_cell{k}=diff_alpha_mthd_2;  %Emily added 20221125
        proc_alpha_per_cell{k}=proc_alpha;
        diff_alpha_per_cell{k}=diff_alpha; %Emily added 20221124
%         all_proc_traj_per_cell{k}=all_proc_trajs; %Emily added 20221124,removed
%         all_diff_traj_per_cell{k}=all_diff_trajs; %Emily added 20221124,removed
        fraction_proc_per_cell{k}=fraction_proc_time;
        fraction_proc_out_per_cell{k}=fraction_plus_proc_time;
        fraction_proc_in_per_cell{k}=fraction_minus_proc_time;
        fraction_diff_per_cell{k}=fraction_diff_time;
        fraction_diff_out_per_cell{k}=fraction_plus_diff_time;
        fraction_diff_in_per_cell{k}=fraction_minus_diff_time;
        fraction_stat_per_cell{k}=fraction_stat_time;
        num_proc_per_cell{k}=numel(proc_runs); %Emily added 20221116
        num_diff_per_cell{k}=numel(diff_runs); %Emily added 20221116
        num_stat_per_cell{k}=numel(stat_runs); %Emily added 20221116
        
        proc_alpha_mthd_2=[]; %Emily added 20221125
        diff_alpha_mthd_2=[]; %Emily added 20221125
        proc_alpha=[];
        diff_alpha=[]; %Emily added 20221124
%         all_proc_trajs=[]; %Emily added 20221124
%         all_diff_trajs=[]; %Emily added 20221124
        fraction_proc_time=[];
        fraction_plus_proc_time=[];
        fraction_minus_proc_time=[];
        fraction_diff_time=[];
        fraction_plus_diff_time=[];
        fraction_minus_diff_time=[];
        fraction_stat_time=[];
        proc_run=[]; %Emily added 20221114
        proc_time=[]; %Emily added 20221114
        res=[]; %Emily added 20221114
        proc_runs=[]; %Emily added 20221116
        diff_runs=[]; %Emily added 20221116
        stat_runs=[]; %Emily added 20221116
        
    end %End of cell analysis

     results{iteration}.proc_out=fraction_proc_out_per_cell;
     results{iteration}.diff_out=fraction_diff_out_per_cell;    
     results{iteration}.proc_alpha=proc_alpha_per_cell;
     results{iteration}.diff_alpha=diff_alpha_per_cell; %Emily added 20221124
     results{iteration}.num_proc=num_proc_per_cell;
     results{iteration}.num_diff=num_diff_per_cell;
     results{iteration}.num_stat=num_stat_per_cell;
     results{iteration}.proc_alpha_mthd_2=proc_alpha_mthd_2_per_cell; %Emily added 20221125
     results{iteration}.diff_alpha_mthd_2=diff_alpha_mthd_2_per_cell; %Emily added 20221125
%      results{iteration}.proc_trajs_all=all_proc_traj_per_cell; %Emily added 20221124
%      results{iteration}.diff_trajs_all=all_diff_traj_per_cell; %Emily added 20221124
     
    proc_alph_filtered=[];
    for cell=1:numel(proc_alpha_per_cell)
        if ~isempty(proc_alpha_per_cell(cell))
        proc_alph_filtered=[proc_alph_filtered;proc_alpha_per_cell(cell)];   
        else
        end
    end
    
    diff_alph_filtered=[];
    for cell=1:numel(diff_alpha_per_cell)
        if ~isempty(diff_alpha_per_cell(cell))
        diff_alph_filtered=[diff_alph_filtered;diff_alpha_per_cell(cell)];   
        else
        end
    end
    
    proc_alph_mthd_2_filtered=[];
    for cell=1:numel(proc_alpha_mthd_2_per_cell)
        if ~isempty(proc_alpha_mthd_2_per_cell(cell))
        proc_alph_mthd_2_filtered=[proc_alph_mthd_2_filtered;proc_alpha_mthd_2_per_cell(cell)];   
        else
        end
    end
    
    diff_alph_mthd_2_filtered=[];
    for cell=1:numel(diff_alpha_mthd_2_per_cell)
        if ~isempty(diff_alpha_mthd_2_per_cell(cell))
        diff_alph_mthd_2_filtered=[diff_alph_mthd_2_filtered;diff_alpha_mthd_2_per_cell(cell)];   
        else
        end
    end
    
    frac_diff_out_filtered=[];
    for cell=1:numel(fraction_diff_out_per_cell)
        if ~isempty(fraction_diff_out_per_cell)
            frac_diff_out_filtered=[frac_diff_out_filtered;fraction_diff_out_per_cell(cell)];
        else
        end
    end
    
    frac_proc_out_filtered=[];
    for cell=1:numel(fraction_proc_out_per_cell)
        if ~isempty(fraction_diff_out_per_cell)
            frac_proc_out_filtered=[frac_proc_out_filtered;fraction_proc_out_per_cell(cell)]; 
            else
        end
    end
    
    proc_alpha_w_nans=cell2mat(proc_alph_filtered);
    diff_alpha_w_nans=cell2mat(diff_alph_filtered);
    proc_alpha_mthd_2_w_nans=cell2mat(proc_alph_filtered);
    diff_alpha_mthd_2_w_nans=cell2mat(diff_alph_filtered);
    frac_proc_out_mat_w_nans=cell2mat(frac_proc_out_filtered);
    frac_diff_out_mat_w_nans=cell2mat(frac_diff_out_filtered);
    
    proc_alpha_mat=proc_alpha_w_nans(~isnan(proc_alpha_w_nans));
    diff_alpha_mat=diff_alpha_w_nans(~isnan(diff_alpha_w_nans));
    proc_alpha_mthd_2_mat=proc_alpha_w_nans(~isnan(proc_alpha_mthd_2_w_nans));
    diff_alpha_mthd_2_mat=diff_alpha_w_nans(~isnan(diff_alpha_mthd_2_w_nans));
    frac_proc_out_mat=frac_proc_out_mat_w_nans(~isnan(frac_proc_out_mat_w_nans));
    frac_diff_out_mat= frac_diff_out_mat_w_nans(~isnan(frac_diff_out_mat_w_nans));

    avg_proc_alpha=mean(proc_alpha_mat);
    avg_diff_alpha=mean(diff_alpha_mat);
    avg_proc_alpha_mthd_2=mean(proc_alpha_mthd_2_mat);
    avg_diff_alpha_mthd_2=mean(diff_alpha_mthd_2_mat);
    avg_proc_out=mean(frac_proc_out_mat);
    avg_diff_out=mean(frac_diff_out_mat);
    
    sem_proc_alpha=std(proc_alpha_mat)/sqrt(length(proc_alpha_mat));
    sem_diff_alpha=std(diff_alpha_mat)/sqrt(length(diff_alpha_mat));
    sem_proc_alpha_mthd_2=std(proc_alpha_mthd_2_mat)/sqrt(length(proc_alpha_mthd_2_mat));
    sem_diff_alpha_mthd_2=std(diff_alpha_mthd_2_mat)/sqrt(length(diff_alpha_mthd_2_mat));
    sem_frac_proc_out=std(frac_proc_out_mat)/sqrt(length(frac_proc_out_mat));
    sem_frac_diff_out=std(frac_diff_out_mat)/sqrt(length(frac_diff_out_mat));
    
    %     nbins=50;
%     
%     figure(200), hold on,
%     histogram(proc_alpha1,nbins,'Normalization','probability')
%     xlabel('Alpha for Processive Runs');
%     ylabel('Frequency');
%     publication_fig(0,0,1);
%     nbins=50;
%     
%     figure(200), hold on,
%     histogram(diff_alpha1,nbins,'Normalization','probability')
%     xlabel('Alpha for Processive Runs');
%     ylabel('Frequency');
%     publication_fig(0,0,1);
    
    figure(50), hold on,
    plot(min_lr, avg_proc_alpha, 'ro');
    errorbar(min_lr,avg_proc_alpha,sem_proc_alpha,'k-');
    plot(min_lr,avg_diff_alpha,'bo');
    errorbar(min_lr,avg_diff_alpha,sem_diff_alpha,'k-');
    xlabel('Minimum Run Length (\mum)');
    ylabel('Average Run Alpha');
    publication_fig(0,0,1);
    
    figure(55), hold on,
    plot(min_lr, avg_proc_alpha_mthd_2, 'ro');
    errorbar(min_lr,avg_proc_alpha_mthd_2,sem_proc_alpha_mthd_2,'k-');
    plot(min_lr,avg_diff_alpha_mthd_2,'bo');
    errorbar(min_lr,avg_diff_alpha_mthd_2,sem_diff_alpha_mthd_2,'k-');
    xlabel('Minimum Run Length (\mum)');
    ylabel('Average Run Alpha');
    publication_fig(0,0,1);
    
    figure(100), hold on,
    plot(min_lr, avg_proc_out, 'ro');
    plot(min_lr,avg_diff_out,'b*');
    errorbar(min_lr,avg_proc_out,sem_frac_proc_out,'k-');
    errorbar(min_lr,avg_diff_out,sem_frac_diff_out,'k-');
    xlabel('Minimum Run Length (\mum)');
    ylabel('Fraction of Time of Runs Moving Outward');
    publication_fig(0,0,1);
    
    clear proc_alpha_per_cell proc_alph_filtered;
    clear diff_alpha_per_cell diff_alph_filtered;
    clear proc_alpha_mthd_2_per_cell proc_alph_mthd_2_filtered;
    clear diff_alpha_mthd_2_per_cell diff_alph_mthd_2_filtered;
    clear fraction_proc_per_cell fraction_proc_out_per_cell fraction_proc_in_per_cell_min_lrs fraction_diff_per_cell fraction_diff_out_per_cell fraction_diff_in_per_cell fraction_stat_per_cell;
    clear num_proc_per_cell num_diff_per_cell num_stat_per_cell;
end %end of min lr analysis
% save([save_dir_dirbias,'frac_proc_alph_diff_per_min_lr_v2'], 'results');

toc

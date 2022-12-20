function [res] = analyze_run_length_reversals_v4(time_in,position_in,min_lr,plt,k_choose) %Comment out function for testing


% Its 0.03 for Early Endosomes and lysosomes
%% Savitsky Golay smoothing
span=9;    % 10 is good, 25 is good as well
pwr=2;      % 1 is good, 2 is good

% Emily adding for testing
% addpath('/Volumes/Emily_htt/Huntingtin_Project/EE/WT/Rab5/WT_r5_mats/');
% load('1Spots in tracks statistics_20200926_WT_r5_s1_1004-1.mat', 'ab');
% DT=0.12;
% position_in=ab(1).position;
% time_in=(1:length(position_in)).*DT;
% plt=1;
% min_lr=0.4;
% testing section over


position_smooth=smooth(position_in,span,'sgolay',pwr);
timek=time_in;
% if plt==1
% figure(989.*k_choose),hold on,
% plot(position_in)
% hold on, 
% plot(position_smooth)
% ylim([-6 6])
% hold off
% else
% end

[lr_ind_p lr_ind_v]=findinflections(0,position_smooth);
lr_ind = sort([lr_ind_p; lr_ind_v]);
lr_1=position_smooth(lr_ind);
ttim=time_in(lr_ind);
run_bw_reversal=diff(lr_1);
timek_bw_reversal = diff(ttim);
rev_rate=numel(lr_ind)./(timek(end)-timek(1));

if plt==1
figure()
plot(time_in,position_in, 'g-.'), hold on
plot(time_in,position_smooth,'r-'), hold on
plot(ttim,lr_1, 'bo')
hold off
else 
end

proc_run=run_bw_reversal(find(abs(run_bw_reversal)>=min_lr));
diff_run=run_bw_reversal(find(0.03<=abs(run_bw_reversal)& abs(run_bw_reversal)<min_lr)); %Emily modified 20221116 based on matlab warning
stat_run=run_bw_reversal(find(abs(run_bw_reversal)<0.03));

proc_time=timek_bw_reversal(find(abs(run_bw_reversal)>=min_lr));
diff_time=timek_bw_reversal(find(0.03<=abs(run_bw_reversal)& abs(run_bw_reversal)<min_lr)); %Emily modified 20221116 based on matlab warning
stat_time=timek_bw_reversal(find(abs(run_bw_reversal)<0.03));

proc_plus_runs=proc_run(find(proc_run>0));
proc_minus_runs=proc_run(find(proc_run<0));
proc_plus_times=proc_time(find(proc_run>0));
proc_minus_times=proc_time(find(proc_run<0));

%Emily modified function 20221026 to determine alpha of proc runs (getting the
%position values in this function)

proc_run_ind=find(abs(run_bw_reversal)>=min_lr);
run_starts=proc_run_ind(1:(end-1));
run_ends=proc_run_ind(2:end);

count=1;
if numel(proc_run_ind)>1 %&& proc_run_ind(1)+1 ~=proc_run_ind(2)  %Emily removed && condition 20221124
    for p=1:(numel(proc_run_ind)-1)
        res.proc_trajs{p}=position_smooth(run_starts(p):run_ends(p));
        count=p+1;
    end
else
    res.proc_trajs{count}=[];
end

%Emily adding new section to determine alpha of diffusive runs (getting the
%position values in this function)
diff_run_ind=find(0.03<=abs(run_bw_reversal)& abs(run_bw_reversal)<min_lr);
diff_run_starts=diff_run_ind(1:(end-1));
diff_run_ends=diff_run_ind(2:end);
count_2=1;
if numel(diff_run_ind)>1 %&& diff_run_ind(1)+1 ~=diff_run_ind(2) %Emily removed && condition 20221124
    for q=1:(numel(diff_run_ind)-1)
        res.diff_trajs{q}=position_smooth(diff_run_starts(q):diff_run_ends(q));
        count_2=q+1;
    end
else
    res.diff_trajs{count_2}=[];
end


res.run_bw_rev=run_bw_reversal;
res.time_bw_rev=timek_bw_reversal;
res.reversal_rate=rev_rate;
res.proc_run=proc_run;
res.diff_run=diff_run;
res.stat_run=stat_run;
res.proc_time=proc_time;
res.diff_time=diff_time;
res.stat_time=stat_time;
res.avg_velplus=mean(proc_plus_runs./proc_plus_times);
res.avg_velminus=mean(proc_minus_runs./proc_minus_times);
res.avg_lrplus=mean(proc_plus_runs);
res.avg_lrminus=mean(proc_minus_runs);
res.proc_traj=res.proc_trajs;
res.diff_traj=res.diff_trajs;



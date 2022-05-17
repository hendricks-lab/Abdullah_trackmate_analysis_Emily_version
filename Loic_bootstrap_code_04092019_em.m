% This code is written by Loïc Chaubet for Linda, drawing on the
% bstrp_func_loic_improved function that is found in another file of the
% same name as the function. 04 - sept - 2019
% This code does the bootstrapping analysis for comparing the distribution
% and statisical significance of Rg values between two conditions.
%
% all_Rg = [Rg, Itau]
% iterations = number of bootstrap iterations [1000]
% Pbin = the fraction of the data points in each bin [0.1]
% Nwindow = the number of bins included in the sliding window [3]
% alpha = (100-alpha)% is the confidence interval [5]
%nsamp= the number of times to take the sample, based on the minimum number
%of rg trajectories available (the smallest dataset).


function Rg_bstrap = Loic_bootstrap_code_04092019(all_Rg,nsamp,iterations, alpha);

all_Rg=all_Rg.';
zp=0; %make plots?
  %The original binning
  dbin=0.05;
%   jkeep = find(all_Rg(:,2)>0); 
%   bin_quartile = quantile(all_Rg(jkeep,2),[Pbin:Pbin:1]);
%   jb=1;
%   jb_max=Nwindow;
%   bin_index{jb}=find(all_Rg(:,2)==0); %make first bin for I = 0
%   while jb_max <= (numel(bin_quartile)-Nwindow)
%      jb=jb+1;
%      jb_max=jb+Nwindow;
%      if bin_quartile(jb) == bin_quartile(jb_max),
%          bin_index{jb} = find(all_Rg(:,2) == bin_quartile(jb));
%      else
%      bin_index{jb} = find(all_Rg(:,2) > bin_quartile(jb) & all_Rg(:,2) <= bin_quartile(jb_max));
%      end
%   end
%   binMax = numel(bin_index);   
%   for binN = 1:binMax  % finding the index in all_Rg where tau intensities satisfy the given conditions (i.e. finding the binning) 
%       bin_index{binN} = find(all_Rg > binN-1 & all_Rg(:,2) < binN+binSize); % this becomes the indices of all Rg values that have a corresponding tau intensity that falls within the bin range. The range of the bins goes like this: so first bin is 0-5, the 2nd bin is 1-6, etc.
%    end
  
binMax = nsamp;      
for binN = 1:binMax
%     if length(bin_index{binN}) == 0 % if there were no samples in a given bin, set things to NaN
%         Rg_mean_of_means(binN) = NaN;
%         Rg_means{binN} = NaN;
%         x_means_of_means(binN) = NaN;
%         x_means{binN}=NaN;
%     else
        drawn_samples = bstrp_func_loic_improved(1:nsamp,nsamp,iterations); %Emily modified to suit her purpose,no intensity binning 
        % If the number of data points in bin1 is 250 then this will pick 250 random numbers between 1 and 250, "iterations" times.
        % The output will be in [iterationsx250]. Notice that the first input is a
        % line vector, going from 1 to the length of the data
        % Linda's explanation: ...        
        all_Rg_bin = all_Rg;  % This makes a new matrix of Rg sample that only includes the Rg values within the specified bin range
% %         all_x_bin{binN}=all_Rg(bin_index{binN},2);
        for boot_ite = 1:iterations
            Rg_means(boot_ite) = mean(all_Rg(drawn_samples(boot_ite,:),1));  %Rg means are the means of each itinerated sample, so there will be 1000 means if there are 1000 itinerations.
%             x_means{binN}(boot_ite) = mean(all_x_bin(drawn_samples(boot_ite,:),1));   
        end
        Rg_mean_of_means = mean(Rg_means);  % Within the bin range, this is the mean of all 1000 means. The numel(Rg_mean_of_means) here will equal the number of bins.
%         x_mean_of_means(binMax) = mean(x_means{binMax});
%     end
    
end

    CI = prctile(Rg_means(1:1000),[alpha/2  100-alpha/2]);
    Rg_UB_of_means=Rg_mean_of_means(1)-CI(1,2);
    Rg_LB_of_means=CI(1,2)-Rg_mean_of_means(1);
    %errorbar(x_means_of_means(i),Rg_mean_of_means(i),Rg_mean_of_means(i)-CI(i,1),CI(i,2)-Rg_mean_of_means(i),'o','markersize',20,'linewidth',2,'CapSize',0,'color',[0 0 0.7])
    hold on
% 
% if zp==1
% figure
% for i = 1:binMax
%     %plot(i.*ones(length(Rg_means{i}),1),Rg_means{i},'o')
%     plot(Rg_means,'o')
%     hold on
% end

figure
plot(Rg_mean_of_means,'s','markersize',16,'linewidth',2);

% figure
% 
% errorbar(Rg_mean_of_means,Rg_LB_of_means,Rg_UB_of_means,'o','markersize',20,'linewidth',2,'CapSize',0)
% end

Rg_bstrap.bstrap_means=Rg_means;
Rg_bstrap.mean_of_means=Rg_mean_of_means;
Rg_bstrap.LB_of_means=Rg_LB_of_means;
Rg_bstrap.UB_of_means=Rg_UB_of_means;
% Rg_bstrap.x_mean_of_means=x_mean_of_means;

save('bstrap_means.mat','Rg_bstrap')
end
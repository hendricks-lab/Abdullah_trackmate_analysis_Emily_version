%% MSD Function

function[res]=msd_fun_v1(K,xk1,dtt,Ndat1,plt)%calculates x dimension MSD (1D)
%K=1:1:250;
Ndat=Ndat1;
for j=1:numel(K),
    k=K(j);
    dispsqk=[];
    for kd=1:Ndat,
        %clear xk
        xk=xk1;
        [mx,nx]=size(xk); % gives the number of columns then the number of rows (eg. 1x18 is [1 18])
        if mx<nx,
            xk=xk'; %switches the rows and columns
            [mx,nx]=size(xk);
        end
        if k<(mx-1),
            jmax=min(mx,max(K)*4);
            j0=1:(jmax-k);
            j1=(k+1):jmax;
            dkd = (xk(j1)-xk(j0)).^2; %displacement squared
            dispsqk=[dispsqk;dkd];
        end
    end
    jk=find(~isnan(dispsqk)); %taking out the NaNs
    N(j)=numel(dispsqk(jk)); %calculating the number of msds calculated
   msd(j)=sum(dispsqk(jk))/N(j); %calculates the 1d msd in the x axis I think (function is calculating the mean in this step)
   sem(j)=std(dispsqk(jk))./sqrt(N(j)); %calculates the standard error of the mean
end
%calculate diffusion coefficient
pmsd=polyfit(K,msd,1);%fitting a line to the first K frames of the msd
DiffCoef=pmsd(1)/2;

logtime=log10(dtt.*K); 
logmsd=log10(msd);
p=polyfit(logtime,logmsd,1);
tt1=logtime;
mm=polyval(p,tt1);

% figure(1), hold on
% plot(logtime, logmsd,'Color', 'k'), hold on
% plot(logtime, (p(1)*logtime+p(2)), 'Color', 'b')


if plt== 1
    cmap=colormap(jet(10)); % Color lookup table
    kcol=1; % Color code
    figure,hold on,
    plot(dtt.*K,msd,'.','LineWidth',2,'Color','k','LineWidth',2), xlabel('Time Interval (sec)');
    ylabel('Mean Squared Displacement (nm^{2})');
    errorbar(dtt.*K,msd,sem,'LineWidth',1);
    set(gca,'fontsize',16);
    publication_fig(0,0,1);

    axes('Position',[0.174385026737966 0.408450704225352 0.173743315508021 0.484859154929577]);
    plot(logtime,logmsd,'.','MarkerSize',15,'Color','k'), xlabel('log[Time Interval (sec)]'), ylabel('log[Mean Squared Displacement (nm^{2})]'),axis equal
    title(['m = ',num2str(p(1))]);
    hold on, plot(tt1,mm,'-','Color',cmap(kcol,:),'LineWidth',2)
else
end

timek=dtt.*K;
msd_f=msd;
log_tk=logtime;
log_msdf=logmsd;

res.slp=p(1);
res.logtk=log_tk;
res.log_msd=log_msdf;
res.timek=timek;
res.msd=msd_f;

clear xk; clear Ndat; clear k; clear K; clear kd; clear xk1; clear dkd; clear jk; clear N; clear j; clear msd; clear sem; clear pmsd; clear DiffCoeff;  clear logtime; clear logmsd;
clear p; clear tt1; clear mm; 
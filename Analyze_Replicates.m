%% analyze replicates from strr_Vhwidth_sweep

par.Vhwidth=0.175*linspace(0,1,4);
for w=1:4
    figure(w)
    hold on
    for r=1:1:size(strr,1)
        plot(strr{r,w}(:,1), strr{r,w}(:,2), 'k-')
    end
    xlabel('Iteration', 'fontsize',14)
ylabel('Average speed (distance/frame)', 'fontsize',14)
title(['s_{hw} = ',num2str(par.Vhwidth(w)/0.175), 's_{\mu}'],'fontsize', 14)
ylim([0,3.5])
end
%% get median trajectories
velR=[]; asR=[];
vdevR=[]; dmurR=[];
fr=strr{1,1}(:,1);
repdat=zeros(size(strr,1), size(strr,2), 1);
cm=cool(4);
% figure(1)

for w=1:4
    w
for r=1:size(strr,1)
    velR=[velR, strr{r,w}(:,2)];
    vdevR=[vdevR, strr{r,w}(:,3)];
    asR=[asR, strr{r,w}(:,4)];
    dmurR=[dmurR, strr{r,w}(:,6)];
    %get curve integrals--------
% repdat(r,w,1)=sum(strr{r,w}(:,2));
end
velR=median(velR,2);
    vdevR=median(vdevR,2);
    asR=median(asR,2);
    dmurR=median(dmurR,2);
    
plot(fr, velR,'-', 'linewidth',2,'color', cm(w,:)); 
% plot(fr, vdevR,'b-'); 
xlabel('Iteration', 'fontsize',14)
ylabel('Group avg speed (distance/frame)', 'fontsize',14)
% ylabel('mean AS', 'fontsize',14)
hold on
% title(['s_{hw} = ',num2str(par.Vhwidth(w)/0.175), 's_{\mu}'],'fontsize', 14)
end
std(repdat,0,1)
%% BOX PLOT
% plot(par.Vhwidth, mean(repdat,1), 'ko')
% errorbar(par.Vhwidth, mean(repdat,1), std(repdat,0,1), 'ko')
boxplot(repdat, 'Labels',{'0s_mu','0.33s_mu', '0.66s_mu','1.00s_mu'})
% xlim([-0.04, 0.2])
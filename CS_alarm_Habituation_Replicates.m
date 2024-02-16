%%
fps=30;
par.Tfinal=fps*150;%#number frames
Reps=2;
par.Np=60;%20:40:200;
Np=par.Np;
seq_data=1.8739; %mean speed over all ants pxl/frame. used as the equilib of the speed params
par.sd=0.001;%linspace(0,0.002,5);
par.Rcrit=40;%(in pxls) 40*mppx=5.8mm
par.a=0.1062;%linspace(0.01, 0.2, 5);
par.lamstar = linspace(0.5,3.5,3);
mppx=0.145177;%mm/pxl
par.kmax=6.8881;%seq_data*3.67;pxl/frame
strfr=cell(length(par.Np),1);
fdatr=cell(length(par.a),1);
clear fdatr strfr 
Vtracks=[];
VTRACKS_M=cell(Reps,1);
tic
for n=1:length(par.lamstar)
    n
    remaining_time(n,length(par.a),toc)
for r=1:Reps
    r
    [strf, fdat, ~]=CSAHab(par.Tfinal, par.Np, par.kmax, par.sd, par.Rcrit, par.a, par.lamstar(n));
   fdatr{n,1}(:,:,r)=fdat(1:end-1,:);
% strd=strf2strd(strf);
% Vtracks=[];
% for p=1:length(strd)
%     Vtracks=[Vtracks; strd{p}];
% end
%     VTRACKS_M{r,1}=Vtracks;
    
end
end
disp('done reps')

%% get mean trajectories for sweep over [param] - error bars?
whitebg([1 1 1])
fr=1:par.Tfinal;
    kerL=300;
    cm=copper(5);
%     cm=rand(5,3);
for n=1:length(par.lamstar)
    FDAT{n,1}=mean(fdatr{n,1},3);%average time series for all data across all replicates (vertically)
    ma_s=conv(FDAT{n,1}(:,3),ones(kerL,1)/kerL,'valid'); %smooth it
    plot((floor(kerL/2):par.Tfinal-floor(kerL/2))/fps, pf2ms(ma_s,mppx, fps), '-', 'linewidth',2,'color', cm(n,:)); hold on%plot 
%     plot(fr, FDAT{n,1}(:,3),'k-'); hold on %unsmoothed mean speed
%     plot(fr, FDAT{n,1}(:,4),'b-'); 
end
% legend({'\beta = 0', '\beta = 5E-4',  '\beta = 1E-3',  '\beta = 1.5E-3',  '\beta = 2E-3'}, 'Location', 'east','FontSize', 12, 'Box', 'off')
% legend({'N = 20', 'N = 60',  'N = 100',  'N = 140',  'N = 180'}, 'Location', 'northeast','FontSize', 12, 'Box', 'off')
legend({'a = 0.01', 'a = 0.058',  'a = 0.105 ',  'a = 0.153',  'a = 0.2'}, 'Location', 'northeast','FontSize', 12, 'Box', 'off')

xlabel('time (s)', 'FontSize', 15); ylabel('mm/s', 'FontSize', 15);
%% look at variation within replicates with all fixed params. plots smoothed version of mean speed trajectory for all replicates
for r=1:Reps
    ma_s=conv(fdatr{2,1}(:,3,r),ones(kerL,1)/kerL,'valid');
    plot((floor(kerL/2):par.Tfinal-floor(kerL/2))/fps, ma_s, '-', 'linewidth',1,'color', [0 0 0]); hold on
end
%% concat all replicate data
VTRACKS=[];
for r=1:length(VTRACKS_M)
    VTRACKS=cat(1,VTRACKS, VTRACKS_M{r,1});
end
%% collapsed violin plot
%replace speed param by actual speed for stats:
VTRACKS(:,5)=norm_array(VTRACKS(:,7:8));

vio_edg=2:900:par.Tfinal+2; %edges in frames for each segment, nonoverlapping
maxV=max(VTRACKS(:,5));
whitebg([1 1 1])
for t=1:length(vio_edg)-1
    t
    bloklog=(VTRACKS(:,1)>=vio_edg(t))&(VTRACKS(:,1)<vio_edg(t+1));
    [Vfreq, edges]=histcounts(pf2ms(VTRACKS(bloklog,7),mppx,fps), pf2ms(linspace(0,29,30),mppx,fps),'normalization','pdf');
    semilogy((edges(1:end-1)'), (Vfreq), '.-', 'color', [1-t/5 0 t/5],'linewidth', 1.5); hold on
end
ylim([1E-07, 1]); xlim([0, 140])
xlabel('speed (mm/s)', 'fontsize', 14); ylabel('probability', 'fontsize', 14)
legend('t \in [0, 30) s', 't \in [30, 60) s', 't \in [60, 90) s', 't \in [90, 120) s', 't \in [120, 150] s')

%%
for r=1:1:Reps
    
    figure(4)
    hold on
whitebg([1 1 1])
% plot(strr{r,1}(:,1), strr{r,1}(:,4), '-'); 

plot(strr{r,1}(:,1), strr{r,1}(:,2), '-');
% plot(strr{r,1}(:,1), strr{r,1}(:,4), 'r-');
% plot(strr{r,1}(:,4)*0.01, strr{r,1}(:,2),'-');
% ylim([0, 1])
xlabel('Iteration number', 'fontsize',14)
ylabel('Average Velocity (mm/s)', 'fontsize',14)
end
%% plot video
pause_time=.01;%plotting pause (s)
% circle
R=1040/2;%radius of nest. Nest area = circle
circl=[linspace(0, 2*pi, 80)', R*ones(80,1)];
[X, Y]=pol2cart(circl(:,1), circl(:,2));
cm=jet(200);
whitebg([0 0 0])
for t=550:4:par.Tfinal
    lst2=strf{t,1}(:,:);
hh=figure(1);

%       plot(lst2(:,2),lst2(:,3),'ko','MarkerSize',3,'MarkerFaceColor',[.5 .9 .5])
%       scatter(lst2(:,2),lst2(:,3),10,lst2(:,4),'filled');% color by speed
%             scatter(lst2(:,2),lst2(:,3),10,map2colsp(lst2(:,4),cm,[0,6]),'filled');%color by speed
  scatter(lst2(:,3),lst2(:,4),20,map2colsp(lst2(:,5),cm,[0 10]),'filled');%color = binary velocity
  
      hold on
      plot(X,Y,'g-')
% %       if isempty(I)==0
% % %       scatter(lst2(I(trans),2),lst2(I(trans),3),30,repmat([1 1 1],length(I(trans)),1))%plot one of the interacting pair circle
% %        scatter(lst2(I,2),lst2(I,3),30,repmat([1 1 1],length(I),1))%plot one of the interacting pair circle
% % 
% %       end
%       axis([-R R -R R])
      title(['Frame: ' num2str(t)], 'fontsize',13)
axis equal
% %       hold on
%         
% %         fig = gcf;
% %         fig.InvertHardcopy = 'off';
% %         saveas(hh, ['C:\Users\Michael_R_Lin\Documents\KANG-FEWELL LAB\Alarm Propagation\movies\' num2str(t) '.png'])
hold off

% figure(3); whitebg([1 1 1])%PHASE SPACE PLOT
% scatter(lst2(:,4),0,10,'filled'); 
% colormap cool
% scatter(lst2(:,4),yc, 10, 1.5*lst2(:,5)+1, 'filled')
% hold on
% axis([0 5 0 1])
% xlabel('\mu'); ylabel('AS'); zlabel('\Delta \mu')
% camorbit(t/Tfinal*10,0)

pause(pause_time)
end
%% plot colored tracks
strd=strf2strd(strf);
% 
whitebg([0 0 0])
cm=colormap(jet(200));
hold on
a=1;
for p=1:1:length(strd)
    co=rand(1,3);
plot(strd{p,1}(:,3),strd{p,1}(:,4),'-','Color',map2colsp(p,cm,[1 par.Np(1)]),'LineWidth', 1,'MarkerSize',2)%plot tracks
hold on
%  for f=2:length(strd{p,1}(:,1)) %plot velocity vectors
%     plot([strd{p,1}(f,3) strd{p,1}(f,3)+a*strd{p,1}(f,7)*cosd(strd{p,1}(f,6))], [strd{p,1}(f,4) strd{p,1}(f,4)+a*strd{p,1}(f,7)*sind(strd{p,1}(f,6))],'k-','LineWidth', 1);
%  end

axis equal
end
% hold on
% plot(contacts(contacts(:,1),4),contacts(contacts(:,1),5),'ks','MarkerSize',7,'MarkerFaceColor',[1 .8 0])
axis equal
%% plot lamda_t^p and s compare
% strd=strf2strd(strf);
int=2;
% plot(strd{1}(1:int:end,1), strd{1}(1:int:end,5), 'r-','linewidth', 2); hold on
for p=[50]
    strd{p}(:,9)= norm_array(strd{p}(:,7:8));
    ta=find(strd{p}(1:end,6)==1,1);
    plot(strd{p}(1:int:end,1), strd{p}(1:int:end,5), 'k-','linewidth', 2); hold on
    plot(strd{p}(ta:end,1), strd{p}(ta:end,5), 'r-','linewidth', 2);
    plot(strd{p}(1:int:ta,1), strd{p}(1:int:ta,9), 'ko','markersize', 2);
    plot(strd{p}(ta:int:end,1), strd{p}(ta:int:end,9), 'ro','markersize', 2);
end
xlabel('t', 'fontsize', 16)
ylabel('pxl/timestep', 'fontsize', 16)
legend('\lambda_t^p | A_t^p=0', '\lambda_t^p | A_t^p=1', 's_t^p | A_t^p=0',  's_t^p | A_t^p=1')
%%

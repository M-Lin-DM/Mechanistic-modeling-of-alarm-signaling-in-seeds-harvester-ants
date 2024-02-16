%% Get Alarm Replicates for param sweeps and P(s) distribution fitting
fps=30;
par.Tfinal=fps*150;%#number frames
Reps=2;
par.Np=[20, 60, 180];%linspace(60,180,3);%20:40:200;%
Np=par.Np;
seq_data=1.8739; %mean speed over all ants pxl/frame. used as the equilib of the speed params
par.sd=0.001;%linspace(0,0.002,3);%BETA
par.Rcrit=40;%(in pxls) 40*mppx=5.8mm
par.a=0.1;%linspace(0.01, 0.2, 5);
par.lamstar = 1.8739;%linspace(0.5,3.5,3);%
mppx=0.145177;%mm/pxl
par.kmax=6.8881;%seq_data*3.67;pxl/timestep
strfr=cell(length(par.Np),1);
fdatr=cell(length(par.a),1);
clear fdatr strfr 
Vtracks=[];
VTRACKS_n=cell(Reps,1);
vio_edg=2:900:par.Tfinal+2; %edges in frames for each segment, nonoverlapping
Pos = cell(5,3); %we use 5 time windows and 3 param values for each param sweep
%%
tic
for n=1:3 %each parameter value
    n
    remaining_time(n,length(par.a),toc)
    for r=1:Reps
        r
       [strf, fdat, ~]=CSAHab(par.Tfinal, par.Np, par.kmax, par.sd(n), par.Rcrit, par.a, par.lamstar);
       fdatr{n,1}(:,:,r)=fdat(1:end-1,:);
    strd=strf2strd(strf);
    Vtracks=[];
    for p=1:length(strd)
        Vtracks=[Vtracks; strd{p}(:,[1,7:8])]; %grab only speed data for memory saving
    end
        VTRACKS_n{r,1}=Vtracks;

    end%-------------------------------end all replicates--------------
% concat all replicate data
VTRACKS=[];
    for r=1:length(VTRACKS_n)
        VTRACKS=cat(1,VTRACKS, VTRACKS_n{r,1});
    end
    VTRACKS(:,4)=norm_array(VTRACKS(:,2:3));
    
    for t=1:length(vio_edg)-1
    t
    bloklog=(VTRACKS(:,1)>=vio_edg(t))&(VTRACKS(:,1)<vio_edg(t+1));
    [Vfreq, edges]=histcounts(pf2ms(VTRACKS(bloklog,4),mppx,fps), pf2ms(linspace(0,29,30),mppx,fps),'normalization','pdf');
    Pos{t,n} = [edges(1:end-1)', Vfreq'];
    end
    remaining_time(n,3,toc)
end
toc
disp('done reps')
%% fit P(s) for all param indicies n and windows w
bins=(Pos{1,1}(:,1) + [Pos{1,1}(2:end,1); Pos{1,1}(end,1)+4.35531])/2; %get bin centers
for n=1:3
    for w = 1:5
        logp = [bins, log(Pos{w,n}(:,2))];
        logp = logp(~isnan(logp(:,2)), :); logp = logp(~isinf(logp(:,2)), :); %ignore nan and inf
        [m, b, sm, sb] = linearfit(logp(:,1), logp(:,2));
        fitpar{w,n} = [m, b, sm, sb];
    end
end
%% plot
vio_bins = (vio_edg(1:end-1) + vio_edg(2:end))/2;
cm=copper(3);
for n=1:3
    par_n = cat(1, fitpar{:,n});
%     plot(vio_bins, par_n(:,1), 'o-', 'color', cm(n,:))
    errorbar(vio_bins/fps, par_n(:,1), 2*par_n(:,3), 'o-', 'color', cm(n,:), 'linewidth', 2); hold on
end
% ylim([1E-07, 1]); xlim([0, 140])
xlabel('time (s)', 'fontsize', 14); ylabel('\gamma', 'fontsize', 16)
% legend({'N = 60',  'N = 120',  'N = 180'}, 'Location', 'northeast','FontSize', 12, 'Box', 'off')
% legend({'\lambda^* = 0.5',  '\lambda^* = 2',  '\lambda^* = 3.5'}, 'Location', 'northeast','FontSize', 12, 'Box', 'off')
legend({'\beta = 0',  '\beta = 0.001',  '\beta = 0.002'}, 'Location', 'southwest','FontSize', 12, 'Box', 'off')

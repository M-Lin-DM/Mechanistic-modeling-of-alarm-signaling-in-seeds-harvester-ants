%% Get Alarm Replicates for param sweeps and P(s) distribution fitting
fps=30;
Tinject = 2; %= 4500/6
clear par
par.Tfinal=4500;%fps*150 + Tinject;%#number frames 4800=2m 40s = empirical video time. we are using 4500 because we can evenly divide it into 5 segments of 900fr (30s) 
Reps=300
par.Np=[60];%linspace(60,180,3);%20:40:200;%
% Np=par.Np;
fps=30;%frames per second for figures
mppx=0.145177;%mm/pxl
%CHANGE these if running replicates corresponding to 'get_CSA_output' (one
%param setting per colony)
% Colony b:
%0.00022, 60, ms2pf(2,mppx, fps), 0.3, 0.1
par.seq_data = 1.87*1%0.5513;% %LAMBDA
par.Rcrit=58.2;%(in pxls) 40*mppx=5.8mm
par.a=0.1497%0.1%;%linspace(0.01, 0.2, 5);
par.sd=0.00032;%linspace(2E-4,0.002,3);%0.001%BETA

par.anaught = 0.1;%[0.01, 0.06, 0.11];

% seq_data=ms2pf(5,mppx, fps); %mean speed over all ants pxl/frame. used as the equilib of the speed params
mppx=0.145177;%mm/pxl
par.kmax=6.8881*1;%pxl/timestep the value is necessary to max the max 30fps which is the last bin in the s_avg data empirical
% strfr=cell(length(par.Np),1);
% fdatr=cell(length(par.a),1);
% clear fdatr strfr 

window_edge=2:900:par.Tfinal+2; %edges in frames for each segment, nonoverlapping
s_edge = linspace(0,29,31);% global bins for s, for all reps and params
% Pos = cell(5,3); %we use 5 time windows and 3 param values for each param sweep

%% run parameter sweep
tic
for n=1:1% %each parameter value
%     n
    VFREQ =zeros(length(window_edge)-1, length(s_edge) -1, Reps); %P(s) in each window for one param level
    GMS = zeros(Reps, par.Tfinal);%mean(s) over time for one param level, all replicates
%     remaining_time(n,length(par.Np),toc)
    parfor r=1:Reps
        r
        disp(['parameter level: ' num2str(n) ' rep #: ' num2str(r) ' Remaining time in parameter level: '])
%         remaining_time(r,Reps,toc)
       [strf, fdat]=CSAHab_fast(par.Tfinal, par.Np, par.kmax, par.sd, par.Rcrit, par.a, par.seq_data, par.anaught);
        VFREQ(:,:,r) = strf2Vfreq(strf, window_edge, s_edge);
        GMS(r,:)=fdat(:,3)'; %framewise mean speed
        lambda_bar(r,:) = fdat(:,5)'; %framewise mean lambda
        Acount(r,:) = fdat(:,2)'; %total # alarmed

    end%-------------------------------end all replicates--------------
    save(['D:\Dissertation work ASU\Alarm Propagation\Habituation Model Replicates\N_' num2str(par.Np(n)) '_using_fitted_params300.mat'], 'GMS', 'par', 'lambda_bar', 'Acount', 'VFREQ')
%     save(['D:\KANG LAB\Alarm Propagation\Habituation Model Replicates\Colony_cC.mat'], 'GMS', 'par')

end
toc
disp('done reps')
%% LOAD DATA , pool data for all replicates. fit P(s) over time to get gammas. for all param indicies n and windows w
par.Tfinal=4500;
window_edge=2:900:par.Tfinal+2; %edges in frames for each segment, nonoverlapping
s_edge = linspace(0,29,31);% global bins for s, for all reps and params
fps=30;
segment = 1:Reps; %number of reps to sample
for n=1
    %MAKE SURE THE FILE YOU ARE LOADING IS CORRECT
    load(['D:\Dissertation work ASU\Alarm Propagation\Habituation Model Replicates\N_' num2str(par.Np(n)) '_using_fitted_params300.mat'])%MAKE SURE THE FILE YOU ARE LOADING IS CORRECT
    [out, Vfreq_pdf, bins, Vfreq_se]=VFREQ2gamma_of_t(VFREQ(:,:,segment), s_edge);%collapse all reps and fit each of the windows' histogram to a line, then get gamma for each line
    Vfreqpdf{n}= Vfreq_pdf; %VFREQ2gamma_of_t has already normalized this to be a pdf AND units of bins are already in mm/s
    Vfreqse{n} = Vfreq_se;
%     par_n=out;
    GMS_reps{n}=[mean(GMS(segment,:), 1)', std(GMS(segment,:), 0, 1)'/sqrt(length(segment))]; %mean, standard error over reps
    lambda_bar_reps{n} = [mean(lambda_bar(segment,:), 1)', std(lambda_bar(segment,:), 0, 1)'/sqrt(length(segment))];
    Acount_reps{n} = [mean(Acount(segment,:), 1)', std(Acount(segment,:), 0, 1)'/sqrt(length(segment)), quantile(Acount(segment,:), [0.25, 0.5, 0.75])'];
end
%% semilog MODEL vio

window_edge=2:900:par.Tfinal+2; 
%replace speed param by actual speed for stats:
window_bins = (window_edge(1:end-1) + window_edge(2:end))/2;
whitebg([1 1 1])
n=1;
figure(1)
for t=1:1:length(window_bins)
    semilogy(bins, (Vfreqpdf{n}(t,:)), '.-', 'color', [1-t/5 0 t/5],'linewidth', 1); hold on
end
% for t=1:1:length(window_bins)
%         fill([bins'; flipud(bins')], log10([Vfreqpdf{n}(t,:)' - Vfreqse{n}(t,:)'; flipud(Vfreqpdf{n}(t,:)' + Vfreqse{n}(t,:)')]), [1-t/5 0 t/5], 'linestyle','none');  hold on %"BANDZ" error bars for ultimate style
% end
% ylim([1E-08, 1]); 
xlim([0, 50])
xlabel('speed (mm/s)', 'fontsize', 14); ylabel('probability', 'fontsize', 14)
axis square
% legend('t \in [0, 30) s', 't \in [30, 60) s', 't \in [60, 90) s', 't \in [90, 120) s', 't \in [120, 150] s')
%you can put a command for axis labels
% title(['N = ' num2str(par.Np(n)) ' | 300 replicates'])
alpha(0.4)
%% linear vio plot MODEL
figure(2)
subplot(1,2,2)
for t=1:1:length(window_bins)
    plot(bins, (Vfreqpdf{n}(t,:)), '-', 'color', [1-t/5 0 t/5],'linewidth', 1); hold on
end
for t=1:1:length(window_bins)
        fill([bins'; flipud(bins')], ([Vfreqpdf{n}(t,:)' - Vfreqse{n}(t,:)'; flipud(Vfreqpdf{n}(t,:)' + Vfreqse{n}(t,:)')]), [1-t/5 0 t/5], 'linestyle','none');  hold on %"BANDZ" error bars for ultimate style
end
ylim([0, 0.21]); xlim([0, 130])
xlabel('speed (mm/s)', 'fontsize', 14); ylabel('probability', 'fontsize', 14)
legend('t \in [0, 30) s', 't \in [30, 60) s', 't \in [60, 90) s', 't \in [90, 120) s', 't \in [120, 150] s')
%% for empirical plot, make vfreq over windows for all T vid

strf=strf_IT4;
Vfreq_over_windows_emp(:,:,1) = strf2Vfreq(strf, window_edge, s_edge);
strf=strf_IT45;
Vfreq_over_windows_emp(:,:,2) = strf2Vfreq(strf, window_edge, s_edge);
strf=strf_IT71;
Vfreq_over_windows_emp(:,:,3) = strf2Vfreq(strf, window_edge, s_edge);


[~, Vfreq_pdf_emp, bins, Vfreq_se_emp]=VFREQ2gamma_of_t(Vfreq_over_windows_emp, s_edge);
%% semilog plot empirical
figure(4)
for t=1:1:length(window_bins)
    semilogy(bins, (Vfreq_pdf_emp(t,:)), '.-', 'color', [1-t/5 0 t/5],'linewidth', 1); hold on
end
for t=1:1:length(window_bins)
        fill([bins'; flipud(bins')], log10([Vfreqpdf{n}(t,:)' - Vfreqse{n}(t,:)'; flipud(Vfreqpdf{n}(t,:)' + Vfreqse{n}(t,:)')]), [1-t/5 0 t/5], 'linestyle','none');  hold on %"BANDZ" error bars for ultimate style
end
xlim([0, 50])
xlabel('speed (mm/s)', 'fontsize', 14); ylabel('probability', 'fontsize', 14)
axis square
% legend('t \in [0, 30) s', 't \in [30, 60) s', 't \in [60, 90) s', 't \in [90, 120) s', 't \in [120, 150] s')
%you can put a command for axis labels
% title(['N = ' num2str(par.Np(n)) ' | 300 replicates'])
alpha(0.1)
%% linear vio plot EMPirical
figure(2)
subplot(1,2,1)
for t=1:1:length(window_bins)
    plot(bins, (Vfreq_pdf_emp(t,:)), '-', 'color', [1-t/5 0 t/5],'linewidth', 1); hold on
end
for t=1:1:length(window_bins)
        fill([bins'; flipud(bins')], ([Vfreq_pdf_emp(t,:)' - Vfreq_se_emp(t,:)'; flipud(Vfreq_pdf_emp(t,:)' + Vfreq_se_emp(t,:)')]), [1-t/5 0 t/5], 'linestyle','none');  hold on %"BANDZ" error bars for ultimate style
end
alpha(0.05)
ylim([0, 0.21]); xlim([0, 130])
xlabel('speed (mm/s)', 'fontsize', 14); ylabel('probability', 'fontsize', 14)
legend('t \in [0, 30) s', 't \in [30, 60) s', 't \in [60, 90) s', 't \in [90, 120) s', 't \in [120, 150] s')

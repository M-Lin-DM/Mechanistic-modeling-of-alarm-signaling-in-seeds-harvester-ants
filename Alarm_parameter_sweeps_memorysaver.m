%% Get Alarm Replicates for param sweeps and P(s) distribution fitting
fps=30;
Tinject = 2; %= 4500/6
clear par
par.Tfinal=5100;%fps*150 + Tinject;%#number frames
Reps=300;
par.Np=[20, 60, 180];%linspace(60,180,3);%20:40:200;%
% Np=par.Np;
fps=30;%frames per second for figures
mppx=0.145177;%mm/pxl
%CHANGE these if running replicates corresponding to 'get_CSA_output' (one
%param setting per colony)
% Colony b:
%0.00022, 60, ms2pf(2,mppx, fps), 0.3, 0.1
par.sd=0.00022;%linspace(2E-4,0.002,3);%0.001%BETA
par.Rcrit=60;%(in pxls) 40*mppx=5.8mm
par.seq_data = ms2pf(2,mppx, fps);%
par.a=0.3;%linspace(0.01, 0.2, 5);
par.anaught = 0.0;%[0.01, 0.06, 0.11];

% seq_data=ms2pf(5,mppx, fps); %mean speed over all ants pxl/frame. used as the equilib of the speed params
mppx=0.145177;%mm/pxl
par.kmax=6.8881;%pxl/timestep the value is necessary to max the max 30fps which is the last bin in the s_avg data empirical
% strfr=cell(length(par.Np),1);
% fdatr=cell(length(par.a),1);
% clear fdatr strfr 

% window_edge=2:900:par.Tfinal+2; %edges in frames for each segment, nonoverlapping
% s_edge = linspace(0,29,31);% global bins for s, for all reps and params
% Pos = cell(5,3); %we use 5 time windows and 3 param values for each param sweep

%% run parameter sweep
tic
for n=1:1% %each parameter value
%     n
%     VFREQ =zeros(length(window_edge)-1, length(s_edge) -1, Reps); %P(s) in each window for one param level
    GMS = zeros(Reps, par.Tfinal);%mean(s) over time for one param level, all replicates
%     remaining_time(n,length(par.Np),toc)
    parfor r=1:Reps
        r
        disp(['parameter level: ' num2str(n) ' rep #: ' num2str(r) ' Remaining time in parameter level: '])
%         remaining_time(r,Reps,toc)
       [strf, fdat]=CSAHab_fast(par.Tfinal, par.Np, par.kmax, par.sd, par.Rcrit, par.a, par.seq_data, par.anaught);
%         VFREQ(:,:,r) = strf2Vfreq(strf, window_edge, s_edge);
        GMS(r,:)=fdat(:,3)'; %framewise mean speed
%         lambda_bar(r,:) = fdat(:,5)'; %framewise mean lambda
%         Acount(r,:) = fdat(:,2)'; %total # alarmed

    end%-------------------------------end all replicates--------------
    save(['D:\KANG LAB\Alarm Propagation\Habituation Model Replicates\N_' num2str(par.Np(n)) '.mat'], 'GMS', 'par', 'lambda_bar', 'Acount')
%     save(['D:\KANG LAB\Alarm Propagation\Habituation Model Replicates\Colony_cC.mat'], 'GMS', 'par')

end
toc
disp('done reps')
%% LOAD DATA , pool data for all replicates. fit P(s) over time to get gammas. for all param indicies n and windows w
% window_edge=2:900:par.Tfinal+2; %edges in frames for each segment, nonoverlapping
% s_edge = linspace(0,29,31);% global bins for s, for all reps and params
fps=30;
segment = 1:100%300; %number of reps to sample
for n=1:3
    load(['D:\Dissertation work ASU\Alarm Propagation\Habituation Model Replicates\N_' num2str(par.Np(n)) '.mat'])
%     [out, Vfreq_pdf, bins, Vfreq_se]=VFREQ2gamma_of_t(VFREQ(:,:,segment), s_edge);%collapse all reps and fit each of the windows' histogram to a line, then get gamma for each line
    Vfreqpdf{n}= Vfreq_pdf;
    Vfreqse{n} = Vfreq_se;
%     par_n=out;
    GMS_reps{n}=[mean(GMS(segment,:), 1)', std(GMS(segment,:), 0, 1)'/sqrt(length(segment))]; %mean, standard error over reps
    lambda_bar_reps{n} = [mean(lambda_bar(segment,:), 1)', std(lambda_bar(segment,:), 0, 1)'/sqrt(length(segment))];
    Acount_reps{n} = [mean(Acount(segment,:), 1)', std(Acount(segment,:), 0, 1)'/sqrt(length(segment)), quantile(Acount(segment,:), [0.25, 0.5, 0.75])'];
end
%% store the colony ABC ma_s and ma_se
kerL=100;
kernel = ones(kerL,1)/kerL;
ma_s_cC=conv(GMS_reps{n}(:,1), kernel,'valid'); %smooth it, convert to mm/s
ma_se_cC=conv(GMS_reps{n}(:,2), kernel,'valid'); %smooth it, convert to mm/s
% fr = (1:length(ma_s_A))'/fps +offset/fps;
% fr = ((1:length(ma_s)-clip)+offsetb)/fps;
fr=1:length(ma_s_cC);
fr= fr';
%%
n=1
fr=1:length(GMS_reps{n});
errorbar(fr, GMS_reps{n}(:,1),GMS_reps{n}(:,2), '-', 'linewidth',1);
%% MAIN parameter sweep plots
% get mean trajectories for sweep over [param] - error bars?
whitebg([1 1 1])
kerL=100;
kernel = ones(kerL,1)/kerL;
fr=(floor(kerL/2):par.Tfinal-floor(kerL/2))/fps; %fr for "valid" conv
fr=fr';
% fr=1:par.Tfinal;
    cm=copper(3);
%     cm=rand(5,3);
for n=1:3
    %plot mean s
    ma_s=pf2ms(conv(GMS_reps{n}(:,1), kernel,'valid'), mppx, fps); %smooth it, convert to mm/s
    ma_se=pf2ms(conv(GMS_reps{n}(:,2), kernel,'valid'), mppx, fps); %smooth it, convert to mm/s
    %plot mean lambda:
%     ma_s=pf2ms(conv(lambda_bar_reps{n}(:,1), kernel,'valid'), mppx, fps); %smooth it, convert to mm/s
%     ma_se=pf2ms(conv(lambda_bar_reps{n}(:,2), kernel,'valid'), mppx, fps); %smooth it, convert to mm/s
    %plot total alarmed 
%         ma_s=conv(Acount_reps{n}(:,1), kernel,'valid')/par.Np(n); %smooth it
%         ma_se=conv(Acount_reps{n}(:,2), kernel,'valid')/par.Np(n); %smooth it

    switch n
        case 1
                    plot(fr, ma_s, '-', 'linewidth',2,'color', cm(n,:)); hold on
        case 2
                    plot(fr, ma_s, '--', 'linewidth',2,'color', cm(n,:)); hold on
        case 3
                    plot(fr, ma_s, ':', 'linewidth',2,'color', cm(n,:)); hold on
    end
%     plot(fr, ma_s+ma_se, '-', 'linewidth',2,'color', cm(n,:)); hold on%plot 

%get peak heights, times
disp('n = ')
n
max(ma_s) - pf2ms(par.lamstar, mppx, fps)
argmax(ma_s) - 750
%get integral
sum(ma_s(750:end)*(1/fps))
sum(ma_se(750:end)*(1/fps))

end %compute smoothing of mean curve

% for n=1:3%compute smoothing of se curve
%         ma_s=pf2ms(conv(GMS_reps{n}(:,1), kernel,'valid'), mppx, fps); %smooth it, convert to mm/s
%         ma_se=pf2ms(conv(GMS_reps{n}(:,2), kernel,'valid'), mppx, fps); %smooth it, convert to mm/s
%         %plot mean lambda:
% %         ma_s=pf2ms(conv(lambda_bar_reps{n}(:,1), kernel,'valid'), mppx, fps); %smooth it, convert to mm/s
% %         ma_se=pf2ms(conv(lambda_bar_reps{n}(:,2), kernel,'valid'), mppx, fps); %smooth it, convert to mm/s
%         %plot fraction alarmed
% %         ma_s=conv(Acount_reps{n}(:,1), kernel,'valid')/par.Np(n); %smooth it, convert to mm/s
% %         ma_se=conv(Acount_reps{n}(:,2), kernel,'valid')/par.Np(n); %smooth it, convert to mm/s
%         
%         fill([fr; flipud(fr)], [ma_s-ma_se; flipud(ma_s+ma_se)], cm(n,:), 'linestyle','none');  hold on %"BANDZ" error bars for ultimate style
% end
% legend({'N = 20', 'N = 60',  'N = 180'}, 'Location', 'northeast','FontSize', 12, 'Box', 'off')
% legend({'$\rho = 0.01$', '$\rho = 0.06$',  '$\rho = 0.11$'}, 'Location', 'northeast','FontSize', 13, 'Box', 'off', 'Interpreter','latex')
% legend({'$\beta = 2\times 10^{-4}$', '$\beta = 1.1\times 10^{-3}$',  '$\beta = 2 \times 10^{-3}$'}, 'Location', 'northeast','FontSize', 13, 'Box', 'off', 'Interpreter','latex')
% legend({'$a = 0.01$', '$a = 0.058$',  '$a = 0.105$ '}, 'Location', 'northeast','FontSize', 13, 'Box', 'off', 'Interpreter','latex')
% legend({'$\lambda^* = 2.2$', '$\lambda^* = 8.7$',  '$\lambda^* = 15.2$'}, 'Location', 'northeast','FontSize', 13, 'Box', 'off', 'Interpreter','latex')
% legend({'$r = 20$', '$r = 40$',  '$r = 60$'}, 'Location', 'northeast','FontSize', 13, 'Box', 'off', 'Interpreter','latex')


alpha(0.4)
xlabel('time (s)', 'FontSize', 15); 
ylim([6, 18])
ylabel('$\overline{s_t}$ (mm/s)', 'FontSize', 19 , 'Interpreter','latex');
% ylabel('$\overline{\lambda_t}$ (mm/s)', 'FontSize', 19 , 'Interpreter','latex');
% ylabel('fraction of colony alarmed (mm/s)', 'FontSize', 19 , 'Interpreter','latex');
% title(['300 replicates'], 'Interpreter','latex')
%% plot gamma over windows
window_edge=2:900:par.Tfinal+2; 
vio_bins = (window_edge(1:end-1) + window_edge(2:end))/2;
cm=copper(3);
for n=1:3
%     par_n = cat(1, fitpar{:,n});
%     plot(vio_bins, par_n(:,1), 'o-', 'color', cm(n,:))
        switch n
        case 1
                    errorbar(vio_bins/fps, fitpar{n}(:,1), fitpar{n}(:,3), '-', 'color', cm(n,:), 'linewidth', 2); hold on
        case 2
                    errorbar(vio_bins/fps, fitpar{n}(:,1), fitpar{n}(:,3), '--', 'color', cm(n,:), 'linewidth', 2); hold on
        case 3
                    errorbar(vio_bins/fps, fitpar{n}(:,1), fitpar{n}(:,3), ':', 'color', cm(n,:), 'linewidth', 2); hold on
        end
end
errorbar(vio_bins/fps, fitpar_emp(:,1), fitpar_emp(:,3), 'bs-', 'linewidth', 2); %plot empirical gamma
% ylim([1E-07, 1]); xlim([0, 140])
xlabel('time (s)', 'fontsize', 14); ylabel('\gamma', 'fontsize', 16)
legend({'N = 20',  'N = 60',  'N = 180', 'N = 60 (empirical)'}, 'Location', 'northeast','FontSize', 12, 'Box', 'off')
% legend({'\lambda^* = 0.5',  '\lambda^* = 2',  '\lambda^* = 3.5'}, 'Location', 'northeast','FontSize', 12, 'Box', 'off')
% legend({'\beta = 0',  '\beta = 0.001',  '\beta = 0.002'}, 'Location', 'southwest','FontSize', 12, 'Box', 'off')
% title(['300 replicates'])
%% check fits
% plot(bins, log(Vfreqpdf{3}(1,:)), '-'); hold on
% plot(bins, -3.11-0.0583*bins, 'r-')
plot(bins, (Vfreqpdf{3}(1,:)), '-'); hold on
plot(bins, exp(-3.11-0.0583*bins), 'r-')
%% updated violin plot for many replicates
window_edge=2:900:par.Tfinal+2; 
%replace speed param by actual speed for stats:
vio_bins = (window_edge(1:end-1) + window_edge(2:end))/2;
whitebg([1 1 1])
n=3;

for t=1:1:length(vio_bins)
    plot(bins, log10(Vfreqpdf{n}(t,:)), '.-', 'color', [1-t/5 0 t/5],'linewidth', 1); hold on
end
for t=1:1:length(vio_bins)
        fill([bins'; flipud(bins')], log10([Vfreqpdf{n}(t,:)' - Vfreqse{n}(t,:)'; flipud(Vfreqpdf{n}(t,:)' + Vfreqse{n}(t,:)')]), [1-t/5 0 t/5], 'linestyle','none');  hold on %"BANDZ" error bars for ultimate style
end
% ylim([1E-08, 1]); xlim([0, 130])
xlabel('speed (mm/s)', 'fontsize', 14); ylabel('probability', 'fontsize', 14)
legend('t \in [0, 30) s', 't \in [30, 60) s', 't \in [60, 90) s', 't \in [90, 120) s', 't \in [120, 150] s')
%you can put a command for axis labels
title(['N = ' num2str(par.Np(n)) ' | 300 replicates'])
alpha(0.4)
%%
plot(bt,vt,'.-')
%% NO SMOOTHING MAIN parameter sweep plots 
% get mean trajectories for sweep over [param] - error bars?
whitebg([1 1 1])
% kerL=100;
% kernel = ones(kerL,1)/kerL;
% fr=(floor(kerL/2):par.Tfinal-floor(kerL/2))/fps; %fr for "valid" conv
fr=1:par.Tfinal;
fr=fr'/fps;
% fr=1:par.Tfinal;
    cm=copper(4);
%     cm=rand(5,3);
% cm=copper(3);
for n=1:3%fliplr(1:3)% %
    %plot mean s
%     ma_s=GMS_reps{n}(:,1); %smooth it, convert to mm/s
%     ma_se=GMS_reps{n}(:,2); %smooth it, convert to mm/s
    %plot mean lambda:
%     ma_s=lambda_bar_reps{n}(:,1); %smooth it, convert to mm/s
%     ma_se=lambda_bar_reps{n}(:,2); %smooth it, convert to mm/s
    %plot total alarmed 
        ma=Acount_reps{n}/par.Np(n); %smooth it

%     switch n
%         case 1
%             plot_speed_error_bands(fr*fps, ma_s, ma_se, 0, 0, 10, cm(n,:)); hold on
%             l1=plot(fr, pf2ms(ma_s, mppx, fps), '-', 'linewidth',1,'color', cm(n,:)); 
%         case 2
%             plot_speed_error_bands(fr*fps, ma_s, ma_se, 0, 0, 10, cm(n,:)); hold on
%             l2=plot(fr, pf2ms(ma_s, mppx, fps), '+', 'markersize',4,'color', cm(n,:)); 
%         case 3
%             plot_speed_error_bands(fr*fps, ma_s, ma_se, 0, 0, 10, cm(n,:)); hold on
%             l3=plot(fr, pf2ms(ma_s, mppx, fps), 'o', 'markersize',4,'color', cm(n,:)); 
%     end
    %for Acount: quantiles error bands
    switch n
        case 1
            plot_speed_error_bands_quantiles(fr*fps, ma(:,3), ma(:,5), 0, 0, cm(n,:)); hold on
            l1=plot(fr, ma(:,4), '-', 'linewidth',1.5,'color', cm(n,:)); 
        case 2
            plot_speed_error_bands_quantiles(fr*fps, ma(:,3), ma(:,5), 0, 0, cm(n,:)); hold on
            l2=plot(fr, ma(:,4), '--', 'linewidth',1.5,'color', cm(n,:)); 
        case 3
            plot_speed_error_bands_quantiles(fr*fps, ma(:,3), ma(:,5), 0, 0, cm(n,:)); hold on
            l3=plot(fr, ma(:,4), ':', 'linewidth',1.5,'color', cm(n,:)); 
    end

%get peak heights, times
disp('n = ')
n
% max(ma_s) - pf2ms(par.lamstar, mppx, fps)
% argmax(ma_s) - 750
%get integral
% sum(ma_s(750:end)*(1/fps))
% sum(ma_se(750:end)*(1/fps))

end %compute smoothing of mean curve

legend([l1 l2 l3], {'N = 20', 'N = 60',  'N = 180'}, 'Location', 'southeast','FontSize', 13, 'Box', 'off')
% legend([l1 l2 l3], {'$\rho = 0.01$', '$\rho = 0.06$',  '$\rho = 0.11$'}, 'Location', 'northeast','FontSize', 14, 'Box', 'off', 'Interpreter','latex')
% legend([l1 l2 l3], {'$\beta = 2\times 10^{-4}$', '$\beta = 1.1\times 10^{-3}$',  '$\beta = 2 \times 10^{-3}$'}, 'Location', 'northeast','FontSize', 14, 'Box', 'off', 'Interpreter','latex')
% legend([l1 l2 l3], {'$a = 0.01$', '$a = 0.058$',  '$a = 0.105$ '}, 'Location', 'northeast','FontSize', 14, 'Box', 'off', 'Interpreter','latex')
% legend([l1 l2 l3], {'$\lambda^* = 2.2$', '$\lambda^* = 8.7$',  '$\lambda^* = 15.2$'}, 'Location', 'northeast','FontSize', 14, 'Box', 'off', 'Interpreter','latex')
% legend([l1 l2 l3], {'$r = 20$', '$r = 40$',  '$r = 60$'}, 'Location', 'northeast','FontSize', 14, 'Box', 'off', 'Interpreter','latex')


alpha(0.25)
xlabel('time (s)', 'FontSize', 18, 'Interpreter','latex'); 
% ylim([6, 20])
ylim([0, 1.15])
% ylabel('$\overline{s_t}$ (mm/s)', 'FontSize', 18 , 'Interpreter','latex');
% ylabel('$\overline{\lambda_t}$ (mm/s)', 'FontSize', 18 , 'Interpreter','latex');
ylabel('fraction of colony alarmed', 'FontSize', 18 , 'Interpreter','latex');
% title(['300 replicates'], 'Interpreter','latex')
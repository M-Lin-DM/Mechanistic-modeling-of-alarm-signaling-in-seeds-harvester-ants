%% validation plot fitting with replicates data
fps=30;%frames per second for figures
mppx=0.145177;%mm/pxl
offset = 5071;
linewidth = 1;
int=200;
markersize = 3;
clip = 350;

cm=copper(3);
%treatment simulations
% [ma_s_A, ma_se_A, pdat_A, fdat_A] = get_CSA_output(0.001, 40, ms2pf(5,mppx, fps), 0.1, 0.1);% sd, Rcrit, seq_data, a, a0
% [ma_s_B, ma_se_B, pdat_B, fdat_B] = get_CSA_output(0.00025, 60, ms2pf(4,mppx, fps), 0.06, 0.1);% sd, Rcrit, seq_data, a, a0
% [ma_s_C, ma_se_C, pdat_C, fdat_C] = get_CSA_output(0.00022, 60, ms2pf(2,mppx, fps), 0.3, 0.1);% sd, Rcrit, seq_data, a, a0

% control simulations
% [ma_s_cA, ma_se_cA] = get_CSA_output(0.001, 40, ms2pf(5,mppx, fps), 0.1, 0);% sd, Rcrit, seq_data, a, a0
% [ma_s_cB, ma_se_cB] = get_CSA_output(0.00025, 60, ms2pf(4,mppx, fps), 0.06, 0);% sd, Rcrit, seq_data, a, a0
% [ma_s_cC, ma_se_cC] = get_CSA_output(0.00022, 40, ms2pf(2,mppx, fps), 0.25, 0);% sd, Rcrit, seq_data, a, a0


subplot(3,1,1)
plot_speed(ma_s_CA,'b-.', 0,2,1,1, 0) %empirical plots dont modify
plot_speed(ma_s_TA,'b-', 5071,2,1,1, 0)
% 
plot_speed(ma_s_A,'r-o', offset, linewidth, int, markersize, clip)
plot_speed_error_bands(fr, ma_s_A, ma_se_A, offset, clip, 2^.5)
w=2.5;
plot_speed(w*ma_s_cC,'r-o', 0, linewidth, int, markersize, 0)
plot_speed_error_bands(fr, w*ma_s_cC, ma_se_cC, 0, 0, 1)
ylim([0, 17])
xlim([0 330])
xlabel('time (seconds)', 'fontsize', 17, 'interpreter', 'latex')
ylabel('$\overline{s_t}$ (mm/s)',  'fontsize', 18, 'interpreter', 'latex')
legend('$\overline{s_t}$ pre-introduction', '$\overline{s_t}$ post-introduction', '$\overline{s_t}$ model', 'Location', 'northwest','FontSize', 14, 'Box', 'off', 'interpreter', 'latex')


subplot(3,1,2)
plot_speed(ma_s_CB,'b-.', 0,2,1,1, 0)
plot_speed(ma_s_TB,'b-', 5071,2,1,1, 0)

plot_speed(ma_s_B,'r-o', offset, linewidth, int, markersize, clip)
plot_speed_error_bands(fr, ma_s_B, ma_se_B, offset, clip, 1)

w=2;
plot_speed(w*ma_s_cC,'r-o', 0, linewidth, int, markersize, 0)
plot_speed_error_bands(fr, w*ma_s_cC, ma_se_cC, 0, 0, 1)

ylim([0, 17])
xlim([0 330])
xlabel('time (seconds)', 'fontsize', 17, 'interpreter', 'latex')
ylabel('$\overline{s_t}$ (mm/s)',  'fontsize', 18, 'interpreter', 'latex')
% % % 
% % 
% % 
subplot(3,1,3)
plot_speed(ma_s_CC,'b-.', 0,2,1,1, 0)
plot_speed(ma_s_TC,'b-', 5071,2,1,1, 0)

plot_speed(ma_s_C,'r-o', offset, linewidth, int, markersize, clip)
plot_speed_error_bands(fr, ma_s_C, ma_se_C, offset, clip, 1)
plot_speed(ma_s_cC,'r-o', 0, linewidth, int, markersize, 0)
plot_speed_error_bands(fr, ma_s_cC, ma_se_cC, 0, 0, 1)
ylim([0, 17])
xlim([0 330])
xlabel('time (seconds)', 'fontsize', 17, 'interpreter', 'latex')
ylabel('$\overline{s_t}$ (mm/s)',  'fontsize', 18, 'interpreter', 'latex')

alpha(0.4)
%% validation plot fitting colony a b c mean s over t
fps=30;%frames per second for figures
mppx=0.145177;%mm/pxl
offset = 5071;
linewidth = 1;
int=200;
markersize = 3;
clip = 350 

%treatment simulations
% [ma_s_A, ma_se_A, pdat_A, fdat_A] = get_CSA_output(0.001, 40, ms2pf(5,mppx, fps), 0.1, 0.1);% sd, Rcrit, seq_data, a, a0
% [ma_s_B, ma_se_B, pdat_B, fdat_B] = get_CSA_output(0.00025, 60, ms2pf(4,mppx, fps), 0.06, 0.1);% sd, Rcrit, seq_data, a, a0
% [ma_s_C, ma_se_C, pdat_C, fdat_C] = get_CSA_output(0.00022, 60, ms2pf(2,mppx, fps), 0.3, 0.1);% sd, Rcrit, seq_data, a, a0

% control simulations
% [ma_s_cA, ma_se_cA] = get_CSA_output(0.001, 40, ms2pf(5,mppx, fps), 0.1, 0);% sd, Rcrit, seq_data, a, a0
% [ma_s_cB, ma_se_cB] = get_CSA_output(0.00025, 60, ms2pf(4,mppx, fps), 0.06, 0);% sd, Rcrit, seq_data, a, a0
% [ma_s_cC, ma_se_cC] = get_CSA_output(0.00022, 40, ms2pf(2,mppx, fps), 0.25, 0);% sd, Rcrit, seq_data, a, a0



% subplot(3,1,1)
% plot_speed(ma_s_CA,'b-.', 0,2,1,1, 0)
% plot_speed(ma_s_TA,'b-', 5071,2,1,1, 0)
% 
% plot_speed(ma_s_A,'r-o', offset, linewidth, int, markersize, clip)
% plot_speed(ma_s_cA,'r-o', 0, linewidth, int, markersize, 0)
% 
% xlabel('time (s)', 'fontsize', 14)
% ylabel('mm/s',  'fontsize', 14); ylim([0, 15.5])
% legend('$\overline{s_t}$ pre-introduction', '$\overline{s_t}$ post-introduction', '$\overline{s_t}$ model', 'Location', 'northwest','FontSize', 14, 'Box', 'off', 'interpreter', 'latex')


subplot(3,1,2)
plot_speed(ma_s_CB,'b-.', 0,2,1,1, 0)
plot_speed(ma_s_TB,'b-', 5071,2,1,1, 0)

plot_speed(ma_s_B,'r-o', offset, linewidth, int, markersize, clip)
% plot_speed(ma_s_cB,'r-o', 0,linewidth, int, markersize, 0)
ylim([0, 15.5])
xlabel('time (s)', 'fontsize', 14)
ylabel('mm/s',  'fontsize', 14)
% 


% subplot(3,1,3)
% plot_speed(ma_s_CC,'b-.', 0,2,1,1, 0)
% plot_speed(ma_s_TC,'b-', 5071,2,1,1, 0)
% 
% plot_speed(ma_s_C,'r-o', offset, linewidth, int, markersize, clip)
% plot_speed(ma_s_cC,'r-o', 0,linewidth, int, markersize, 0)
% ylim([0, 15.5])
% xlabel('time (s)', 'fontsize', 14)
% ylabel('mm/s',  'fontsize', 14)
%% correlations s vs std s
fdat_mod = [fdat_A(:,3:4); fdat_B(:,3:4); fdat_C(:,3:4)];
fdat_emp = [fdat_IT4; fdat_IT45; fdat_IT71];
plot(fdat_mod(:,1), fdat_mod(:,2), 'k+','markersize',1)
hold on
plot(fdat_emp(:,2), fdat_emp(:,3), 'bo','markersize',1)
%% get gammas data 
Tfinal = 4860;
window_edge=2:900:Tfinal+2; %edges in frames for each segment, nonoverlapping
s_edge = linspace(0,29,31);% global bins for s, for all reps and params
fps=30;
VFREQ_emp= strf2Vfreq(STRF_T, window_edge, s_edge);
[fitpar_emp, Vfreq_pdf, bins, Vfreq_se]=VFREQ2gamma_of_t(VFREQ_emp(:,:), s_edge);%collapse all reps and fit each of the windows' histogram to a line, then get gamma for each line
%% overlay empirical data
errorbar(vio_bins/fps, fitpar_emp(:,1), fitpar_emp(:,2), 's-', 'linewidth', 2); hold on
%% plot s_t_bar savg PDAT for empirical and model
fps=30;%frames per second for figures
mppx=0.145177;%mm/pxl

% [~, ~, pdat] = get_CSA_output(0.001, 40, ms2pf(5,mppx, fps), 0.1, 0.1);% sd, Rcrit, seq_data, a, a0
Reps=1;
STBAR_T = zeros(60, Reps);
STBAR_C = zeros(60, Reps);
for r = 1:Reps
    r
    [~, ~, pdat] = get_CSA_output(0.0002, 20, ms2pf(1,mppx, fps), 0.001, 0.01);% sd, Rcrit, seq_data, a, a0
    STBAR_T(:,r) = pdat(:,9);
    [~, ~, pdat] = get_CSA_noAlarm_output(0.0002, 40, ms2pf(5,mppx, fps), 0.001, 0.01);% sd, Rcrit, seq_data, a, a0
    STBAR_C(:,r) = pdat(:,9);
end
%%
figure(1)
histogram(pf2ms(STBAR_T(:),mppx, fps), 0:1:30, 'normalization', 'count')
ylim([0 25])
figure(2)
histogram(pf2ms(STBAR_C(:),mppx, fps), 0:1:30, 'normalization', 'count')
ylim([0 25])
%% run replicates for each plot
fps=30;%frames per second for figures
mppx=0.145177;%mm/pxl

%treatment simulations
[s_A] = get_CSA_output_reps(0.001, 40, ms2pf(5,mppx, fps), 0.1, 0.1);% sd, Rcrit, seq_data, a, a0
disp('done s_A')
[s_B] = get_CSA_output_reps(0.00025, 60, ms2pf(4,mppx, fps), 0.06, 0.1);% sd, Rcrit, seq_data, a, a0
disp('done s_B')
[s_C] = get_CSA_output_reps(0.00022, 60, ms2pf(2,mppx, fps), 0.3, 0.1);% sd, Rcrit, seq_data, a, a0
disp('done s_C')

% control simulations
% [s_cA] = get_CSA_output_reps(0.001, 40, ms2pf(5,mppx, fps), 0.1, 0);% sd, Rcrit, seq_data, a, a0
% [s_cB] = get_CSA_output_reps(0.00025, 60, ms2pf(4,mppx, fps), 0.06, 0);% sd, Rcrit, seq_data, a, a0
[s_cC] = get_CSA_output_reps(0.00022, 40, ms2pf(2,mppx, fps), 0.25, 0);% sd, Rcrit, seq_data, a, a0
disp('done all')

%% plot resplicates result
Tfinal=5100;
fr=(1:Tfinal)';
offset = 5071;
linewidth = 1;
int=1;
markersize = 3;
clip = 350 

subplot(3,1,1)
plot_speed(ma_s_CA,'b-.', 0,2,1,1, 0) %empirical plots dont modify
plot_speed(ma_s_TA,'b-', 5071,2,1,1, 0)
% 
plot_speed(s_A(:,1),'r-', offset, linewidth, int, markersize, clip)
plot_speed_error_bands(fr, s_A(:,1), s_A(:,2), offset, clip, 10)
w=2.5;
plot_speed(w*s_cC(:,1),'r-', 0, linewidth, int, markersize, 0)
plot_speed_error_bands(fr, w*s_cC(:,1), s_cC(:,2), 0, 0, 10)
ylim([0, 17])
xlim([0 330])
xlabel('time (seconds)', 'fontsize', 17, 'interpreter', 'latex')
ylabel('$\overline{s_t}$ (mm/s)',  'fontsize', 18, 'interpreter', 'latex')
legend('$\overline{s_t}$ pre-introduction', '$\overline{s_t}$ post-introduction', '$\overline{s_t}$ model', 'Location', 'northwest','FontSize', 14, 'Box', 'off', 'interpreter', 'latex')


subplot(3,1,2)
plot_speed(ma_s_CB,'b-.', 0,2,1,1, 0)
plot_speed(ma_s_TB,'b-', 5071,2,1,1, 0)

plot_speed(s_B(:,1),'r-', offset, linewidth, int, markersize, clip)
plot_speed_error_bands(fr, s_B(:,1), s_B(:,2), offset, clip, 10)
w=2;
plot_speed(w*s_cC(:,1),'r-', 0, linewidth, int, markersize, 0)
plot_speed_error_bands(fr, w*s_cC(:,1), s_cC(:,2), 0, 0, 10)

ylim([0, 17])
xlim([0 330])
xlabel('time (seconds)', 'fontsize', 17, 'interpreter', 'latex')
ylabel('$\overline{s_t}$ (mm/s)',  'fontsize', 18, 'interpreter', 'latex')
% % % 
% % 
% % 
subplot(3,1,3)
plot_speed(ma_s_CC,'b-.', 0,2,1,1, 0)
plot_speed(ma_s_TC,'b-', 5071,2,1,1, 0)

plot_speed(s_C(:,1),'r-', offset, linewidth, int, markersize, clip)
plot_speed_error_bands(fr, s_C(:,1), s_C(:,2), offset, clip, 10)
w=1;
plot_speed(w*s_cC(:,1),'r-', 0, linewidth, int, markersize, 0)
plot_speed_error_bands(fr, w*s_cC(:,1), s_cC(:,2), 0, 0, 10)
ylim([0, 17])
xlim([0 330])
xlabel('time (seconds)', 'fontsize', 17, 'interpreter', 'latex')
ylabel('$\overline{s_t}$ (mm/s)',  'fontsize', 18, 'interpreter', 'latex')

alpha(0.4)

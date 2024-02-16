%% Get Alarm Replicates for param sweeps and P(s) distribution fitting
fps=30;
Tinject = 2; %= 4500/6 frame to introdue alarmed ants
clear par
par.Tfinal=5100;%fps*150 + Tinject;%#number frames
Reps=1;
par.Np=60;%[20, 60, 180];%linspace(60,180,3);%20:40:200;%
% Np=par.Np;
fps=30;%frames per second for figures
mppx=0.145177;%mm/pxl
%CHANGE these if running replicates corresponding to 'get_CSA_output' (one
%param setting per colony)

% Set main params---------
%LAMBDA-STAR:
par.seq_data = 1.87;%0.5513;% %mean speed over all ants pxl/frame. used as the equilib of the speed params
%R
par.Rcrit=58.2;%(in pxls) 40*mppx=5.8mm
%a
par.a=0.1497;%0.1%;%linspace(0.01, 0.2, 5);
%BETA
par.sd=0.00032;%linspace(2E-4,0.002,3);%0.001
%rho
par.anaught = 0.1;%[0.01, 0.06, 0.11];

% seq_data=ms2pf(5,mppx, fps); 
mppx=0.145177;%mm/pxl

%LAMDA_max
par.kmax=6.8881;%pxl/timestep the value is necessary to max the max 30fps which is the last bin in the s_avg data empirical
% strfr=cell(length(par.Np),1);
% fdatr=cell(length(par.a),1);
% clear fdatr strfr 

% window_edge=2:900:par.Tfinal+2; %edges in frames for each segment, nonoverlapping
% s_edge = linspace(0,29,31);% global bins for s, for all reps and params
% Pos = cell(5,3); %we use 5 time windows and 3 param values for each param sweep

%Run model. strf is the set of tracks, organized by frame number
[strf, fdat]=CSAHab_fast(par.Tfinal, par.Np, par.kmax, par.sd, par.Rcrit, par.a, par.seq_data, par.anaught);

%Output looks like this:

% fdat:[frame, number alarmed, mean speed, std speed, mean speed param]
% fdat(t,:)=[t, sum(lst(:,5)==1), mean(norm_array(lst2(:,6:7))), std(norm_array(lst2(:,6:7))), mean(lst2(:,4))];%[frame, number alarmed, mean speed, std speed, mean speed param]
% strf: [t id x y velocity_param (lambda_t^p), alarm state(0/1), vx, vy]
% strf{t,1}=[t*ones(Np,1), lst2]; %=[t id x y velocity_param, alarm state(0/1), vx, vy]
%% plot group mean speed. no smoothing
plot(fdat(:,1), fdat(:,3), '-')
xlabel('frame'); ylabel('$\overline{s_t}$','interpreter', 'latex','fontsize', 16)

%% to do smoothing then plot
kerL=100;
kernel = ones(kerL,1)/kerL;
fr=(floor(kerL/2):par.Tfinal-floor(kerL/2))/fps; %fr for "valid" conv
fr=fr';
ma_s=pf2ms(conv( fdat(:,3), kernel,'valid'), mppx, fps); %smooth it, convert to mm/s
plot(fr, ma_s, '-')
xlabel('frame'); ylabel('$\overline{s_t}$','interpreter', 'latex','fontsize', 16)
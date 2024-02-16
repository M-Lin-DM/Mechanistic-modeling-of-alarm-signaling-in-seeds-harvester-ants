function [strf, fdat]=CSAHab_fast(Tfinal,Np,kmax, sd, Rcrit, a, lamstar, a0)%Number frames, number ants, max oc kick distribution, decay rate of speed param toward Veq, interaction radius, bug Lambda param
clear lst lst2 cind pop contact colsp 

fps=30;%frames per second for figures
mppx=0.145177;%mm/pxl
R=1040/2;%radius of nest. Nest area = circle

%build configuration for first frame
lst=2*R*(rand(Np,2)-0.5);
d2o= @(xy) sqrt(xy(:,1).^2+xy(:,2).^2);%measures distance to origin
rho=d2o(lst);
lst=lst(rho<R,:);
while length(lst)<Np %add in points til desired number reached
    lst=[lst; 2*R*(rand(1,2)-0.5)];
    lst=lst(d2o(lst)<R,:);
end

%Velocity modulation params---------------
% lamstar=1.8739; %mean speed over all ants pxl/frame. used as the equilib of the speed params
Veq=lamstar;%equilibrium velocity, (pxl/frame)
V1=kmax+lamstar; %velocity of Alarmed seeds = 8.7500 (pxl/frame) 

b=a/(1-exp(-a*kmax));%NOTE: using a=0.1 is only approriate when speed is measured in mm/s NOT pxl frame. thus we must either 1) convert to pxl/fram after the fact or 2) convert s_Avg data to pxl/fram and then fit 'a'
v=(0:0.1:pf2ms(kmax,mppx, fps))';%values of big lambda
w=b*exp(-a*v);%probability weights for those values
kick=ms2pf(sample(v,w,Np), mppx, fps);%convert to pxl/fram since im using param 'a' in terms of mm/s fig

theta_bias=0.369;%omega=0.0473 deg^-1 -->21.1416 deg=0.369rad=mean of exponential distrib. 0.6699<--std of two-sided distrib ~= mean of one sided(flipped)
delta_d_s=@(s, veq, beta) -(s-veq)*beta; %gravitational or natural force of decay of speed towards equilibirum speed. the faster they move, more they decay in a timestep

%everyone begins with the same speed param, Veq and is unalarmed:
lst=[(1:Np)' lst, ones(Np,1)*Veq, zeros(Np,1)];%[ID x y speed_param signal_carring_indicator..] mean of exponential distribution with param lambda is 1/lambda
%set initial actual SPEEDS (this info only contained in vx vy)
lst=[lst, random('exp', Veq, [Np,1]).*randvect(Np,2)];%[ID x y speed_param signal_carring_indicator vx vy]. we sample random orientations. AND speeds from exponential distribution
%set alarmed fraction A0. they will be injected at t=Tinject:
% a0=0.1; 
A0=ceil(Np*a0); %changed from 1% to 10

%TIME PARAMETERS
% Tfinal=fps*100;
% plot_decimation=1000;
% pause_time=.01;%plotting pause (s)
% store_decimation=1;%30fps means store one point each second
Tinject=2;%750;%time of injection of alarmed ants. only used with Tinject=750 in one result figure
% contact_pause=15;%frames an ant pauses after contact

% fdat=zeros(floor(Tfinal/store_decimation)+1, 5);%frame data: [frame, #alarmed, TRUE mean speed, TRUE std of speeds, mean of speed-params]%note: in previous code, taking the mean of the speed params as the mean speed is not accurate as speed params change less noisily than actual speeds........ come on
fdat=zeros(Tfinal, 5);
fdat(1,:)=[1, sum(lst(:,5)==1), mean(norm_array(lst(:,6:7))), std(norm_array(lst(:,6:7))), mean(lst(:,4))];

% c=1;%counter for data storage
strf=cell(length(fdat),1);%the main data structure. Each cell stores all workers locations, state, walking style
strf{1,1}=[1*ones(Np,1), lst];%store frame number on left col

% START SIMULATION   %---finished intitialization, contact data
for t=2:Tfinal
%     t
%     strf{t,1}=strf{t-1,1};
%inject alarmed ants:
    if t==Tinject
        if A0>0
            lst(1:A0,5)=1;%an alarm state of 1 means agent has been alarmed and this is irreversible
            lst(1:A0,4)=V1;%set initial speed params of seed ant(s)
        end
    end
    lst2=lst;

    %MOVEMENT STEP
    for p=1:Np%move all workers in this frame. interaction comes later
     [TH, ~]=cart2pol(lst(p,6), lst(p,7));%measure angle of previous velocity vector
     TH2=TH+(round(rand(1))*2-1)*random('exp', theta_bias);%sample random exponential deviation to the left or right.  (round(rand(1))*2-1) samples from the set {-1,1}
     [vx, vy]=pol2cart(TH2, 1);%vector of unit length with new theta
     pos2=lst2(p,2:3)+random('exp', lst2(p,4), [1,1])*[vx vy];%sample speed vector magnitude from exponential distribution with mean lamba^p=lst2(p,4)
     
         while d2o(pos2)>=R%if and while particle is beyond boundary, redo the randvect (velocity vector)
         TH2=TH+(round(rand(1))*2-1)*random('exp', theta_bias);%sample delta theta again. random deviation from straight. 
         [vx, vy]=pol2cart(TH2, 1);%vector of unit length with new theta (in general direction of last v vector)
         pos2=lst2(p,2:3)-random('exp', lst2(p,4), [1,1])*[vx vy];%fixed speed set at speed param + 
         end
         
      lst2(p,2:3)=pos2;%update the position in current frame. 
      lst2(p,6:7)=pos2-lst(p,2:3);%computing xv yv is the only reason we need two versions of lst
    end
    
    %add random step vector to each speed: + pull speeds toward equilibrium
lst2(:,4)=lst2(:,4) + delta_d_s(lst2(:,4), Veq, sd);%PULL/PUSH speed param toward Veq    
%     disp('completed movement step')
%MEASURE CONTACTS, KICK SPEEDS---------------
dist=tril(squareform(pdist(lst2(:,2:3),'euclidean')));%compute distance matrix and take lower triangular part
cind=find((dist<Rcrit)&(dist~=0));%find pairs with separation < Rcrit
if isempty(cind)==1%store all contacts from this frame
    I=[]; J=[];
else
    [I,J]=ind2sub(size(dist),cind);%list of pairs. I and J are row-col indicies within dist
    IDpairs=[lst(I,1) lst(J,1), lst(I,5), lst(J,5)];%[ID of neib i, ID of neib j, alarm/info state of I, info state of J]
    activ_log_i=sum(IDpairs(:,3:4)==repmat([0 1], size(IDpairs,1), 1),2)==2;%logical the size of IDpairs where agent i is being activated by j
    activ_log_j=sum(IDpairs(:,3:4)==repmat([1 0], size(IDpairs,1), 1),2)==2;%logical the size of IDpairs where agent j is being activated by i
    lst2(IDpairs(activ_log_i,1), 4) = lst2(IDpairs(activ_log_i,1), 4) + kick(IDpairs(activ_log_i,1));%kick speeds of contacting agents by thier respective kick amounts
    lst2(IDpairs(activ_log_j,2), 4) = lst2(IDpairs(activ_log_j,2), 4) + kick(IDpairs(activ_log_j,2));
    
    lst2(IDpairs(activ_log_i,1), 5) = 1; lst2(IDpairs(activ_log_j,2), 5) = 1;%change alarm status to activated of contacted ants
    %record directed links A(i ,j) =1 if i excited j. edgelist stored
end

lst=lst2;
%(simulation could end here without plotting)
%record data
% if mod(t,store_decimation)==0
%     c=c+1;
fdat(t,:)=[t, sum(lst(:,5)==1), mean(norm_array(lst2(:,6:7))), std(norm_array(lst2(:,6:7))), mean(lst2(:,4))];%[frame, number alarmed, mean speed, std speed, mean speed param]
strf{t,1}=[t*ones(Np,1), lst2]; %=[t id x y velocity_param, alarm state(0/1), vx, vy]
% end

end
end

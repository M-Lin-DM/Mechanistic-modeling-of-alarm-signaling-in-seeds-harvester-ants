function [ma_s, ma_se, pdat, fdat] = get_CSA_output(sd, Rcrit, seq_data, a, a0)
fps=30;
% Tfinal=4561;%#number frames
Tfinal=5100;
Np=60;

mppx=0.145177;%mm/pxl
kmax=6.8881;%pxl/timestep the value is necess

[strf, fdat]=CSAHab_fast(Tfinal,Np,kmax, sd, Rcrit, a, seq_data, a0);%Number frames, number ants, max oc kick distribut
% smooth it
kerL=100;
kernel = ones(kerL,1)/kerL;
fr=(floor(kerL/2):Tfinal-floor(kerL/2))/fps; %fr for "valid" conv
fr=fr';
ma_s=conv(fdat(:,3), kernel,'valid'); %smooth it, convert to mm/s
ma_se=conv(fdat(:,4), kernel,'valid'); %smooth it, convert to mm/s

%get pdat
strd = strf2strd(strf);
pdat=zeros(Np,9);
for j=1:Np
    
    pdat(j,:)=[mean(strd{j}(:,:),1), mean(norm_array(strd{j}(:,7:8)), 1) ];
end


%%
[XX, YY]=meshgrid(-4:0.01:4,-4:0.01:4);
[TH, R]=cart2pol(XX,YY);
%%
% v=exp(-R./2);
v=R>(3-0.1)&R<(3+0.1);
imagesc(v); axis equal
%
th=exp(-(TH-pi/4).^2/2/(pi/12)^2);
imagesc(th); axis equal
%
VD=v.*th;
s=surf(XX,YY,VD,VD);
s.EdgeColor='none';
% imagesc(-4,-4,VD)
% imagesc(VD); 

axis equal
hold on
plot3(0,0,.1,'wo','linewidth', 1.5)
%% scatter plot VD habituation
T=10000;
lamb=1.8739;%pxl/frame
theta_bias=0.2;
     TH2=pi/4+(round(rand(T,1))*2-1).*random('exp', theta_bias, [T,1]);%random deviation from  with slight bias to the right
     [vx, vy]=pol2cart(TH2, 1);%vector of unit length with new theta
     pos2=random('exp', lamb, [T,1]).*[vx vy];%
% dat=[random('exp',lamb)];
plot(pos2(:,1), pos2(:,2),'k.','markersize', 1); axis equal; xlabel('pxls', 'fontsize', 14); ylabel('pxls', 'fontsize', 14)
axis([-2 4 -2 4]); 
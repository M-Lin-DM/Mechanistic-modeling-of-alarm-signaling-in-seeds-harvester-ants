%% forced based plotting of growing information flow subnetwork
%create growing edgelist
% edglst=[];
% for t=1:length(pfsub(1,:))%cycle through topological frames. at each frame one new edge is added
%     edg1=[nodeIDs_sub(find(pfsub(:,t)==-1)) nodeIDs_sub(find(pfsub(:,t)==1))];%edge points from edglst(:,1) --> edglst(:,2)
%     edglst=[edglst; edg1 1];
%     net{t,1}=edglst;
% end
%% alternatively, create net{t} from contacts
clear net
t_list=unique(contacts(:,1));%list of unique timestamps
net=cell(length(t_list),1);
edglst=[];
max_nodes=50;
for t=1:length(t_list)
    edglst=[edglst; contacts(contacts(:,1)==t_list(t),2:3)];
    if size(edglst,1)>max_nodes;
        edglst=edglst(end-max_nodes+1:end,:);
    end
    net{t,1}=unique(edglst,'rows');
end
%% Force based graph drawing for expanding primary form subnetwork. (From
% info seed) INPUT is a growing edgelist, not an adj matrix
scale=10;
% ns=2;%number of nodes in first network
%Nframes=400;
dtstart=0.01;
p1=15;%repultion force parameter
k=20;%spring const. This generally needs to be strong
m1=1;%assign mass to runners
Rspring=2;%spring relaxed length. If too small relative to p1, repulsion forces will dominate and pull the spring far past Rspring. 
%R=.5;
Vmax=100;
xf=0.8*scale;
C=0.7;%1.129drag coefficient
f_int=5;%program will spend f_int frames equilibrating the current mass-spring system before adding new nodes. (ie f_int simulation frames per each cell in net)!!!
clear trks_list
trks=[];
for p=1:numel(net{1,1})
trks=[trks; net{1,1}(p), scale*rand(1,3)];%=[ID of first node, x, y, z]
end

trks(1:numel(net{1,1}),12)=m1;%assign mass to nodes
trks(:,5:7)=zeros(numel(net{1,1}),3);%initial velocity
trks(:,8:10)=zeros(numel(net{1,1}),3);%FORCE, NOT acceleration
trks_list{1}=trks;%main data structure
extraframes=0;
Nframes=ceil(f_int*length(net))+extraframes;
% run computation
clear com trks_list
for f=2:Nframes
    f%note t is equvalent to the main index of net
    t0=ceil((f-1)/f_int);%t of previous frame
    t=ceil(f/f_int);%t=frame within the network's topological evolution. MUST ADJUST f_int SO THAT T REACHES LENGTH(NET) BY THE TIME f=Nframes. otherwise youre not seeing the whole network! 
    %If you have T topological frames, then Nframes must = T*f_int OR f_int=Nframes/T
    if t<=length(net)%then add a new node if necessary. used before the 'extrafrmes' stage
        tnodes=unique(net{t,1}(:,1:2));%=list of unique nodes in this topological frame
        tnodes0=unique(net{t0,1}(:,1:2));
        lostnodes=setdiff(union(tnodes,tnodes0), tnodes);%IDs of nodes lost since last frame
        newnodes=setdiff(union(tnodes,tnodes0), tnodes0);%IDs of new nodes in this frame. it may be one or two
        if (t>t0)&&((isempty(lostnodes)==0)||(isempty(newnodes)==0))%if at frame f, the current t > t(at f-1), and the current network has more OR LESS nodes than at f-1
            trks(ismember(trks(:,1),lostnodes),:)=[];%delete lost nodes from trks
            trks=[trks; newnodes zeros(size(newnodes,1),3), zeros(size(newnodes,1),7), m1*ones(size(newnodes,1),1)];%add row for each new node
        end
        %trim edges that happened sufficiently far in the past
        if size(trks,1)>max_nodes;
            trks=trks(end-max_nodes+1:end,:);
        end

    %     [C,IA,~] = unique(trks(:,1));%TESTTTT 
    %     trks=trks(IA,:);

    else t=length(net);%only look at the last topological frame when plotting "extra" frames.
    end

    clear Dists
     tnodes2=trks(:,1);%node [ID] list
     n=length(tnodes2);
     trks(:,8:10)=zeros(n,3);%reset forces to zero for each frame before redefining them!
     dt=dtstart;
    Dists=squareform(pdist(trks(:,2:4),'euclidean'));%distance matrix
    mvect=[0 0 0];
    for p=1:n %for particle p
%             if (adj(p,:)==zeros(1,n))||(isempty(find(leafnodes(:)==p, 1))==1)%for now,do nothing to disjoint nodes. HOW TO NOT ACT ON DISCONNECTED SUBGRAPHS??
%                 continue
%             end

        %find CENTER OF MASS
        mvect=mvect+[trks(p,2),trks(p,3),trks(p,4)].*trks(p,12);
        %INTERPARTICLE FORCES MACHINERY
        for d=1:n %cycles through all other particles for current particle p
         if d~=p
            p2d=((trks(d,2:4)-trks(p,2:4))/norm(trks(d,2:4)-trks(p,2:4)));%a vector (d-p) points from p to d.
            d2p=-p2d;%((trks(p,2:4)-trks(d,2:4))/norm(trks(p,2:4)-trks(d,2:4)));
            if (ismember([tnodes2(p) tnodes2(d)],net{t,1}(:,1:2),'rows')==0)&&(ismember([tnodes2(d) tnodes2(p)],net{t,1}(:,1:2),'rows')==0)%if particles are not linked
                trks(p,8:10)=trks(p,8:10)-p1*trks(p,12)*trks(d,12)*p2d./(Dists(p,d)+1);%(unit vect from p to d)/R ie repultion is proportional to 1/R
                trks(d,8:10)=trks(d,8:10)-p1*trks(p,12)*trks(d,12)*d2p./(Dists(p,d)+1);%(unit vect from d to p)/R ie repultion is proportional to 1/R
            %trks(p,8:10)=trks(p,8:10)-p1*trks(p,12)*trks(d,12)*p2d;% repultion is proportional to constant factor
            %trks(d,8:10)=trks(d,8:10)-p1*trks(p,12)*trks(d,12)*d2p;% repultion is proportional to constant factor
            
            end
            if ((ismember([tnodes2(p) tnodes2(d)],net{t,1}(:,1:2),'rows')==1)||(ismember([tnodes2(d) tnodes2(p)],net{t,1}(:,1:2),'rows')==1))%if particles are linked, create spring force between them
                trks(p,8:10)=trks(p,8:10)+ trks(p,12)*trks(d,12)*p2d*k*(Dists(p,d)-Rspring);%attraction proportional to k*(x - ideal spring length)
                trks(d,8:10)=trks(d,8:10)+ trks(p,12)*trks(d,12)*d2p*k*(Dists(p,d)-Rspring);%attraction proportional to k*(x - ideal spring length)
            end

%                 if Dists{p,2}(d)<Rspring%add momentum change to collide particles
%                     trks(p,9:11)=[0 0 0];
%                     trks(p,6:8)=trks(p,6:8)+(trks(p,6:8)*(trks(p,12)-trks(d,12))+2*trks(d,12)*trks(d,6:8))./(trks(p,12)+trks(d,12));%adds change in momentum vector;
%                 end
         end

        if sum(trks(p,5:7)==[0 0 0],2)~=3;%DAMPING
        trks(p,8:10)=trks(p,8:10)-C*(trks(p,5:7)./norm(trks(p,5:7)))*norm(trks(p,5:7))^2;%add drag force f=-cv^2
        end
    %KEY STEP: change positions
        trks(p,5:7)=trks(p,5:7)+dt*trks(p,8:10)./trks(p,12);%velocity changed by acceleration component
        trks(p,2:4)=trks(p,2:4)+dt*trks(p,5:7);%position changed by velocity component
        end
    end
    
    if sum(sum(isnan(trks)))>0
        disp('forces too high, tune parameters')
        return
    end
com{f,1}=mvect/sum(trks(:,12));
trks_list{f,1}=trks;
end
%% plot
c=1;

%a=xf;
xf=scale;
% for t=1:length(net)
% net{t}=unique(net{t},'rows');
% end
clf
for f=2:1:Nframes-1
    t=ceil(f/f_int)%t=frame within the network's topological evolution
    if t>length(net)%use for plotting extraframes
        t=length(net);
    end
    n=length(net{t}(:,1));%=number of edges in this frame
    
    figure(1)
    h=plot3(com{f}(1),com{f}(2),com{f}(3),'ko','MarkerSize',3,'MarkerFaceColor',[.1 1 .1]);
    hold on
    %title(num2str(f),'FontSize', 20)
    plot3(trks_list{f}(:,2),trks_list{f}(:,3),trks_list{f}(:,4),'o','MarkerSize',5,'MarkerFaceColor',[.9 .5 .1],'MarkerEdgeColor',[.9 .5 .1]);
%     for l=1:length(trks_list{f}(:,1))
%         text(trks_list{f}(l,2),trks_list{f}(l,3),trks_list{f}(l,4),num2str(trks_list{f}(l,1)),'HorizontalAlignment','right','FontSize',10,'FontWeight','bold','Color',[1 .1 0]);
%     end
    
        for p=1:n%index of edge within net's edgelist
            hold on
%             if isempty(net{t,3}{p})==1
%                 continue
%             end
%         h= plot3([trks_list{f}(ismember(trks_list{f}(:,1),net{t}(p,1)),2) trks_list{f}(ismember(trks_list{f}(:,1),net{t}(p,2)),2)],...
%                         [trks_list{f}(ismember(trks_list{f}(:,1),net{t}(p,1)),3) trks_list{f}(ismember(trks_list{f}(:,1),net{t}(p,2)),3)],...
%                         [trks_list{f}(ismember(trks_list{f}(:,1),net{t}(p,1)),4) trks_list{f}(ismember(trks_list{f}(:,1),net{t}(p,2)),4)],'k-');
        arrow3d(trks_list{f}(trks_list{f}(:,1)==net{t}(p,1),2:4), trks_list{f}(trks_list{f}(:,1)==net{t}(p,2),2:4),...
                            20,'line',[0.17,.1],[5, 3],[0 .7 .5],[0 .7 .5]);%arrow3d(start,stop,ang,line type,p,n,headcol,linecol)
        end
    hold off
 %plot growing track of a single node
    %plot3([trks_list{f-1}([4],2) trks_list{f}([4],2)],[trks_list{f-1}([4],3) trks_list{f}([4],3)],[trks_list{f-1}([4],4) trks_list{f}([4],4)],'k-','MarkerSize',3,'MarkerFaceColor',[.1 1 .1]);%chaser
%     axis([0 xf 0 xf 0 xf])
%     camzoom(1+1.1*f/Nframes)
    camtarget([com{f}(1) com{f}(2) com{f}(3)])
%     axis([-a-c*(f/Nframes) a+c*(f/Nframes) -a-c*(f/Nframes) a+c*(f/Nframes) -a-c*(f/Nframes) a+c*(f/Nframes)])%growing cube
 title(['Time = ' num2str(t_list(ceil(f/f_int)))])
    axis equal
     camorbit(-360*f/Nframes,0*f/Nframes)
    box on
    
    %campos([xf xf xf])
%      pause(0.01)
%    saveas(h, ['C:\Users\m\Documents\KANG LAB\Interacting Particle System\Movies\frame_' num2str(f) '.tif'])
 saveas(h, ['/Users/M.lin/Documents/KANG LAB/Movies/frame_' num2str(f) '.tif'])%mac
  
end
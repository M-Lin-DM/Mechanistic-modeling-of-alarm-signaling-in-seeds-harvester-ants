%% initial cond
%clear all
% adj=edgeL2adj([EL, ones(length(EL),1)], (1:Np)');
Nframes=200;
trks=zeros(Np,13);
trks(:,1)=(1:Np)';


dtstart=0.03;
p1=10;%interaction force parameter
k=10;%spring const
m1=1;%assign mass to runners
Rspring=1;%spring relaxed length
Vmax=100;
xf=0.8*scale;
C=0.29;%drag
clear trks_list
trks(1:n,12)=m1;%assign mass to nodes
trks(:,5:7)=zeros(n,3);%initial velocity
trks(:,8:10)=zeros(n,3);%FORCE, NOT acceleration
trks_list{1}=trks;
%% RUN Computation
for f=2:Nframes
    clear Dists
    trks(:,8:10)=zeros(n,3);%reset forces to zero for each frame before redefining them!
     dt=dtstart;
    for p=1:n%creates list of distance from p to all other particles
        %Dists{p,1}=trks(:,1);%struct with node [ID]%not needed b/c the
        %index p in trks IS the ID
        Dists{p,1}=dist_pt_to_vector(trks(p,2:4),trks(:,2:4));%2nd row of cell array is dist of all pts to point p
        trks(p,13)=norm(trks(p,6:8));%p's speed in this frame before changes
    end
    
   mvect=[0 0 0];
    for p=1:n %for particle p
%             if (adj(p,:)==zeros(1,n))||(isempty(find(leafnodes(:)==p, 1))==1)%for now,do nothing to disjoint nodes. HOW TO NOT ACT ON DISCONNECTED SUBGRAPHS??
%                 continue
%             end
        %find CENTER OF MASS
        mvect=mvect+[trks(p,2),trks(p,3),trks(p,4)].*trks(p,12);
        %INTERPARTICLE FORCES MACHINERY
        for d=1:n %cycles through all other particles for current particle p
            
            if (d~=p)&&(adj(p,d)==0)%if particles are not linked
                trks(p,8:10)=trks(p,8:10)-p1*trks(p,12)*trks(d,12)*((trks(d,2:4)-trks(p,2:4))/norm(trks(d,2:4)-trks(p,2:4)))./Dists{p}(d);
            end
            if (d~=p)&&(adj(p,d)==1)%if particles are linked, create spring force between them
                trks(p,8:10)=trks(p,8:10)+((trks(d,2:4)-trks(p,2:4))/norm(trks(d,2:4)-trks(p,2:4)))*k*(Dists{p}(d)-Rspring);
            end

%                 if Dists{p,2}(d)<Rspring%add momentum change to collide particles
%                     trks(p,9:11)=[0 0 0];
%                     trks(p,6:8)=trks(p,6:8)+(trks(p,6:8)*(trks(p,12)-trks(d,12))+2*trks(d,12)*trks(d,6:8))./(trks(p,12)+trks(d,12));%adds change in momentum vector;
%                 end
        end

        if trks(p,5:7)~=[0 0 0]%DAMPING
        trks(p,8:10)=trks(p,8:10)-C*(trks(p,5:7)./norm(trks(p,5:7)))*norm(trks(p,5:7))^2;%add drag force f=-cv^2
        end
    
         trks(p,5:7)=trks(p,5:7)+dt*trks(p,8:10)./trks(p,12);%velocity changed by acceleration component
         trks(p,2:4)=trks(p,2:4)+dt*trks(p,5:7);%position changed by velocity component
   
    end
com{f}=mvect/sum(trks(:,12));
trks_list{f}=trks;
end
%% plot
c=1;
a=xf;
for f=2:1:Nframes
    figure(1)
    plot3(com{f}(1),com{f}(2),com{f}(3),'ko','MarkerSize',3,'MarkerFaceColor',[1 .1 .1])
    hold on
    plot3(trks_list{f}(:,2),trks_list{f}(:,3),trks_list{f}(:,4),'ko','MarkerSize',7,'MarkerFaceColor',[.1 1 .1]);%chaser
    
        for p=1:n
            if isempty(adjList{p})==1
                continue
            end
            for a=1:length(adjList{p})%
                if adjList{p}(a)>=p %only plots each connection ONCE
                plot3([trks_list{f}(p,2) trks_list{f}(adjList{p}(a),2)],[trks_list{f}(p,3) trks_list{f}(adjList{p}(a),3)],[trks_list{f}(p,4) trks_list{f}(adjList{p}(a),4)],'k-')
                end
            end
        end
    hold off
 %plot growing track of a single node
    %plot3([trks_list{f-1}([4],2) trks_list{f}([4],2)],[trks_list{f-1}([4],3) trks_list{f}([4],3)],[trks_list{f-1}([4],4) trks_list{f}([4],4)],'k-','MarkerSize',3,'MarkerFaceColor',[.1 1 .1]);%chaser
    camzoom(1+1.5*f/Nframes)
    camtarget([com{f}(1) com{f}(2) com{f}(3)])
    AXIS([0 xf 0 xf 0 xf])
    axis([-a-c*(f/Nframes) a+c*(f/Nframes) -a-c*(f/Nframes) a+c*(f/Nframes) -a-c*(f/Nframes) a+c*(f/Nframes)])%growing cube
    axis equal
    camorbit(90*f/Nframes,0*f/Nframes)
    box on
    %campos([xf xf xf])
    
   %saveas(h, ['C:\Users\The Owner Michael\Documents\Experiments\Networks\movies\frame_' num2str(f) '.tif'])
  
end
%% plot com
 figure(1)
%camorbit(0,45)
for f=2:1:Nframes-1
    figure(1)
    plot3(com{f}(1),com{f}(2),com{f}(3),'ko','MarkerSize',3,'MarkerFaceColor',[1 .1 .1])
%     plot3(trks_list{f}([4,6,8],2),trks_list{f}([4,6,8],3),trks_list{f}([4,6,8],4),'ko','MarkerSize',3,'MarkerFaceColor',[.1 1 .1]);%chaser
    %hold on%CAMORBIT AND A BUNCH OF OTHER PLOT TOOLS WONT WORK WITH THIS ON. ONLY USE IF PLOTTING A GROWING TRACK
    
    %plot little two-frame line(looks like a flying bug!)
    hold on
    plot3([trks_list{f-1}([1],2) trks_list{f+1}([1],2)],[trks_list{f-1}([1],3) trks_list{f+1}([1],3)],[trks_list{f-1}([1],4) trks_list{f+1}([1],4)],'-','LineWidth',3,'MarkerEdgeColor',[1 0 .1]);%
    plot3([trks_list{f-1}([2],2) trks_list{f+1}([2],2)],[trks_list{f-1}([2],3) trks_list{f+1}([2],3)],[trks_list{f-1}([2],4) trks_list{f+1}([2],4)],'-','LineWidth',3,'MarkerEdgeColor',[0 1 .1]);%
    plot3([trks_list{f-1}([3],2) trks_list{f+1}([3],2)],[trks_list{f-1}([3],3) trks_list{f+1}([3],3)],[trks_list{f-1}([3],4) trks_list{f+1}([3],4)],'-','LineWidth',3,'MarkerEdgeColor',[.1 0 1]);%
    hold off
    %plot growing track of a single node
    %plot3([trks_list{f-1}([4],2) trks_list{f}([4],2)],[trks_list{f-1}([4],3) trks_list{f}([4],3)],[trks_list{f-1}([4],4) trks_list{f}([4],4)],'k-','MarkerSize',3,'MarkerFaceColor',[.1 1 .1]);%chaser
    
    camzoom(1+1.5*f/Nframes)
    camtarget([com{f}(1) com{f}(2) com{f}(3)])
    %AXIS([0 xf 0 xf 0 xf])
    axis([-a-c*(f/Nframes) a+c*(f/Nframes) -a-c*(f/Nframes) a+c*(f/Nframes) -a-c*(f/Nframes) a+c*(f/Nframes)])%growing cube
    %axis equal
    camorbit(90*f/Nframes,0*f/Nframes)
    box on
    %campos([xf xf xf])
    pause(0.03)
end
%%
%%
%%
%%
%%
%%


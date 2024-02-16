%% initial cond
%clear all
% Np=60;
n=Np;
adj=edgeL2adj([EL, ones(length(EL),1)], (1:n)');
adjList=adj2adjL(adj);
Nframes=800;
trks=zeros(n,10);
trks(:,1)=(1:n)';
trks(:,2:3)=10*rand(n,2)-5;
trks(:,4:5)=zeros(n,2);%initial velocity
trks(:,6:7)=zeros(n,2);%initial force

dtstart=0.06;
p1=@(t) 4*(t/Nframes)^1;%euclidean distance repultion interaction force parameter
% Rrepl=3;%radius under which node repultion occurs
Rrepl =@(t) 5;%5*(t/Nframes)^1;%gradualling increase the repulsion radius

% p2=@(t) 2*(1-(t/Nframes));%graph distance repulsion parameter
p3=@(t) 0.2*(1-(t/Nframes));% repulsion based on absolute difference between betweeness centralities
DRrepl =@(t) ceil(40*(1-(t/Nframes))+4); %topological distance above which repulsion occurs
k=5;%spring const
m1=1;%assign mass to runners
Rspring=1;%spring relaxed length
C=0.1;%drag

clear trks_list
trks(1:n,10)=m1;%assign mass to nodes
trks_list{1}=trks;
%%
%get path lengths
D=zeros(n,n);
for p=1:n
D(p,:) = simple_dijkstra(symmetrize(adj),p);%path lengths in undirected network. you can get lengths through directed graph too!!
end
histogram(D)%check distrib of path lengths then decide cutoff for repulsion
% %get leaves
% betw = node_betweenness_slow(symmetrize(adj));
betw = node_betweenness_faster(symmetrize(adj));

%% RUN Computation
for f=2:Nframes
    f
    clear Dists
    trks(:,6:7)=zeros(n,2);%reset forces to zero for each frame before redefining them!
    dt=dtstart;
    Dists=squareform(pdist(trks(:,2:3)));
    
   mvect=[0 0];
    for p=1:n %for particle p
%             if (adj(p,:)==zeros(1,n))||(isempty(find(leafnodes(:)==p, 1))==1)%for now,do nothing to disjoint nodes. HOW TO NOT ACT ON DISCONNECTED SUBGRAPHS??
%                 continue
%             end
        trks(p,8)=norm(trks(p,6:8));%p's speed in this frame before changes

        %find CENTER OF MASS
        mvect=mvect+[trks(p,2),trks(p,3)].*trks(p,10);
        %INTERPARTICLE FORCES MACHINERY
        for d=1:n %cycles through all other particles for current particle p
            
            if (d~=p)&&(adj(p,d)==0)&&(Dists(p,d)<Rrepl(f))%if particles are not close, repel ('close' based first on euclidean distance for small scale structures, then based on graph distance for long scale network structure)
                trks(p,6:7)=trks(p,6:7)-p1(f)*((trks(d,2:3)-trks(p,2:3))/norm(trks(d,2:3)-trks(p,2:3)))./Dists(p,d)^0.5 ;%- p2(f)*((trks(d,2:3)-trks(p,2:3))/norm(trks(d,2:3)-trks(p,2:3)))*D(p,d)^0.5
            end
            if (d~=p)&&((adj(p,d)==1)||(adj(d,p)==1))%if particles are directly linked, create spring force between them
                trks(p,6:7)=trks(p,6:7)+((trks(d,2:3)-trks(p,2:3))/norm(trks(d,2:3)-trks(p,2:3)))*k*(Dists(p,d)-Rspring);
            end
            if (d~=p)&&(adj(p,d)==0)&&(D(p,d)>DRrepl(f))&&(D(p,d)<Inf)%if particles are not close topologically repel based on difference between betweeness centralities
                trks(p,6:7)=trks(p,6:7)- p3(f)*((trks(d,2:3)-trks(p,2:3))/norm(trks(d,2:3)-trks(p,2:3)))*abs(betw(d)-betw(p))^1;
            end

%                 if Dists{p,2}(d)<Rspring%add momentum change to collide particles
%                     trks(p,9:11)=[0 0 0];
%                     trks(p,6:8)=trks(p,6:8)+(trks(p,6:8)*(trks(p,10)-trks(d,10))+2*trks(d,10)*trks(d,6:8))./(trks(p,10)+trks(d,10));%adds change in momentum vector;
%                 end
        end

        if trks(p,4:5)~=[0 0]%DAMPING
        trks(p,6:7)=trks(p,6:7)-C*trks(p,4:5)./norm(trks(p,4:5))*norm(trks(p,4:5))^2;%add drag force f=-cv^2
        end
    
         trks(p,4:5)=trks(p,4:5)+dt*trks(p,6:7)./trks(p,10);%velocity changed by acceleration component
         trks(p,2:3)=trks(p,2:3)+dt*trks(p,4:5);%position changed by velocity component
   
    end
com{f}=mvect/sum(trks(:,10));
trks_list{f,1}=trks;
end
%% plot
% drawArrow = @(x,y) quiver( x(1),y(1),x(2)-x(1),y(2)-y(1),0,'color', [0 0 0],'AutoScaleFactor', 2);
seed=EL(1,1);
whitebg([1 1 1])
cm=jet(200);
c=1;
for f=200:5:Nframes%2:1:Nframes
    figure(1)
    plot(com{f}(1),com{f}(2),'.','MarkerSize',1,'MarkerFaceColor',[1 .1 .1])
    hold on
    quiver(trks_list{f}(EL(:,1), 2), trks_list{f}(EL(:,1), 3), trks_list{f}(EL(:,2), 2)-trks_list{f}(EL(:,1), 2), trks_list{f}(EL(:,2), 3)-trks_list{f}(EL(:,1), 3),0,'color', [0 0 0],'markersize', 4)
%         for p=1:n
%             if isempty(adjList{p})==1
%                 continue
%             end
%             for a=1:length(adjList{p})%
% %                 if adjList{p}(a)>=p %only plots each connection ONCE
% %                 plot([trks_list{f}(p,2) trks_list{f}(adjList{p}(a),2)],[trks_list{f}(p,3) trks_list{f}(adjList{p}(a),3)],'k-')
%                 drawArrow([trks_list{f}(p,2) trks_list{f}(adjList{p}(a),2)], [trks_list{f}(p,3) trks_list{f}(adjList{p}(a),3)])
% %                 end
%             end
%         end

%     scatter(trks_list{f}(:,2),trks_list{f}(:,3), 20,  map2colsp(betw,cm,[0,1.5]),'filled')
            plot(trks_list{f}(:,2),trks_list{f}(:,3),'ko','MarkerSize',5,'MarkerFaceColor',[.1 1 .1]);%chaser
%     plot(trks_list{f}(seed,2),trks_list{f}(seed,3),'ko','MarkerSize',5,'MarkerFaceColor',[1 0 .2]);

    hold off
 %plot growing track of a single node
    %plot3([trks_list{f-1}([4],2) trks_list{f}([4],2)],[trks_list{f-1}([4],3) trks_list{f}([4],3)],[trks_list{f-1}([4],4) trks_list{f}([4],4)],'k-','MarkerSize',3,'MarkerFaceColor',[.1 1 .1]);%chaser
%     camzoom(1+1.5*f/Nframes)
%     AXIS([0 xf 0 xf])
%     axis([-a-c*(f/Nframes) a+c*(f/Nframes) -a-c*(f/Nframes) a+c*(f/Nframes) -a-c*(f/Nframes) a+c*(f/Nframes)])%growing cube
    axis equal
    %campos([xf xf xf])
    title(num2str(f))
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


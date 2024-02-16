function Vfreq_over_windows = strf2Vfreq(strf, window_edge, s_edge)
%input: window_edge, s_edge must be in units of pxl/fram
Vfreq_over_windows = zeros(length(window_edge)-1, length(s_edge) -1);
 strd=strf2strd(strf);
    Vtracks=[];
    for p=1:length(strd)
        Vtracks=[Vtracks; strd{p}(:,[1,7:8])]; %[t, vx, vy] grab only speed data for memory saving
    end
    Vtracks(:,4)=norm_array(Vtracks(:,2:3));%compute step size or speed
%     if any(Vtracks(:,4)>=s_edge(end))
%         disp('increase max speed bin edge')
%         Vtracks(Vtracks(:,4)>=s_edge(end), 4)
%     end
    
    for t=1:length(window_edge)-1
    window=(Vtracks(:,1)>=window_edge(t))&(Vtracks(:,1)<window_edge(t+1));
    [Vfreq, ~]=histcounts(Vtracks(window,4), s_edge, 'normalization','count'); %we will normalize to a pdf later after all reps. it can cause problems to estimate pdf now
    Vfreq_over_windows(t,:) = Vfreq;
    end
end
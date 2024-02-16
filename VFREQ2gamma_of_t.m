function [out, VFREQ_pdf, bins, VFREQ_se] = VFREQ2gamma_of_t(VFREQ, s_edge)
%input: VFREQ =[#windows x #speed bins x #reps]
%s_edge (in units of pxl/frame)
%input a stack of replicates-->sum and convert to single pdf (over time)-->fit each
%time window-->output fit parameters over time window
fps=30;%frames per second for figures
mppx=0.145177;%mm/pxl
bins=(s_edge(1:end-1) + s_edge(2:end))/2; %get bin centers
bin_width = bins(2)-bins(1);
bins=pf2ms(bins, mppx, fps); bin_width=pf2ms(bin_width, mppx, fps); %convert bin width in order to convert to pdf over the mm/s space
reps= size(VFREQ,3);
%gen mean frequencies and se of mean
VFREQ_pdf = sum(VFREQ,3)/reps;%sum over all replicates in the stack
VFREQ_se = std(VFREQ,0,3)/sqrt(reps); 
%now convert from mean frequency to pdf
window_sum = sum(VFREQ_pdf,2);
VFREQ_pdf= VFREQ_pdf./(window_sum*bin_width); %get pdf estimate, OVER the mm/s space
VFREQ_se = VFREQ_se./(window_sum*bin_width);
%
out= zeros(size(VFREQ,1), 4);
for w = 1:size(VFREQ,1)
    logp = [bins; log(VFREQ_pdf(w,:))];
    logp = logp(:, ~isnan(logp(2,:))); logp = logp(:, ~isinf(logp(2,:))); %ignore nan and inf
    [m, b, sm, sb] = linearfit(logp(1,:)', logp(2,:)');
    out(w,:) = [m, b, sm, sb];
end

end
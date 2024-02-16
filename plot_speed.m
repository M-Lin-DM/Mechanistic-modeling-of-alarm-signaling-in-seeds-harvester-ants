function [] = plot_speed(ma_s, linestyle, offset, linewidth, int, markersize, clip)
fps=30;%frames per second for figures
mppx=0.145177;%mm/pxl

plot(((1:int:length(ma_s)-clip)+offset)/fps, pf2ms(ma_s(1:int:length(ma_s)-clip), mppx, fps), linestyle,'linewidth', linewidth, 'markersize', markersize); hold on
end
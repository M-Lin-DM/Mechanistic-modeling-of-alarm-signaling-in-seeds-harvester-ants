function [] = plot_speed_error_bands_quantiles(frames, pt25, pt75, offset, clip, RGB)
fps=30;%frames per second for figures
mppx=0.145177;%mm/pxl
% cm=copper(3);

frames = frames(1:end-clip);
frames = (frames+offset)/fps;


fill([frames; flipud(frames)], [pt25(1:end-clip); flipud(pt75(1:end-clip))], RGB, 'linestyle','none');  hold on %"BANDZ" error bars for ultimate style
alpha(0.4)
end
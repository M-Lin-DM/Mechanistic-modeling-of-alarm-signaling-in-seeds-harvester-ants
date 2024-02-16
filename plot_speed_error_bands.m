function [] = plot_speed_error_bands(frames, ma_s, ma_se, offset, clip, scale, RGB)
fps=30;%frames per second for figures
mppx=0.145177;%mm/pxl
% cm=copper(3);

frames = frames(1:end-clip);
frames = (frames+offset)/fps;

ma_s = pf2ms(ma_s, mppx, fps);
ma_se = pf2ms(ma_se, mppx, fps);

fill([frames; flipud(frames)], [ma_s(1:end-clip)-scale*ma_se(1:end-clip); flipud(ma_s(1:end-clip)+scale*ma_se(1:end-clip))], RGB, 'linestyle','none');  hold on %"BANDZ" error bars for ultimate style
alpha(0.4)
end
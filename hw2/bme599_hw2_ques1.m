% BME 599, HW2
% before running, must add epg_funcs_hargreaves to path
addpath('C:\Users\amaya\OneDrive\Documents\MATLAB\bme599\amurguia\hw2\epg_funcs_hargreaves\');

% code in epg_funcs_hargreaves file is from here: https://web.stanford.edu/~bah/software/epg/

% 1ai

seq1 = epg_cpmg(pi, 64, 200, 50, 5);
seq2 = epg_cpmg(pi, 64, 600, 150, 5);
seq3 = epg_cpmg(pi, 64, 900, 200, 5);
seq4 = epg_cpmg(pi, 64, 1200, 250, 5);
seq5 = epg_cpmg(pi, 64, 1500, 300, 5);

TE_ms_vec = 5:5:5*64;

figure
plot(TE_ms_vec,abs(seq1), TE_ms_vec, abs(seq2), TE_ms_vec,abs(seq3), TE_ms_vec,abs(seq4), TE_ms_vec,abs(seq5));
legend("T1 200ms, T2 50ms", "T1 600ms, T2 150ms", "T1 900ms, T2 200ms", "T1 1200ms, T2 250ms", "T1 1500ms, T2 300ms");

%% 1aii

seq1 = epg_cpmg(120*pi/180, 64, 200, 50, 5);
seq2 = epg_cpmg(120*pi/180, 64, 600, 150, 5);
seq3 = epg_cpmg(120*pi/180, 64, 900, 200, 5);
seq4 = epg_cpmg(120*pi/180, 64, 1200, 250, 5);
seq5 = epg_cpmg(120*pi/180, 64, 1500, 300, 5);

figure
plot(TE_ms_vec,abs(seq1), TE_ms_vec, abs(seq2), TE_ms_vec,abs(seq3), TE_ms_vec,abs(seq4), TE_ms_vec,abs(seq5));
legend("T1 200ms, T2 50ms", "T1 600ms, T2 150ms", "T1 900ms, T2 200ms", "T1 1200ms, T2 250ms", "T1 1500ms, T2 300ms");

%% 1aiii

seq1 = epg_cpmg(60*pi/180, 64, 200, 50, 5);
seq2 = epg_cpmg(60*pi/180, 64, 600, 150, 5);
seq3 = epg_cpmg(60*pi/180, 64, 900, 200, 5);
seq4 = epg_cpmg(60*pi/180, 64, 1200, 250, 5);
seq5 = epg_cpmg(60*pi/180, 64, 1500, 300, 5);

figure
plot(TE_ms_vec,abs(seq1), TE_ms_vec, abs(seq2), TE_ms_vec,abs(seq3), TE_ms_vec,abs(seq4), TE_ms_vec,abs(seq5));
legend("T1 200ms, T2 50ms", "T1 600ms, T2 150ms", "T1 900ms, T2 200ms", "T1 1200ms, T2 250ms", "T1 1500ms, T2 300ms");

%% 1b

T1_vec_ms = 200:100:1500;
T2_vec_ms = 50:30:300;
flip_angle_rad = [pi, 120*pi/180, 60*pi/180];
TE_idx_plot_vec = [6,16,32,48];
sig = zeros(length(flip_angle_rad),length(T1_vec_ms),length(T2_vec_ms),length(TE_idx_plot_vec));


for i = 1:length(flip_angle_rad)
    for j = 1:length(T1_vec_ms)
        for k = 1:length(T2_vec_ms)
            seq = epg_cpmg(flip_angle_rad(i), 64, T1_vec_ms(j), T1_vec_ms(k), 5);
            for m = 1:length(TE_idx_plot_vec)
                sig(i,j,k,m) = seq(TE_idx_plot_vec(m));
            end
        end
    end
end

cnt=1;
for i = 1:length(flip_angle_rad)
    for m = 1:length(TE_idx_plot_vec)
        subplot(3,4,cnt)
        contour(T2_vec_ms,T1_vec_ms,squeeze(abs(sig(i,:,:,m))));
        flip = num2str(flip_angle_rad(i));
        echo = num2str(TE_idx_plot_vec(m));
        title(sprintf('Flip = %s, echo %s',flip,echo));
        xlabel("T2(ms)")
        ylabel("T1(ms)")
        cnt=cnt+1;
    end
end

% cnt=1;
% for i = 1:length(flip_angle_rad)
%     for m = 1:length(TE_idx_plot_vec)
%         subplot(3,4,cnt)
%         %contour(T2_vec_ms,T1_vec_ms,squeeze(abs(sig(i,:,:,m))));
%         imagesc(squeeze(abs(sig(i,:,:,m))));
%         flip = num2str(flip_angle_rad(i));
%         echo = num2str(TE_idx_plot_vec(m));
%         title(sprintf('Flip = %s, echo %s',flip,echo));
%         xlabel("T2(ms)")
%         ylabel("T1(ms)")
%         cnt=cnt+1;
%     end
% end

disp("done")
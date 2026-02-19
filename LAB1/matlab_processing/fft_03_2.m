%% import
clear all
close all
freq_spec_1 = readtable("FNM\output_fft\fft_decode_freq_spectrum_1.txt");
freq_spec_1 = table2array(freq_spec_1);

freq_spec_2 = readtable("FNM\output_fft\fft_decode_freq_spectrum_2.txt")
freq_spec_2 = table2array(freq_spec_2);

filt_time_series = readtable("FNM\output_fft\fft_filtered_time_series.txt")
filt_time_series = table2array(filt_time_series);

%% Parameters
N = length(freq_spec_1);
delta = 1/N;

n_arr = -N/2:1:N/2-1;
f_arr = n_arr/(N*delta);

T_arr = linspace(0,1,8192);

time_ser_amp = sqrt(filt_time_series.^2);

%% Plot
close all
f_c = 1024;
f_m = 128;
tiledlayout(1,2,"TileSpacing","compact")
nexttile(1)
plot(f_arr, freq_spec_1, LineWidth=2)
xlim([0 2000])
xline(f_c - f_m, '--r', 'f_c - f_m (f = 896)', 'LabelHorizontalAlignment','left', FontSize=12, Color='black');
xline(f_c + f_m, '--r', 'f_c + f_m (f = 1152)', 'LabelHorizontalAlignment','right', FontSize=12, Color='black');
set(gca, "FontSize", 16)
xlabel('Frequency', FontSize=14)
box on
grid on
title('a)', FontSize=24)

nexttile(2)
plot(f_arr, freq_spec_2, LineWidth=2)
xlim([0 2000])
set(gca, "FontSize", 16)
xlabel('Frequency', FontSize=14)
box on
grid on
title('b)', FontSize=24)
set(gcf,'PaperPositionMode','auto');        % preserves figure size
print(gcf,'C:\Dokument\4an\FNM\plots_fft\fft_decode','-depsc2')

% figure(3)
% plot(T_arr, time_ser_amp)
% 
% max(filt_time_series)
% 
% filt_over = filt_time_series > 2.0;


%%
figure(4)
bit_one = time_ser_amp(1:1024);

plot(T_arr(1:1024),bit_one)

figure(5)

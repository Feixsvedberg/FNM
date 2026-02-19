close all
clear all

%% Import
fft_real = readtable("FNM\output_fft\fft_gauss_real.txt");
fft_real = table2array(fft_real);
fft_complex = readtable("FNM\output_fft\fft_gauss_complex.txt");
fft_complex = table2array(fft_complex);

fft_exact = readtable("FNM\output_fft\fft_exact_trans.txt");
fft_exact = table2array(fft_exact);

%% Plotting preprocessing
close all
N = length(fft_real);
delta = (1/3-(-1/3))/N;
n_arr = -N/2+1:1:N/2;
f_arr = n_arr/(N*delta);

diff = abs(fft_real-fft_exact);

tiledlayout(3,1,"TileSpacing","loose")


nexttile(1)
scatter(f_arr,fft_real,Marker="*")
hold on
plot(f_arr, fft_exact, LineWidth=1)
xlim([-300 300])
box on
grid on
set(gca, FontSize = 12)
title('a)', FontSize=24)
legend('FFT Real Numerical', 'FFT Analytical')
ylim([0 1.1])

nexttile(2)
plot(f_arr, fft_complex, LineWidth=2)
box on
grid on
set(gca, FontSize = 12)
title('b)', FontSize=24)
xlim([-1500 1500])

nexttile(3)
plot(f_arr, diff, LineWidth=2)
box on
grid on
set(gca, FontSize = 12)
title('c)', FontSize=24)
xlabel("Frequency")
xlim([-1500 1500])

set(gcf,'PaperPositionMode','auto');        % preserves figure size
print(gcf,'C:\Dokument\4an\FNM\plots_fft\fft_01','-depsc2')

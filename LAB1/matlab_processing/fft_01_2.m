close all
fft_trans_data = readtable("FNM\output_fft\fft_gauss_2.txt");
fft_trans_data = table2array(fft_trans_data);
N = length(fft_trans_data)/2;

delta = (1/3-(-1/3))/N;
sigma = 1/64;
%Matlab 1 to N in index grrr
N_array = 1:1:N*2;

%Hämta rätt element (varannat).
real_index = mod(N_array(:),2) ~= 0;
complex_index = real_index == 0;

real_part = delta*fft_trans_data(real_index);
%correct time shift

complex_part = delta*fft_trans_data(complex_index);


fft_exact = readtable("FNM\output_fft\fft_exact_trans.txt");
fft_exact = table2array(fft_exact);
%% Krux
n_arr = -N/2:1:N/2-1;

f_arr = n_arr/(N*delta);


figure(1)
plot(f_arr,real_part);
hold on
scatter(f_arr, fft_exact);
legend('numerical','analytical')

figure(2)
plot(f_arr, complex_part);
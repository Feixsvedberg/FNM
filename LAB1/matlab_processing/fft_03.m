close all
fft_trans_data = readtable("FNM\output_fft\fft_decode_1.txt");
fft_trans_data = table2array(fft_trans_data);
N = length(fft_trans_data)/2;

delta = 1/N;

N_array = 1:1:N*2;
%Hämta rätt element (varannat).
real_index = mod(N_array(:),2) ~= 0;
complex_index = real_index == 0;

real_part = delta*fft_trans_data(real_index);
%correct time shift

complex_part = delta*fft_trans_data(complex_index);


fft_filtered = readtable("FNM\output_fft\fft_filtered_signal.txt");
fft_filtered = table2array(fft_filtered);
N = length(fft_filtered)/2;

real_filtered = delta*fft_filtered(real_index);
complex_filtered = delta*fft_filtered(complex_index);
%% Krux
n_arr = -N/2:1:N/2-1;

f_arr = n_arr/(N*delta);

transf_mag = sqrt(real_part.^2+complex_part.^2);
filter_mag = sqrt(real_filtered.^2+complex_filtered.^2);

f_spectrum = transf_mag.^2;
filtered_spectrum = filter_mag.^2;


figure(1)
plot(f_arr,f_spectrum)


figure(2)
plot(f_arr, filtered_spectrum)
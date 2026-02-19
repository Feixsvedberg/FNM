close all
fft_trans_data = readtable("FNM\output_fft\fft_trans_data.txt");
fft_trans_data = table2array(fft_trans_data);
N = length(fft_trans_data)/2;

delta = (1/3-(-1/3))/N;
sigma = 1/64;
%Matlab 1 to N in index grrr
N_array = 1:1:N*2;

%Hämta rätt element (varannat).
real_index = mod(N_array(:),2) ~= 0;
complex_index = real_index == 0;

real_part = fft_trans_data(real_index);
%correct time shift

complex_part = fft_trans_data(complex_index);

%dela upp listan av frekvenser för att matcha gsl output layout
k = 0:N-1;
f = (k .* 1 /(delta*N));
f(k >= N/2) = (k(k >= N/2) - N) .* 1 /(delta*N);
real_part_corrected = real_part .* transpose(exp((1i*2*pi*f*1/3))); 

trans_anal = @(f) exp(-0.25*(2*pi*sigma*f).^2);


f_anal = linspace(-1/(2*delta), 1/(2*delta), 10000);



figure(1)
%lägg märke till omskalning med delta, se sida 1 i labbinstruktioner.
scatter(f,real_part_corrected*delta,Color = 'blue')
hold on
plot(f_anal,trans_anal(f_anal))

figure(2)
plot(f,complex_part)
hold on
plot(f_anal,trans_anal(f_anal))

diff = abs(real_part_corrected*delta-(trans_anal(f))');

figure(3)
plot(f, diff)
xlim([-100 100])














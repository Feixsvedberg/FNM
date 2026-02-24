frequencies_100 = table2array(readtable("FNM\LAB2\output_phonons\frequencies100.txt"));
c_100 = table2array(readtable("FNM\LAB2\output_phonons\c100.txt"));

%branch indices
b1_index = 1:3:length(frequencies_100);
b2_index = 2:3:length(frequencies_100);
b3_index = 3:3:length(frequencies_100);

freq_b1_100 = frequencies_100(b1_index);
freq_b2_100 = frequencies_100(b2_index);
freq_b3_100 = frequencies_100(b3_index);

close all

figure(1)
plot(c_100,freq_b1_100, LineWidth=2, Color='Blue')
hold on
plot(c_100, freq_b2_100, LineWidth=2)
plot(c_100, freq_b3_100, LineWidth=2, Color='Blue')
set(gca, "LineWidth", 2)
box on 
grid on
set(gca, "FontSize", 18)
title("q = (1 0 0)", FontSize=32)
xlabel("c", FontSize=28)
ylabel("Hz", FontSize=28)
legend("\omega_{1}/\omega_{2}", "\omega_{3}", FontSize = 26)

%%
frequencies_110 = table2array(readtable("FNM\LAB2\output_phonons\frequencies110.txt"));
c_110 = table2array(readtable("FNM\LAB2\output_phonons\c110.txt"));

freq_b1_110 = frequencies_110(b1_index);
freq_b2_110 = frequencies_110(b2_index);
freq_b3_110 = frequencies_110(b3_index);


figure(2)
plot(c_110, freq_b1_110, LineWidth=2)
hold on
plot(c_110, freq_b2_110, LineWidth=2)
plot(c_110, freq_b3_110, LineWidth=2)
set(gca, "LineWidth", 2)
box on 
grid on
set(gca, "FontSize", 18)
title("q = (1 1 0)", FontSize=32)
xlabel("c", FontSize=28)
ylabel("Hz", FontSize=28)
legend("\omega_{1}", "\omega_{2}", "\omega_{3}", FontSize = 26)

%%
frequencies_111 = table2array(readtable("FNM\LAB2\output_phonons\frequencies111.txt"));
c_111 = table2array(readtable("FNM\LAB2\output_phonons\c111.txt"));

freq_b1_111 = frequencies_111(b1_index);
freq_b2_111 = frequencies_111(b2_index);
freq_b3_111 = frequencies_111(b3_index);

figure(3)
plot(c_111, freq_b1_111, LineWidth=2, Color="Blue")
hold on
plot(c_111, freq_b2_111, LineWidth=2)
plot(c_111, freq_b3_111, LineWidth=2, Color="Blue")
set(gca, "LineWidth", 2)
box on 
grid on
set(gca, "FontSize", 18)
title("q = (1 1 1)", FontSize=32)
xlabel("c", FontSize=28)
ylabel("Hz", FontSize=28)
legend("\omega_{1}/\omega_{2}", "\omega_{3}", FontSize = 26)
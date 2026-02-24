clear all
close all
cv_Ar = table2array(readtable("FNM\LAB2\output_phonons\cv_Ar.txt"));
T_Ar = table2array(readtable("FNM\LAB2\output_phonons\T_Ar.txt"));

cv_Kr = table2array(readtable("FNM\LAB2\output_phonons\cv_Kr.txt"));
T_Kr = table2array(readtable("FNM\LAB2\output_phonons\T_Kr.txt"));

cv_Xe = table2array(readtable("FNM\LAB2\output_phonons\cv_Xe.txt"));
T_Xe = table2array(readtable("FNM\LAB2\output_phonons\T_Xe.txt"));

cv_Ne = table2array(readtable("FNM\LAB2\output_phonons\cv_Ne.txt"));
T_Ne = table2array(readtable("FNM\LAB2\output_phonons\T_Ne.txt"));
%% Plot
figure(1)
plot(T_Ar, cv_Ar, LineWidth=2)
hold on
plot(T_Kr, cv_Kr, LineWidth=2)
plot(T_Xe, cv_Xe, LineWidth=2)
plot(T_Ne, cv_Ne, LineWidth=2)
set(gca, "LineWidth", 2)
box on 
grid on
set(gca, "FontSize", 18)
title("", FontSize=32)
xlabel("T(K)", FontSize=28)
ylabel("C_{V}/V", FontSize=28)
legend("Argon", "Krypton", "Xenon", "Neon", Orientation = 'horizontal', FontSize = 26)
ylim([0 2.5e6])
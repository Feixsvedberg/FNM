close all
clear all

V_morse = table2array(readtable("FNM\LAB4\output_numerov\V_morse.txt"));



psi_0 = table2array(readtable("FNM\LAB4\output_numerov\psi_E0_morse.txt"));
psi_1 = table2array(readtable("FNM\LAB4\output_numerov\psi_E1_morse.txt"));
psi_2 = table2array(readtable("FNM\LAB4\output_numerov\psi_E2_morse.txt"));
psi_3 = table2array(readtable("FNM\LAB4\output_numerov\psi_E3_morse.txt"));
psi_4 = table2array(readtable("FNM\LAB4\output_numerov\psi_E4_morse.txt"));
psi_5 = table2array(readtable("FNM\LAB4\output_numerov\psi_E5_morse.txt"));

R_morse =table2array(readtable("FNM\LAB4\output_numerov\R_morse.txt"));

N = length(psi_0);
R0 = 0.75e-10;
RN = 3e-10;

grid = linspace(R0, RN, N);

%%
figure(1)
tiledlayout(3,2, TileSpacing="compact")
nexttile(1)
plot(grid, psi_0)

nexttile(2)
plot(grid, psi_1)

nexttile(3)
plot(grid, psi_2)

nexttile(4)
plot(grid, psi_3)

nexttile(5)
plot(grid, psi_4)

nexttile(6)
plot(grid, psi_5)




figure(2)
plot(grid, V_morse)




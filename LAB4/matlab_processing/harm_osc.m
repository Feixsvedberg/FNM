close all
clear all
Y_left = table2array(readtable('FNM\LAB4\output_numerov\Y_left.txt'));
Y_right = table2array(readtable('FNM\LAB4\output_numerov\Y_right.txt'));
psi_E0 = table2array(readtable('FNM\LAB4\output_numerov\psi_E0.txt'));
psi_E0_exact = table2array(readtable('FNM\LAB4\output_numerov\psi_E0_exact.txt'));
psi_E0_diff = table2array(readtable('FNM\LAB4\output_numerov\psi_E0_diff.txt'));

psi_E1 = table2array(readtable('FNM\LAB4\output_numerov\psi_E1.txt'));
psi_E1_exact = table2array(readtable('FNM\LAB4\output_numerov\psi_E1_exact.txt'));
psi_E1_diff = table2array(readtable('FNM\LAB4\output_numerov\psi_E1_diff.txt'));

psi_E2 = table2array(readtable('FNM\LAB4\output_numerov\psi_E2.txt'));
psi_E2_exact = table2array(readtable('FNM\LAB4\output_numerov\psi_E2_exact.txt'));
psi_E2_diff= table2array(readtable('FNM\LAB4\output_numerov\psi_E2_diff.txt'));

psi_E3 = table2array(readtable('FNM\LAB4\output_numerov\psi_E3.txt'));
psi_E3_exact = table2array(readtable('FNM\LAB4\output_numerov\psi_E3_exact.txt'));
psi_E3_diff= table2array(readtable('FNM\LAB4\output_numerov\psi_E3_diff.txt'));
%% 
close all
N = 1000;
left_grid = linspace(-8, 0, N/2);
right_grid = linspace(0, 8, N/2);
grid = linspace(-8, 8, N);

% figure(1)
% plot(left_grid, Y_left)
% hold on
% plot(right_grid, Y_right)

figure(1)
tiledlayout(2,2,TileSpacing="compact")

nexttile(1)
plot(grid, psi_E0);
hold on 
plot(grid, psi_E0_exact, Marker="v")
legend("psi", "psi_exact")

nexttile(2)
plot(grid, psi_E1, Marker="x",Color='blue')
hold on
plot(grid, psi_E1_exact, Marker="v", Color='red')
legend("psi", "psi_exact")

nexttile(3)
plot(grid, psi_E2, Marker="x",Color='blue')
hold on
plot(grid, psi_E2_exact, Marker="v", Color='red')
legend("psi", "psi_exact")

nexttile(4)
plot(grid, psi_E3, Marker="x",Color='blue')
hold on
plot(grid, psi_E3_exact, Marker="v", Color='red')
legend("psi", "psi_exact")

figure(2)
tiledlayout(2,2,TileSpacing="compact")
nexttile(1)
plot(grid, psi_E0_diff)

nexttile(2)
plot(grid, psi_E1_diff)

nexttile(3)
plot(grid, psi_E2_diff)

nexttile(4)
plot(grid, psi_E3_diff)

close all
clear all
Y_left = table2array(readtable('FNM\LAB4\output_numerov\Y_left.txt'));
Y_right = table2array(readtable('FNM\LAB4\output_numerov\Y_right.txt'));
psi_E0 = table2array(readtable('FNM\LAB4\output_numerov\psi_E0.txt'));
psi_E0_exact = table2array(readtable('FNM\LAB4\output_numerov\psi_E0_exact.txt'));
psi_E0_diff = table2array(readtable('FNM\LAB4\output_numerov\psi_E0_diff.txt'));

psi_E0_bc = table2array(readtable('FNM\LAB4\output_numerov\psi_E0_bc.txt'));

psi_E1 = table2array(readtable('FNM\LAB4\output_numerov\psi_E1.txt'));
psi_E1_exact = table2array(readtable('FNM\LAB4\output_numerov\psi_E1_exact.txt'));
psi_E1_diff = table2array(readtable('FNM\LAB4\output_numerov\psi_E1_diff.txt'));

psi_E2 = table2array(readtable('FNM\LAB4\output_numerov\psi_E2.txt'));
psi_E2_exact = table2array(readtable('FNM\LAB4\output_numerov\psi_E2_exact.txt'));
psi_E2_diff= table2array(readtable('FNM\LAB4\output_numerov\psi_E2_diff.txt'));

psi_E3 = table2array(readtable('FNM\LAB4\output_numerov\psi_E3.txt'));
psi_E3_exact = table2array(readtable('FNM\LAB4\output_numerov\psi_E3_exact.txt'));
psi_E3_diff= table2array(readtable('FNM\LAB4\output_numerov\psi_E3_diff.txt'));

psi_all_left = table2array(readtable('FNM\LAB4\output_numerov\psi_all_left.txt'));
%% 
close all
N = length(psi_E0);
left_grid = linspace(-8, 0, N/2);
right_grid = linspace(0, 8, N/2);
grid = linspace(-8, 10, N);

% figure(1)
% plot(left_grid, Y_left)
% hold on
% plot(right_grid, Y_right)

figure(1)

t = tiledlayout(4,1,TileSpacing="compact")
title(t, "a)", FontSize = 24)
nexttile(1)

plot(grid, psi_E0, Color='blue', LineWidth = 2);
%hold on 
%plot(grid, psi_E0_exact, Color = 'red')
legend("\psi_{0}")
box on
grid on
set(gca, "FontSize", 12)
ylim([-0.2 1])
legend("\psi_{0}", FontSize = 18)


nexttile(2)
plot(grid, psi_E1, Color='blue', LineWidth = 2);
%hold on
%plot(grid, psi_E1_exact, Color='red')
legend("\psi_{1}", FontSize = 18)
box on
grid on
set(gca, "FontSize", 12)
ylim([-0.8 0.8])

nexttile(3)
plot(grid, psi_E2, Color='blue', LineWidth = 2);
%hold on
%plot(grid, psi_E2_exact, Color='red')
legend("\psi_{2}")
box on
grid on
set(gca, "FontSize", 12)
ylim([-0.8 0.8])
legend("\psi_{2}", FontSize = 18)

nexttile(4)
plot(grid, psi_E3, Color='blue', LineWidth = 2);
%hold on
%plot(grid, psi_E3_exact, Color='red')
legend("\psi_{3}")
box on
grid on
set(gca, "FontSize", 12)
%ylim([-0.8 0.8])
legend("\psi_{3}", FontSize = 18)

%set(gcf,'PaperPositionMode','auto');        % preserves figure size
%print(gcf,'C:\Dokument\4an\FNM\plots_numerov\norm_psi','-depsc2')

figure(2)

t = tiledlayout(4,1,TileSpacing="compact")
title(t, "b)", FontSize = 24)
nexttile(1)

plot(grid, psi_E0_diff, LineWidth=2)
box on
grid on
set(gca, "FontSize", 12)
%ylim([-6 10]*1e-13)
legend("n = 0", FontSize = 18)




nexttile(2)
plot(grid, psi_E1_diff, LineWidth=2)
box on
grid on
set(gca, "FontSize", 12)
%ylim([-6 6]*1e-13)
legend("n = 1", FontSize = 18)



nexttile(3)
plot(grid, psi_E2_diff, LineWidth=2)
box on
grid on
set(gca, "FontSize", 12)
%ylim([-10 2]*1e-13)
legend("n = 2", FontSize = 18)

nexttile(4)
plot(grid, psi_E3_diff, LineWidth=2)
box on
grid on
set(gca, "FontSize", 12)
legend("n = 3", FontSize = 18)
%set(gcf,'PaperPositionMode','auto');        % preserves figure size
%print(gcf,'C:\Dokument\4an\FNM\plots_numerov\diff_psi','-depsc2')

figure(3)
plot(grid, psi_E0-psi_E0_bc, LineWidth=2)
grid on
box on
set(gca, "LineWidth", 2)
set(gca, "FontSize", 16)
xlim([-8 8])
%set(gcf,'PaperPositionMode','auto');        % preserves figure size
%print(gcf,'C:\Dokument\4an\FNM\plots_numerov\different_ini','-depsc2')


f = @(x) exp(x.^2/2);
anal_unphys = f(grid);
figure(4)
plot(grid, log(abs(psi_all_left)), LineWidth=2)
hold on
plot(grid, log(anal_unphys/(norm(anal_unphys))), LineWidth=2) 
grid on
box on
set(gca, "FontSize", 16)
set(gca, "LineWidth", 2)
xlim([5 8])
legend("ln(\psi_{left})", "ln(exp(x^{2}/2)", fontsize=22, Location="southeast")
%set(gcf,'PaperPositionMode','auto');        % preserves figure size
%print(gcf,'C:\Dokument\4an\FNM\plots_numerov\log','-depsc2')


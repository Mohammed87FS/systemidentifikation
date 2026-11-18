%% Bilder fuers README
%  Voraussetzung: aufgabe_type1_stoerung_entfernen.m wurde einmal ausgefuehrt
%  (type1_bereinigt_workspace.mat liegt im Ordner)

clear; clc; close all;

thisDir = fileparts(mfilename('fullpath'));
if isempty(thisDir), thisDir = pwd; end
L = load(fullfile(thisDir, 'type1_bereinigt_workspace.mat'), 'Ts', 'u_clean', 'y_clean');
outDir = fullfile(thisDir, 'docs');
if ~exist(outDir, 'dir'), mkdir(outDir); end

z   = iddata(L.y_clean(:), L.u_clean(:), L.Ts);
sys = tfest(z, 2, 1, 0, tfestOptions('Display', 'off'));

% Ausschnitt Messung vs. Modell, ueber 10 s sieht man nichts
t    = (0:numel(L.y_clean)-1)' * L.Ts;
ysim = lsim(idtf(sys), L.u_clean(:), t);
idx  = t >= 2 & t <= 3;

figure('Color', 'w', 'Position', [100 100 900 360]);
plot(t(idx), L.y_clean(idx), 'k-', t(idx), ysim(idx), 'r--', 'LineWidth', 1.1);
grid on
xlabel('t (s)'); ylabel('y');
legend('y bereinigt (Messung)', 'Modell G(s), simuliert', 'Location', 'best');
title('Ausschnitt 2 s bis 3 s, Fit ueber alle 10 s: 96,95 %');
exportgraphics(gcf, fullfile(outDir, 'modell_vs_messung.png'), 'Resolution', 130);

% Sprungantwort und Pole/Nullstellen vom identifizierten Modell
figure('Color', 'w', 'Position', [100 100 900 360]);
subplot(1, 2, 1);
step(sys, 0.08); grid on
title('Sprungantwort G(s)');
subplot(1, 2, 2);
pzmap(sys);
xlim([-120 10]); ylim([-160 160]);
title('Pole und Nullstelle');
exportgraphics(gcf, fullfile(outDir, 'sprung_pzmap.png'), 'Resolution', 130);

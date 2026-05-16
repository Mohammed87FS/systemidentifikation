%% Typ 1 — Systemidentifikation auf bereinigten Daten
%  Voraussetzung: type1_bereinigt_workspace.mat im gleichen Ordner.

clear; clc; close all;

%% --- Daten laden ---
thisDir = fileparts(mfilename('fullpath'));
if isempty(thisDir), thisDir = pwd; end
L = load(fullfile(thisDir, 'type1_bereinigt_workspace.mat'), ...
         'Ts', 'u_clean', 'y_clean');

z = iddata(L.y_clean(:), L.u_clean(:), L.Ts);

%% --- Ordnungs-Sweep: RSS nach 1-Schritt-Prädiktion ---
pairs = [2 1; 2 2; 3 2];   % [np nz]
opt   = tfestOptions('Display', 'off');

fprintf('\n=== tfest Ordnungsvergleich (kleinere RSS = besser) ===\n');
fprintf('%-4s %-4s %-12s\n', 'np', 'nz', 'RSS_pred');
y = z.OutputData;
for r = 1:size(pairs, 1)
    np = pairs(r, 1);  nz = pairs(r, 2);
    sys = tfest(z, np, nz, 0, opt);
    yhat = predict(sys, z).OutputData;
    rss  = sum((y - yhat).^2);
    fprintf('%-4d %-4d %-12.4g\n', np, nz, rss);
end

%% --- Hauptmodell: np=2, nz=1 (Parsimonie-Wahl) ---
np = 2;  nz = 1;
sys = tfest(z, np, nz, 0, opt);

% Übertragungsfunktion anzeigen
fprintf('\n--- Geschätzte Übertragungsfunktion ---\n');
sys

% Pole/Nullstellen mit omega_n und zeta
fprintf('\n--- Pole ---\n');
p = pole(sys);
for k = 1:numel(p)
    if imag(p(k)) > 0
        wn   = abs(p(k));
        zeta = -real(p(k)) / wn;
        fprintf('  s = %.3g %+.3gj  -> omega_n=%.2f rad/s (%.2f Hz), zeta=%.4f\n', ...
                real(p(k)), imag(p(k)), wn, wn/(2*pi), zeta);
    end
end
fprintf('\n--- Nullstelle(n) ---\n');
disp(zero(sys));

%% --- Plots ---
figure('Name', 'compare', 'Color', 'w');
compare(z, sys);

figure('Name', 'Pole/Nullstellen', 'Color', 'w');
pzmap(sys); grid on;
title(sprintf('tfest np=%d, nz=%d', np, nz));
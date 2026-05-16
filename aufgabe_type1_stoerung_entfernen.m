%% Typ 1 — Überlagerung entfernen: u(t), y(t) aus 2026_Type1_nr1.mat
%  f1(t) = A*cos(w*t + phi),  w = 4 rad/s
%  f2(t) = a*exp(-t^2/b)
%  f3(t) = a*ln(sqrt(t) + 0.7)
%  Methode 1: KQ / Pseudoinverse  |  Methode 2: fminsearch
%  Modellwahl: kleinste RSS

clear; clc; close all;

%% --- Daten laden ---
thisDir = fileparts(mfilename('fullpath'));
if isempty(thisDir), thisDir = pwd; end
S = load(fullfile(thisDir, '2026_Type1_nr1.mat'));
t = S.t(:);  u = S.u(:);  y = S.y(:);
Ts = mean(diff(t));
w  = 4;   % rad/s
fprintf('N=%d, Ts=%.4g s\n', numel(t), Ts);

%% --- Modelldefinitionen (Prediktor-Funktionen) ---
% c bzw. p = Parametervektor, tt = Zeitvektor
f1 = @(c, tt) c(1)*cos(w*tt) + c(2)*sin(w*tt);          % linear in c
f2 = @(p, tt) p(1)*exp(-(tt.^2)/p(2));                   % nicht-linear in p(2)
f3 = @(a, tt) a*log(sqrt(max(tt,0)) + 0.7);              % linear in a

%% --- Fits für beide Signale (u und y), beide Methoden ---
results = struct();
for sig = {'u', 'y'}
    s = sig{1};
    d = eval(s);   % Messsignal

    % METHODE 1: KQ / Pseudoinverse
    % f1: linear, Regressor [cos(wt) sin(wt)]
    X = [cos(w*t), sin(w*t)];
    c_lin = pinv(X) * d;
    results.(s).m1.f1 = make_fit(f1, c_lin, t, d);

    % f3: linear, eine Spalte
    X3 = log(sqrt(max(t,0)) + 0.7);
    a_lin = pinv(X3) * d;
    results.(s).m1.f3 = make_fit(f3, a_lin, t, d);

    % f2: b-Raster + KQ für a (nichtlineares b, lineares a)
    b_grid = logspace(log10(Ts^2), log10((t(end)-t(1))^2), 400);
    bestRSS = inf;  bestParams = [];  bestPred = [];
    for b = b_grid
        xe = exp(-(t.^2)/b);
        a  = pinv(xe) * d;
        pr = a * xe;
        r  = sum((d - pr).^2);
        if r < bestRSS
            bestRSS = r;  bestParams = [a; b];  bestPred = pr;
        end
    end
    results.(s).m1.f2 = struct('pred', bestPred, 'params', bestParams, 'RSS', bestRSS);

    % METHODE 2: fminsearch auf RSS
    % f2 in log(b) parametrisieren, damit b > 0 garantiert
    f2_log = @(p, tt) p(1)*exp(-(tt.^2)/exp(p(2)));

    [c_opt, rss1] = fminsearch(@(c) sum((d - f1(c,t)).^2),       c_lin);
    [p_opt, rss2] = fminsearch(@(p) sum((d - f2_log(p,t)).^2),   [results.(s).m1.f2.params(1); log(results.(s).m1.f2.params(2))]);
    [a_opt, rss3] = fminsearch(@(a) sum((d - f3(a,t)).^2),       a_lin);

    results.(s).m2.f1 = make_fit(f1, c_opt, t, d);
    results.(s).m2.f2 = struct('pred', f2_log(p_opt,t), 'params', [p_opt(1); exp(p_opt(2))], 'RSS', rss2);
    results.(s).m2.f3 = make_fit(f3, a_opt, t, d);
end

%% --- RSS-Tabellen ---
fprintf('\n=== RSS Methode 1 (KQ) ===\n');
print_rss(results, 'm1');
fprintf('\n=== RSS Methode 2 (fminsearch) ===\n');
print_rss(results, 'm2');

%% --- Modellwahl pro Kanal (kleinstes RSS aus Methode 2) ---
models = {'f1', 'f2', 'f3'};
[u_hat, iu, ~] = pick_best(results.u.m2, models, t);
[y_hat, iy, ~] = pick_best(results.y.m2, models, t);
fprintf('\nGewählt: u → %s, y → %s\n', models{iu}, models{iy});

u_clean = u - u_hat;
y_clean = y - y_hat;

%% --- Plots ---
plot_signal('Eingang u', t, u, u_hat, u_clean, models{iu});
plot_signal('Ausgang y', t, y, y_hat, y_clean, models{iy});

% Methodenvergleich für f1 auf u (KQ vs fminsearch)
figure('Name', 'f1 auf u: beide Methoden', 'Color', 'w');
plot(t, u, 'k-'); hold on;
plot(t, results.u.m1.f1.pred, 'b--', 'LineWidth', 1.2);
plot(t, results.u.m2.f1.pred, 'r-',  'LineWidth', 1.0);
grid on; legend('u', 'f1 KQ', 'f1 fminsearch', 'Location', 'best');
xlabel('t (s)'); title('Methodenvergleich');

%% --- Speichern ---
fit_methode1 = collect_method(results, 'm1');
fit_methode2 = collect_method(results, 'm2');
save(fullfile(thisDir, 'type1_bereinigt_workspace.mat'), ...
     't', 'Ts', 'u', 'y', 'u_clean', 'y_clean', 'u_hat', 'y_hat', ...
     'fit_methode1', 'fit_methode2', 'iu', 'iy', 'models');
fprintf('\nGespeichert: type1_bereinigt_workspace.mat\n');

%% ===== Helferfunktionen =====
function fit = make_fit(predFun, params, t, d)
    pr = predFun(params, t);
    fit = struct('pred', pr, 'params', params, 'RSS', sum((d - pr).^2));
end

function print_rss(R, method)
    for sig = {'u', 'y'}
        s = sig{1};
        fprintf('Kanal %s: f1=%.4g  f2=%.4g  f3=%.4g\n', s, ...
            R.(s).(method).f1.RSS, R.(s).(method).f2.RSS, R.(s).(method).f3.RSS);
    end
end

function [best_pred, idx, best_rss] = pick_best(method_struct, models, ~)
    rss = cellfun(@(m) method_struct.(m).RSS, models);
    [best_rss, idx] = min(rss);
    best_pred = method_struct.(models{idx}).pred;
end

function plot_signal(name, t, raw, hat, clean, modelName)
    figure('Name', name, 'Color', 'w');
    plot(t, raw, 'k-'); hold on;
    plot(t, hat,   'r-', 'LineWidth', 1.2);
    plot(t, clean, 'b-');
    grid on; legend('Messung', sprintf('Überlagerung (%s)', modelName), 'bereinigt', 'Location', 'best');
    xlabel('t (s)'); title(name);
end

function out = collect_method(R, method)
    out = struct('u', R.u.(method), 'y', R.y.(method));
end
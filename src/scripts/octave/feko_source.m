close all
clear
clc

%%
% 光速m/s
C0 = 299792458;
MUE0 = 4e-7 * pi;
EPS0 = 1 / (MUE0 * C0 ^ 2);
Z0 = sqrt(MUE0 / EPS0);
% 点数
N = 401;
% 频率/Hz
freq = 94.05e9;
% 波长/m
wavelength = C0 / freq;
% 束腰半径/m
w0 = 4e-3;
% 采样间隔/m
ds = 1.5e-3;
% 计算区域尺寸
length = ds * (N - 1);

% 束腰面
z0 = wavelength * 0;
Fdata = gauss_source(freq, w0, ds, N, z0);
AG = 20 * log10(abs(Fdata) + eps);
Fdata = reshape(Fdata, N * N, 1);

F_abs = abs(Fdata);
F_arg = angle(Fdata) * 180 / pi; % 角度要用角度制

% Ey
Edata = zeros(N * N, 4);
Edata(:, 3) = F_abs;
Edata(:, 4) = F_arg;
save('E_field.txt', 'Edata', '-ascii');

% Hx
Hdata = zeros(N * N, 4);
Hdata(:, 1) = F_abs / Z0;
Hdata(:, 2) = F_arg;
save('H_field.txt', 'Hdata', '-ascii');

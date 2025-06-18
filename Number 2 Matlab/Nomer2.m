% TUGAS TAKE-HOME EXAM - WF2202 - SOAL 2
% NAMA: [Muhamad Hanif Hafizhan] 13123069
%       [Mochamad Arkan Nugraha] 13123007

clear all;
clc;
close all;

L = 1;
alpha = 0.1;
Ti = 100; %T initial
Ts = 300; %T surface
dx = 0.05;
t_final = 2; %asumsi saja untuk bisa diperlihatkan di plotnya

dt_array = [0.005, 0.01, 0.05]; 


%Solusi
num_cases = length(dt_array);
x = (0:dx:L)';
T_final_profiles = zeros(length(x), num_cases);

for k = 1:num_cases
    dt = dt_array(k);
    lambda = alpha * dt / dx^2; %Mengecek lambda karena jika lebih dari 0.5 visualisasi dari numeriknya gagal

    if lambda <= 0.5, fprintf('Hasil dt = %.3f: STABIL \n\n', dt); 
    else, fprintf('Hasil dt = %.3f: TIDAK STABIL \n\n', dt); 
    end

    T = ones(length(x), 1) * Ti; %buat vektor kolom yang isinya 1 semua
    T(1) = Ts; %membuat ujung awal dan akhir T nya sama sama 300
    T(end) = Ts;

    T_old = T;
    t_steps = round(t_final / dt);
    for p = 1:t_steps
        for i = 2:length(x)-1
            T(i) = T_old(i) + lambda * (T_old(i+1) - 2*T_old(i) + T_old(i-1)); %Turunan dari PDE heat equation menggunakan FTCS
        end
        T_old = T;
    end
    T_final_profiles(:, k) = T;
end


%PLOT nya (ada dua yang kasus stabil dan tidak stabil)
figure('Name', 'Analisis Hasil Simulasi Lengkap', 'NumberTitle', 'off', 'Position', [100, 100, 800, 700]);


subplot(2, 1, 1);
hold on;

% Plot kasus stabil (dt=0.005)
plot(x, T_final_profiles(:, 1), 'b-o', 'LineWidth', 1, 'MarkerSize', 5, 'DisplayName', 'dt=0.005 (λ=0.2)');
% Plot kasus stabil (dt=0.01)
plot(x, T_final_profiles(:, 2), 'g--x', 'LineWidth', 1, 'MarkerSize', 7, 'DisplayName', 'dt=0.01 (λ=0.4)');

% Hitung dan plot solusi analitis (sudah dihitung secara manual bag. a)
n_terms = 100;
n_vec = (1:2:2*n_terms)';
C_n = -800 ./ (pi * n_vec);
sin_terms = sin(n_vec * pi * x' / L);
exp_terms = exp(-alpha * (n_vec * pi / L).^2 * t_final);
transient_sum = sum(C_n .* sin_terms .* exp_terms, 1);
T_analytical = Ts + transient_sum';
plot(x, T_analytical, 'k-', 'LineWidth', 3, 'DisplayName', 'Solusi Analitis');

% Pengaturan Plot Atas
title('Panel A: Perbandingan Solusi Stabil');
ylabel('Suhu (T) [°F]');
legend('show', 'Location', 'southeast', 'FontSize', 11);
grid on;
box on;
axis([0 L 100 310]);


subplot(2, 1, 2);
% Plot kasus tidak stabil (dt=0.05)
plot(x, T_final_profiles(:, 3), 'r:d', 'LineWidth', 1, 'MarkerSize', 5, 'DisplayName', 'dt=0.05 (λ=2.0)'); % 'd' untuk diamond

% Pengaturan Plot Bawah
title('Panel B: Visualisasi Osilasi Liar (Tidak Stabil)');
xlabel('Posisi (x) [ft]');
ylabel('Suhu (T) [°F]');
legend('show', 'Location', 'northeast', 'FontSize', 11);
grid on;
box on;

hold off;
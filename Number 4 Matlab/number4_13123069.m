% TUGAS TAKE-HOME EXAM - WF2202 - SOAL 3
% NAMA: [Mochaman Arkan Nugraha] 13123007
%       [Muhamad Hanif Hafizhan] 13123069

clear;clc;close all;

% Inisalisasi parameter (13123091)
E = 23069;      
M = 3069;       
I = 1500;       
L = 200;        
h = 20;         

% konstanta C untuk menyederhanakan persamaan
C = M / (E * I);

%ini function handle jadi saya bisa langsung recall fy dan fz nya langsung
fy = @(z) z; 
fz = @(z, C) C * (1 + z^2)^(3/2); 

%inisialisasi interval 0 sampai L dengan step h
interval = 0:h:L;

n = length(interval); %berapa kali harus melangkah


%Membuat matrix 0 untuk menyimpan hasil dan membuat kondisi awalnya
y_rk2 = zeros(1, n);
z_rk2 = zeros(1, n);

y_rk4 = zeros(1, n);
z_rk4 = zeros(1, n);


%Perhitungan prediktor korektornya
for i = 1:(n - 1)
    %RK orde 2 sesuai formula bagian a
    k1_rk2 = h * fy(z_rk2(i));
    l1_rk2 = h * fz(z_rk2(i), C);
    
    k2_rk2 = h * fy(z_rk2(i) + l1_rk2);
    l2_rk2 = h * fz(z_rk2(i) + l1_rk2, C);
    
    y_rk2(i+1) = y_rk2(i) + 0.5 * (k1_rk2 + k2_rk2);
    z_rk2(i+1) = z_rk2(i) + 0.5 * (l1_rk2 + l2_rk2);
    
    %RK orde 4
    k1_rk4 = h * fy(z_rk4(i));
    l1_rk4 = h * fz(z_rk4(i), C);
    
    k2_rk4 = h * fy(z_rk4(i) + 0.5*l1_rk4);
    l2_rk4 = h * fz(z_rk4(i) + 0.5*l1_rk4, C);
    
    k3_rk4 = h * fy(z_rk4(i) + 0.5*l2_rk4);
    l3_rk4 = h * fz(z_rk4(i) + 0.5*l2_rk4, C);
    
    k4_rk4 = h * fy(z_rk4(i) + l3_rk4);
    l4_rk4 = h * fz(z_rk4(i) + l3_rk4, C);
    
    y_rk4(i+1) = y_rk4(i) + (k1_rk4 + 2*k2_rk4 + 2*k3_rk4 + k4_rk4) / 6;
    z_rk4(i+1) = z_rk4(i) + (l1_rk4 + 2*l2_rk4 + 2*l3_rk4 + l4_rk4) / 6;
    
end



fprintf('\n\n--- Nilai y1 dan z1 dengan RK Orde 2 ---\n');
fprintf('   Pada x = %.f mm:\n', interval(2));
fprintf('   Defleksi (y1) = %.6f mm\n', y_rk2(2));
fprintf('   Kemiringan (z1) = %.6f\n', z_rk2(2));
fprintf('\n\n--- Tabel Hasil y dan z dengan RK Orde 4 ---\n');
hasil_tabel_rk4 = table(interval', y_rk4', z_rk4', 'VariableNames', {'x_mm', 'y_defleksi_mm', 'z_kemiringan'});
disp(hasil_tabel_rk4);

%PLOT nya
figure;
hold on;

% Plot RK4
plot(interval, y_rk4, '-', 'LineWidth', 2.5, 'DisplayName', 'RK Orde 4');
% Plot RK2
plot(interval, y_rk2, '--o', 'LineWidth', 1, 'Color', 'g', 'MarkerSize', 6, 'DisplayName', 'RK Orde 2');

hold off;
grid on;
title('Perbandingan Defleksi Baloks');
xlabel('Jarak dari Ujung Jepit, x (mm)');
ylabel('Defleksi, y (mm)');
legend('show', 'Location', 'northwest');
set(gca, 'FontSize', 12);
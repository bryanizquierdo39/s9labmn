%% Laboratorio: Trabajo de compresión W = integral P(V) dV
clc; clear;

% Datos del ejercicio: n = 12 subintervalos, 13 puntos, h = 0.25 L
V = [1, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 2.75, 3, 3.25, 3.5, 3.75, 4];
P = [300.2, 242.1, 201.5, 169.8, 151, 134.3, 120.8, 107.9, 99.1, 93.4, 86.2, 79.5, 75.1];

a = V(1); b = V(end); n = 12; h = (b-a)/n;

fprintf('Datos: V de %.2f a %.2f L, n = %d, h = %.2f L\n', a, b, n, h);

% 1. Regla del Trapecio Compuesta - Polinomio grado 1
% ---------------------------------------------------------------------------------------------------------------------------------------------
% I ≈ h/2 [f(x0) + 2*sum(f(xi)) + f(xn)]
S_trap = P(1) + P(end);
for i = 2:n
    S_trap = S_trap + 2*P(i);
end
W_trap = h/2 * S_trap;

% 2. Regla de Simpson 1/3 Compuesta - Polinomio grado 2, n debe ser
% par------------------------------------------------------------------------------------------------------------------------------------------
% I ≈ h/3 [f(x0) + 4*sum(f(impares)) + 2*sum(f(pares)) + f(xn)]
if mod(n,2) ~= 0
    error('Simpson 1/3 requiere n par');
end
S_simp13 = P(1) + P(end);

% x1,x3,x5,... => peso 4
for i = 2:2:n
    S_simp13 = S_simp13 + 4*P(i);
end

% x2,x4,x6,... => peso 2
for i = 3:2:n-1
    S_simp13 = S_simp13 + 2*P(i);
end
W_simp13 = h/3 * S_simp13;

% 3. Regla de Simpson 3/8 Compuesta - Polinomio grado 3, n múltiplo de
% 3--------------------------------------------------------------------------------------------------------------------------------------------
% I ≈ 3h/8 [f(x0) + 3*sum(i≠3k) + 2*sum(i=3k) + f(xn)]
if mod(n,3) ~= 0
    error('Simpson 3/8 requiere n múltiplo de 3');
end
S_simp38 = P(1) + P(end);

for i = 2:n

    j = i - 1;  % subíndice matemático

    if mod(j,3) == 0
        S_simp38 = S_simp38 + 2*P(i);
    else
        S_simp38 = S_simp38 + 3*P(i);
    end

end
W_simp38 = 3*h/8 * S_simp38;

% 4. Cuadratura de Gauss-Legendre m=3
% Filosofía: nodos óptimos en [-1,1] mapeados a [a,b]
% Para m=3: t = [-sqrt(3/5), 0, sqrt(3/5)], w = [5/9, 8/9, 5/9]
% Como solo tenemos datos discretos, creamos P(V) con spline cúbico
pp = spline(V, P); % interpolación cúbica como pide el ejercicio

t = [-sqrt(3/5), 0, sqrt(3/5)];
w = [5/9, 8/9, 5/9];

W_gauss = 0;
for i = 1:3
    % Mapeo lineal [t] -> [x]: x = (b-a)/2 * t + (b+a)/2
    x_i = (b-a)/2 * t(i) + (b+a)/2;
    P_i = ppval(pp, x_i); % evaluamos el spline en el nodo Gauss
    W_gauss = W_gauss + w(i) * P_i;
end
W_gauss = (b-a)/2 * W_gauss; % factor de escala del cambio de variable

% Tabla de resultados - Gauss como referencia para error relativo
W_ref = integral(@(x) ppval(pp,x),a,b);
metodos = {'Trapecio'; 'Simpson 1/3'; 'Simpson 3/8'; 'Gauss-Legendre'};
valores = [W_trap; W_simp13; W_simp38; W_gauss];
error_rel = abs((valores - W_ref)./W_ref) * 100;

fprintf('MÉTODO | VALOR APROXIMADO [kPa·L] | ERROR RELATIVO %%\n');
fprintf('----------------------------------------------------------------\n');
for i = 1:4
    fprintf('%-15s | %20.4f | %15.4f\n', metodos{i}, valores(i), error_rel(i));
end

% Gráfico para visualizar el ruido y la interpolación
figure;
plot(V, P, 'ro', 'MarkerSize', 8, 'DisplayName', 'Datos con ruido');
hold on;
V_fino = linspace(a, b, 200);
P_fino = ppval(pp, V_fino);
plot(V_fino, P_fino, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Spline cúbico P(V)');
xlabel('Volumen V [L]'); ylabel('Presión P [kPa]');
title('Interpolación cúbica para Gauss-Legendre');
legend; grid on;
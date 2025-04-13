function x = fkine(q)
% fkine calcola la cinematica diretta del manipolatore
% Input:
%   q: vettore delle variabili articolari [theta1; theta2; theta3]
% Output:
%   x: vettore di posa [x; y; z] dove:
%      x = cos(theta1)+cos(theta1+theta2)+cos(theta1+theta2+theta3)
%      y = sin(theta1)+sin(theta1+theta2)+sin(theta1+theta2+theta3)
%      z = theta1+theta2+theta3

theta1 = q(1);
theta2 = q(2);
theta3 = q(3);

x = [ cos(theta1) + cos(theta1+theta2) + cos(theta1+theta2+theta3);
      sin(theta1) + sin(theta1+theta2) + sin(theta1+theta2+theta3);
      theta1 + theta2 + theta3 ];
end


function Jinv = J_inv(q)
% J_inv calcola l'inversa della Jacobiana su q
% Input:
%   q: vettore delle variabili articolari [theta1; theta2; theta3]
% Output:
%   J: matrice Jacobiana 3x3 inversa

theta1 = q(1);
theta2 = q(2);
theta3 = q(3);

det_J = 1/sin(theta2);

Jinv = zeros(3,3);
Jinv(1,1) = cos(theta1+theta2);
Jinv(1,2) = sin(theta1+theta2);
Jinv(1,3) = sin(theta3);

Jinv(2,1) = -cos(theta1) -cos(theta1+theta2);
Jinv(2,2) = -sin(theta1) -sin(theta1+theta2);
Jinv(2,3) = -sin(theta2+theta3) -sin(theta3);

Jinv(3,1) = cos(theta1);
Jinv(3,2) = sin(theta1);
Jinv(3,3) = sin(theta2+theta3) + sin(theta2);

Jinv = det_J*Jinv;
end

function J = jacobianIK(q)
% jacobianIK calcola il Jacobiano della cinematica diretta del manipolatore
% Input:
%   q: vettore delle variabili articolari [theta1; theta2; theta3]
% Output:
%   J: matrice Jacobiana 3x3

theta1 = q(1);
theta2 = q(2);
theta3 = q(3);

J = zeros(3,3);
J(1,1) = -sin(theta1) - sin(theta1+theta2) - sin(theta1+theta2+theta3);
J(1,2) = -sin(theta1+theta2) - sin(theta1+theta2+theta3);
J(1,3) = -sin(theta1+theta2+theta3);

J(2,1) = cos(theta1) + cos(theta1+theta2) + cos(theta1+theta2+theta3);
J(2,2) = cos(theta1+theta2) + cos(theta1+theta2+theta3);
J(2,3) = cos(theta1+theta2+theta3);

J(3,:) = [1 1 1];
end



%% Gradient method
function [q, err, iter] = gradientIK(q0, x_d, alpha, tol, max_iter)
    % gradientIK esegue il metodo del gradiente per la J inversa.
    % Input:
    %   q0      - vettore iniziale di configurazione [theta1; theta2; theta3]
    %   x_d     - posa desiderata, vettore [x; y; z]
    %   alpha   - fattore di aggiornamento
    %   tol     - tolleranza per il criterio di arresto (norma dell'errore)
    %   max_iter- numero massimo di iterazioni
    % Output:
    %   q       - configurazione finale approssimata
    %   err     - errore finale (norma di x_d - fkine(q))
    %   iter    - numero di iterazioni eseguite
    
    % Preallocazione degli array per migliorare l'efficienza
    err_history = zeros(max_iter, 1);           % vettore per memorizzare gli errori
    q_history = zeros(length(q0), max_iter);      % matrice per memorizzare i valori dei 3 angoli

    q = q0;
    iter = 1;

    while iter < max_iter
        % Riavvolgi il q in un intervallo [-π, π] per evitare multipli di
        % giri interi
        q = mod(q + pi, 2*pi) - pi;

        % Calcola posa attuale
        x = fkine(q);
        % Calcola errore
        e = x_d - x;
        err = norm(e);

        % Salva i dati per i plot (utilizzando la preallocazione)
        err_history(iter) = err;
        q_history(:, iter) = q;

        % Visualizza errore ed aggiornamento (come richiesto)
        fprintf('Iterazione %d: q = [%f, %f, %f], Errore = %e\n', iter, q(1), q(2), q(3), err);

        % Verifica criterio di arresto
        if err < tol
            fprintf('Tolleranza raggiunta!!!! \n');
            break;
        end


        %AGGIORNAMENTO
        % Calcola il Jacobiano
        J = jacobianIK(q);
        % Aggiorna il vettore q
        q = q + alpha * (J') * e;

        iter = iter + 1;
    end
    
    if iter == max_iter
        fprintf('---DIVERGE!!!!---\n');
    end
    % Se l'algoritmo ha eseguito meno di max_iter iterazioni, trimmare gli array
    err_history = err_history(1:iter);
    q_history = q_history(:, 1:iter);

    % Plot errore
    figure('Name','Gradient Method_err','NumberTitle','off');
    plot(1:length(err_history), err_history, 'LineWidth', 2);
    xlabel('Iteration');
    ylabel('Error');
    title('Gradient Method - Error evolution during iterations');
    grid on;

    % Plot angoli
    figure('Name','Gradient Method_q','NumberTitle','off');
    plot(1:size(q_history, 2), q_history(1,:), '-r', 'LineWidth', 2); hold on;
    plot(1:size(q_history, 2), q_history(2,:), '-g', 'LineWidth', 2);
    plot(1:size(q_history, 2), q_history(3,:), '-b', 'LineWidth', 2);
    xlabel('Iteration');
    ylabel('q variation');
    legend('q(1)', 'q(2)', 'q(3)');
    title('Gradient Method - q angles evolution during iterations');
    grid on;
end


%% Newton Method
function [q, err, iter] = newtonIK(q0, x_d, tol, max_iter)
    % newtonIK esegue il metodo di Newton per l'inversa cinematica.
    % Input:
    %   q0      - vettore iniziale di configurazione [theta1; theta2; theta3]
    %   x_d     - posa desiderata, vettore [x; y; z]
    %   tol     - tolleranza per il criterio di arresto (norma dell'errore)
    %   max_iter- numero massimo di iterazioni
    % Output:
    %   q       - configurazione finale approssimata
    %   err     - errore finale (norma di x_d - fkine(q))
    %   iter    - numero di iterazioni eseguite
    
    q = q0;
    iter = 1;

    err_history = zeros(max_iter, 1);
    q_history = zeros(length(q0) ,max_iter);

    while iter < max_iter
        % Riavvolgi il q in un intervallo [-π, π] per evitare multipli di
        % giri interi
        q = mod(q + pi, 2*pi) - pi;

        % Calcola posa attuale
        x = fkine(q);
        % Calcola errore
        e = x_d - x;
        err = norm(e);

        % Salva err e q
        err_history(iter) = err;
        q_history(:, iter) = q;

        % Visualizza errore ed aggiornamento (come richiesto)
        fprintf('Iterazione %d: q = [%f, %f, %f], Errore = %e\n', iter, q(1), q(2), q(3), err);
        
        % Verifica criterio di arresto
        if err < tol
            fprintf('Tolleranza raggiunta \n')
            break;
        end

        %AGGIORNAMENTO
        % Calcola il Jacobiano
        Jinv = J_inv(q);
        % Aggiorna il vettore q: risolvi Jinv * delta_q = e
        delta_q = Jinv * e;  
        q = q + delta_q;
        
        if any(isnan(q))
            fprintf('Singular Point!!\n');
            fprintf('Exiting...\n');
            return;
        end

        iter = iter + 1;
    end
    
    % Verifica se abbiamo trovato una sol prima di max_iter
    if iter == max_iter
        fprintf('---DIVERGE!!!!---\n');
    end

    % Se l'algoritmo ha eseguito meno di max_iter iterazioni, trimmare gli array
    err_history = err_history(1:iter);
    q_history = q_history(:, 1:iter);

    % Plot errore
    figure('Name','Newton Method_err','NumberTitle','off');
    plot(1:length(err_history), err_history, 'LineWidth', 2);
    xlabel('Iteration');
    ylabel('Error');
    title('Newton Method - Error evolution during iterations');
    grid on;

    % Plot angoli
    figure('Name','Newton Method_q','NumberTitle','off');
    plot(1:size(q_history, 2), q_history(1,:), '-r', 'LineWidth', 2); hold on;
    plot(1:size(q_history, 2), q_history(2,:), '-g', 'LineWidth', 2);
    plot(1:size(q_history, 2), q_history(3,:), '-b', 'LineWidth', 2);
    xlabel('Iteration');
    ylabel('q variation');
    legend('q(1)', 'q(2)', 'q(3)');
    title('Newton Method - q angles evolution during iterations');
    grid on;
end

%% Initialization
% mainIK.m
clc; clear; close all;

% Parametri iniziali e definitivi
%q0 = [0; 0; 0];             % configurazione iniziale
q0 = [pi/2; pi/2; pi/2];
x_d = [2; 1; 0];            % posa desiderata

% Parametri per il metodo del gradiente
alpha = 1/10;
%alpha = 1/2;
tol = 1e-6;
max_iter = 1000;

%% Analitic Solution

fprintf('---  Analitic Solution ---\n');
syms sym_theta1 sym_theta2 sym_theta3 real;
eqns = [ cos(sym_theta1) + cos(sym_theta1+sym_theta2) + cos(sym_theta1+sym_theta2+sym_theta3) == x_d(1);
         sin(sym_theta1) + sin(sym_theta1+sym_theta2) + sin(sym_theta1+sym_theta2+sym_theta3) == x_d(2);
         sym_theta1 + sym_theta2 + sym_theta3 == x_d(3) ];

sol = solve(eqns, [sym_theta1, sym_theta2, sym_theta3]);

n = max([length(sol.sym_theta1), length(sol.sym_theta2), length(sol.sym_theta3)]);

for i = 1:1:n
    theta1_sol = sol.sym_theta1(i);
    theta2_sol = sol.sym_theta2(i);
    theta3_sol = sol.sym_theta3(i);

    fprintf('Soluzione analitica: theta 1 = %f, theta 2 = %f, theta 3 = %f \n\n', theta1_sol, theta2_sol, theta3_sol)
    fprintf('OVVERO: \n q = [%f π, %f π, %f π] \n\n\n', theta1_sol/pi, theta2_sol/pi, theta3_sol/pi)
end

%% Gradient Method

fprintf('--- Gradient Method ---\n');
[q_grad, err_grad, iter_grad] = gradientIK(q0, x_d, alpha, tol, max_iter);
fprintf('Soluzione Gradient Method: q = [%f, %f, %f] in %d iterazioni, Errore finale = %e\n\n', ...
        q_grad(1), q_grad(2), q_grad(3), iter_grad, err_grad);
fprintf('OVVERO: \n q = [%f π, %f π, %f π] \n\n\n', q_grad(1)/pi, q_grad(2)/pi, q_grad(3)/pi)


%% Newton Method
fprintf('--- Newton Method ---\n');
[q_newton, err_newton, iter_newton] = newtonIK(q0, x_d, tol, max_iter);
fprintf('Soluzione Newton Method: q = [%f, %f, %f] in %d iterazioni, Errore finale = %e\n', ...
        q_newton(1), q_newton(2), q_newton(3), iter_newton, err_newton);
fprintf('OVVERO: \n q = [%f π, %f π, %f π] \n\n\n', q_newton(1)/pi, q_newton(2)/pi, q_newton(3)/pi)

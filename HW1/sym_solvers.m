
theta1 = 0;
theta2 = 0;
theta3 = 0;

% Definizione di Jinv
Jinv = [ cos(theta1+theta2),         sin(theta1+theta2),        sin(theta3);
        -cos(theta1)-cos(theta1+theta2), -sin(theta1)-sin(theta1+theta2), -sin(theta2+theta3)- sin(theta3);
         cos(theta1),                sin(theta1),               sin(theta2+theta3)+ sin(theta2) ];

Jinv = 1/sin(theta2) * Jinv
% Aggiorna il vettore q: risolvi Jinv * delta_q = e
delta_q = Jinv * [2;3; 0.2];  
q = [0;0;0] + delta_q

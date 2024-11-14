clear global; clc;
% Introducing a function "elestiff_key" to calculate stiffness matrix (6*6)
function key = elestiff_key(E, A, I, x)
    lxx = x(3) - x(1);
    lyy = x(4) - x(2);
    lee = sqrt(lxx^2 + lyy^2);
    C = lxx / lee; 
    S = lyy / lee; 
    CS = C * S;
    key = [(A*E/lee)*C^2,  (A*E/lee)*CS,      -(6*E*I/lee^2)*S,  -(A*E/lee)*C^2,   -(A*E/lee)*CS,       -(6*E*I/lee^2)*S;
          (A*E/lee)*CS,    (12*E*I/lee^3)*C^2, (6*E*I/lee^2)*C,  -(A*E/lee)*CS,    -(12*E*I/lee^3)*C^2, (6*E*I/lee^2)*C;
         -(6*E*I/lee^2)*S, (6*E*I/lee^2)*C,    (4*E*I/lee),       (6*E*I/lee^2)*S, -(6*E*I/lee^2)*C,    (2*E*I/lee);
         -(A*E/lee) * C^2,-(A*E/lee)*CS,       (6*E*I/lee^2)*S,   (A*E/lee)*C^2,   (A*E/lee)*CS,        (6*E*I/lee^2)*S;
         -(A*E/lee)*CS,   -(12*E*I/lee^3)*C^2, -(6*E*I/lee^2)*C,   (A*E/lee)*CS,   (12*E*I/lee^3)*C^2,  -(6*E*I/lee^2)*C;
         -(6*E*I/lee^2)*S, (6*E*I/lee^2)*C,    (2*E*I/lee),      (6*E*I/lee^2)*S,  -(6*E*I/lee^2)*C,    (4*E*I/lee)];   
end
clear global; clc;
% Introducing a function "elestiff_kee" to calculate stiffness matrix (6*6)
function kee = elestiff_kee(E, A, I, x)
    lxx = x(3) - x(1);
    lyy = x(4) - x(2);
    lee = sqrt(lxx^2 + lyy^2);
    C = lxx / lee; 
    S = lyy / lee; 
    CS = C * S;
    kee = [(12*E*I/lee^3)*S^2,   -(12*E*I/lee^3)*CS,  -(6*E*I/lee^2)*S, -(12*E*I/lee^3)*S^2,  (12*E*I/lee^3)*CS,   -(6*E*I/lee^2)*S;
          -(12*E*I/lee^3)*CS,    (A*E/lee)*S^2,       (6*E*I/lee^2)*C,  (12*E*I/lee^3)*CS,    -(A*E/lee)*S^2,      (6*E*I/lee^2)*C;
          -(6*E*I/lee^2)*S,      (6*E*I/lee^2)*C,     (4*E*I/lee),      (6*E*I/lee^2)*S,      -(6*E*I/lee^2)*C,    (2*E*I/lee);
          -(12*E*I/lee^3)*S^2,   12*E*I/lee^3*CS,     (6*E*I/lee^2)*S,  (12*E*I/lee^3)*S^2,   -(12*E*I/lee^3)*CS,  (6*E*I/lee^2)*S;
          (12*E*I/lee^3)*CS,     -(A*E/lee)*S^2,      -(6*E*I/lee^2)*C, -(12*E*I/lee^3)*CS,   (A*E/lee)*S^2,       -(6*E*I/lee^2)*C;
          -(6*E*I/lee^2)*S,      (6*E*I/lee^2)*C,     (2*E*I/lee),      (6*E*I/lee^2)*S,      -(6*E*I/lee^2)*C,    (4*E*I/lee)];   
end
E = 200e+09; I = 450e-09; A = 0.012;
% Member 1 
x1 = [0, 0, 4.05, 0];
k1 = elestiff_key(E, A, I, x1);
% Member 2 
x2 = [4.05, 0, 4.05, 2.9];
k2 = elestiff_kee(E, A, I, x2);
% Assembly
K = zeros(9,9);
F = zeros(9,1);
K(1:6,1:6) = k1(1:6, 1:6);
K(4:9,4:9) =  K(4:9,4:9) + k2(1:6, 1:6);
F(2) = -101.762e+03;
F(3) = 75.355e+03;
F(5) = -101.76e+03;
F(6) = -75.355e+03;
% making matrix small (3*3) to calculate displacements
Ksmal = K(4:6,4:6);
Fsmall = F(4:6);
% Equations to get the values of displacements
usmall = inv(Ksmal)* Fsmall;
uno = zeros(9,1);
uno([4:6]) = usmall;
% Equations to get the values of support reactions
Ksmall = K([1,2,3,7,8,9],[4:6]);
Frf = Ksmall*usmall;

clear global; clc;
% Introducing a function "elestiff_key" to calculate stiffness matrix (4*4)
function key = elestiff_key(E, A, x)
    lxx = x(3) - x(1);
    lyy = x(4) - x(2);
    lee = sqrt(lxx^2 + lyy^2);
    C = lxx / lee; 
    S = lyy / lee; 
    CS = C * S;
    fa = A * E / lee;
    key = [C^2, CS, -C^2, -CS;
           CS, S^2, -CS, -S^2;
           -C^2, -CS, C^2, CS;
           -CS, -S^2, CS, S^2] * fa;
end
% Introducing a function "elestiff_Q" to calculate stiffness matrix (1*4)
function Q = elestiff_Q(E, A, x)
    lxx = x(3) - x(1);
    lyy = x(4) - x(2);
    lee = sqrt(lxx^2 + lyy^2);
    C = lxx / lee; 
    S = lyy / lee; 
    fa = A * E / lee;
    Q = [-C, -S, C, S] * fa;
end
% Given values
E = 200e9; A = 0.0015;
% Member 1 
x1 = [0, 0, 2.9, 0];
k1 = elestiff_key(E, A, x1);
% Member 2
x2 = [2.9, 0, 5.8, 0];
k2 = elestiff_key(E, A, x2);
% Member 3
x3 = [5.8, 0, 2.9, 2.9];
k3 = elestiff_key(E, A, x3);
% Member 4 
x4 = [2.9, 0, 2.9, 2.9];
k4 = elestiff_key(E, A, x4);
% Member 5
x5 = [0, 0, 2.9, 2.9];
k5 = elestiff_key(E, A, x5);
% Member 6
x6 = [0, 2.9, 2.9, 2.9];
k6 = elestiff_key(E, A, x6);
% Assembly
K = zeros(10,10);
F = zeros(10,1);

K(1:4,1:4) = k1(1:4, 1:4);
K(3:6,3:6) = K(3:6,3:6) + k2(1:4, 1:4);
K(5:8, 5:8) = K(5:8, 5:8) + k3(1:4, 1:4);
K([3, 4, 7, 8],[3, 4, 7, 8]) = K([3, 4, 7, 8], [3, 4, 7, 8])+ k4(1:4, 1:4);
K([1, 2, 7, 8],[1, 2, 7, 8]) = K([1, 2, 7, 8], [1, 2, 7, 8]) + k5(1:4, 1:4);
K(7:10, 7:10) = K(7:10, 7:10)+ k6(1:4, 1:4);
% Given Values
F(5) = 7900;
F(6) = -36000;
%using boundary condition and making matrix small (6*6) to calculate displacements
Ksmall = K([3:8],[3:8]);
Fsmall = F([3:8]);
% Equations for getting the values of all unknown displacements
usmall = inv(Ksmall)* Fsmall;
uno = zeros(10,1);
uno([3:8]) = usmall;
% Equations for getting the values of support reactions
Frf = K([1,2,9,10],[3:8])*usmall;
E = 200e9; A = 0.0015;
% Member 1
x1 = [0, 0, 2.9, 0];
Q1 = elestiff_Q(E, A, x1);
% Member 2
x2 = [2.9, 0, 5.8, 0];
Q2 = elestiff_Q(E, A, x2);
% Member 3
x3 = [5.8, 0, 2.9, 2.9];
Q3 = elestiff_Q(E, A, x3);
% Member 4
x4 = [2.9, 0, 2.9, 2.9];
Q4 = elestiff_Q(E, A, x4);
% Member 5
x5 = [0, 0, 2.9, 2.9];
Q5 = elestiff_Q(E, A, x5);
% Member 6
x6 = [0, 2.9, 2.9, 2.9];
Q6 = elestiff_Q(E, A, x6);
% caculating force value for member 1 
q1 = Q1*uno([1:4]);
% caculating force value for member 2 
q2 = Q2*uno([3:6]);
% caculating force value for member 3 
q3 = Q3*uno([5:8]);
% caculating force value for member 4
q4 = Q4*uno([3, 4, 7, 8]);
% caculating force value for member 5 
q5 = Q5*uno([1, 2, 7, 8]);
% caculating force value for member 6
q6 = Q6*uno([7:10]);
QQQ = [q1;q2;q3;q4;q5;q6]
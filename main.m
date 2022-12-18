clear all

T = 318.15
psat1 = 0.665206 % Presion de saturación acetona
psat2 = 0.291383 % Presion de saturación benceno
XA = [0.04700 0.09630 0.22070 0.29360 0.40110 0.47590 0.61250 0.70450 0.80810 0.90840 0.95290]
XB = 1 - XA
n = length(XA)
P = zeros(1,n); % Memory preallocation
Y1 = zeros(1,n); % Memory preallocation
for i=1:n
    [gamma1, gamma2] = UNIFAC(XA(i), XB(i));
    P(i) = ((psat1 * gamma1 * XA(i)) + (psat2 * gamma2 * XB(i)))*100; % a kilo pascales
    Y1(i) = (psat1 * gamma1 * XA(i))/(P(i)/100); % se divide por 100 para no afectar al vector XB
end

figure
title("Grafica Experimental ELV Acetona(1)/Benceno(2) T=318.15K")
yyaxis left
xlabel('x,y')
ylabel('Presión(KPa)')
plot(XA,P)
hold
yyaxis left
plot(Y1,P)
ylabel('Presión(KPa)')

XA = [0.04700 0.09630 0.22070 0.29360 0.40110 0.47590 0.61250 0.70450 0.80810 0.90840 0.95290]
XB = [0.14440 0.25740 0.44170 0.52040 0.61390 0.66970 0.76140 0.82010 0.88050 0.94180 0.96990]
y = [33.428,36.666,43.230,46.450,50.647,53.293,57.722,60.527,63.380,66.037,67.189]
figure
title("Grafica de Brown I. Smith")
yyaxis left
xlabel('x,y')
ylabel('Presión(KPa)')
plot(XA,y)
hold
yyaxis left
plot(XB,y)
ylabel('Presión(KPa)')

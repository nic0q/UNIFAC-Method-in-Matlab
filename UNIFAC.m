function [gamma1,gamma2] = UNIFAC(xAa, xBb)
% Información del sistema
T=318.15;                   % Temperatura (K)
xA=xAa;                     % Fracción molar del componente A
xB=xBb;                    % Fracción molar del componente B
xi=[xA xB];                 % Registro de matriz de concentración de componentes
NC=2;                       % Número de componentes
NG=3;                       % N° de grupos funcionales necesarios para describir el sistema
SubG=[1 10 19];             % Matriz de subgrupos (1=CH3; 10=ACH; 19=CH3CO), ver tabla H.1 apéndice texto de Smith, Van Ness, Abbott
%CH3OHCH3
%C6H6
% Parámetros de sustancias puras (leidos de tabla H.1, Apéndice H de Smith, Van Ness, Abbott)
Rk=[0.9011 0.5313 1.67];  % Parámetros Rk de los subgrupos funcionales
Qk=[0.848 0.4 1.488];     % Parámetros Qk de los subgrupos funcionales
nA=[1 0 1];                % Frecuencia de los subgrupos en el compuesto A
nB=[0 6 0];                % Frecuencia de los subgrupos en el compuesto B

% Parámetros de mezcla e interacción (leidos de tabla H.2, Apéndice H de Smith, Van Ness, Abbott)
amk=[0 61.13 476.40;-11.12 0 25.77;26.76 140.1 0]; % Matriz de parámetros amk (K)

% Evaluación de los parámetros de mezcla
rA=sum(nA.*Rk);            % Evaluación del parámetro de tamaño, especie A
rB=sum(nB.*Rk);            % Evaluación del parámetro de tamaño, especie B
qA=sum(nA.*Qk);            % Evaluación del parámetro de superficie, especie B
qB=sum(nB.*Qk);            % Evaluación del parámetro de superficie, especie B

ri=[rA rB];                % Registro de matriz de parámetros de tamaño
qi=[qA qB];                % Registro de matriz de parámetros de superficie

ekA=nA.*Qk/qA;             % Evaluación de los parámetros de superficie de subgrupos, especie A
ekB=nB.*Qk/qB;             % Evaluación de los parámetros de superficie de subgrupos, especie B

eki=[ekA; ekB];             % Registro de matriz de parámetros superficiales de subgrupo en ambos compuestos   

Taumk=exp(-amk/T);         % Evaluación de la matriz Taumk

% Evaluación de los parámetros beta por especie y subgrupo funcional
BetaA1=sum(ekA*Taumk(:,1));  % Especie A, subgrupo 1
BetaA2=sum(ekA*Taumk(:,2));  % Especie A, subgrupo 2
BetaA33=sum(ekA*Taumk(:,3)); % Especie A, subgrupo 33
BetaB1=sum(ekB*Taumk(:,1));  % Especie B, subgrupo 1
BetaB2=sum(ekB*Taumk(:,2));  % Especie B, subgrupo 2
BetaB33=sum(ekB*Taumk(:,3)); % Especie B, subgrupo 33

Betaik=[BetaA1 BetaA2 BetaA33;BetaB1 BetaB2 BetaB33]; % Registro como matriz

% Evaluación de los parámetros Theta de cada subgrupo funcional
for i=1:NG
    Thetak(i)=(xA*qA*ekA(i)+xB*qB*ekB(i))/(xA*qA+xB*qB);
end

% Estimación del parámetro sk 
for j=1:NG
    sk(j)=sum(Thetak*Taumk(:,j));
end

% Estimación de los parámetros Ji
Ji=ri./(sum(xi.*ri));

% Estimación de los parámetros Li
Li=qi./(sum(xi.*qi));

% Estimación del Logaritmo de gammaC
for z=1:NC
    LngammaCi(z)=1-Ji(z)+log(Ji(z))-5*qi(z)*(1-Ji(z)/Li(z)+log(Ji(z)/Li(z)));
end

% Estimación del logaritmo de gammaR
for w=1:NC
    LngammaRi(w)=qi(w)*(1-sum(Thetak.*Betaik(w,:)./sk-eki(w,:).*log(Betaik(w,:)./sk)));
end 

% Estimación de Ln(gamma) de cada componente
Lngamma=LngammaCi+LngammaRi;

% Estimación del coeficiente de actividad (gamma) para cada componente
gamma=exp(Lngamma);

% Imprime valores de coefiencientes de actividad de componentes
gamma1 = gamma(1)
gamma2 = gamma(2)
end

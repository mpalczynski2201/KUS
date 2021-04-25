%% KUS - Projektowanie kaskadowej struktury regulacji napędem prądu stałego
% autor: Przemysław Cencora, Mateusz Pałczyński 
% AGH - EAIiIB - III rok, 
% Data modyfikacji: 24.04.2021

%% Clean
%
clear all; close all; clc

%% 1. Dane:
% Dane silnika dla nr.11:
%
Pn = 27*1e3; % moc znamionowa [W]
Un = 230;   % znamionowe napięcie zasilania uzwojenia twornika [V]
In = 141;   % znamionowy prąd twornika [A]
nn = 725;   % prędkość obrotowa [obr/min]
Rt = 0.239; % Rezystancja twornika [ohm]
Lt = 1.63*1e-3; % Indukcyjność twornika [H]
Js = 0.97;  % moment bezwładności silnika [kgm^2]


% Obliczenia na podsrawie danych silnika:
%
JMR = 10*Js; % moment bezwładności maszyny roboczej [kgm^2]
J = Js %+ JMR; % całkowity moment bezwładności [kgm^2]
T = Lt/Rt; % elektormagnetyczna stała czasowa [s]
wn = nn*2*pi/60; % prędkość kątowa [rad/s]
psien = (Un-In*Rt)/wn; % znamionowy strumień skojarzony rotacyjnie z uzw. twornika [Wb]
B = J*Rt/(psien^2); % elektromechaniczna stała czasowa silnika [s]
Mn = In * psien; % znamionowy moment obciążenia [Nm]

% Dane projektowe:
%
lambdan = 2; % ograniczenie prądowe
p = 50; % dopuszczalna krotność prądu
Y = 10 / (2.5*In); % wsp. wzmocnienia toru sprzężenia zwrotnego od prądu twornika
Kt = 10 / (1.2*wn); % wsp. wzmocnienia toru sprzężenia zwrotnego od prędkości kątowej
Kp = 1.5*Un / 10; % Wzmocnienie regulatora tyrystorowego
% Tc = 0.05*Tp; % krok całkowania



%% Ad1
% Obliczenie transmitancji (8):

% GwU
LwU = [0 0 1/psien];
MwU = [B*T B 1];
sys = tf(LwU,MwU)

%GwM
LwM = Rt/psien^2 * [0 T 1];
MwM = [B*T B 1];
sys = tf(LwM,MwM)

%GIU
LIU = [0 B/Rt 0];
MIU = [B*T B 1];
sys = tf(LIU,MIU)

%GIM
LIM = [0 0 1/psien];
MIM = [B*T B 1];
sys = tf(LIM,MIM)


%% Ad2
% Regulator prądu PI - kryterium modułowe:
%
tau0 = 3.3*1e-3; % stała czasowa regulatora tyrystorowego [s]
kz = 1/Y; % zastępczy współczynnik wzmocnienia
beta = 2*tau0; % stała czasowa przebiegu prądu twornika [s]

Kri = T*Rt / (2*Kp*Y*tau0); % wzmocnienie regulatora prądu 
Tri = T; % stała czasowa regulatora prądu [s];
uz0 = lambdan*In/kz;


%Regulator prędkości PI - kryterium symetryczne:
%
Krw = J*Y / (4*Kt*tau0*psien); % wzmocnienie regulatora prędkości
Trw = 8*tau0; % stała czasowa regulatora prędkości [s]

% %Regulator prędkości P:
% %
% domega = 0.02*wn;
% Krw = Mn / (psien*kz*Kt*domega); % wzmocnienie regulatora prędkości




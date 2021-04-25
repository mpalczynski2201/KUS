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



%% Ad1 - Wyznaczyć transmitancje Gwu,Gwm,Giu,Gim i narysować odpowiedzi skokowe prądu twornika I, jego pochodnej i prędkości kątowej omega
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


%% Ad2 - Wyznaczyć nastawy regulatorów i wykonać dla nich symulacje.
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

%% Ad 3 - Wyznaczyć transmitancje układu otwartego. Wykreślić ploty Nyquista i Bodego i na ich podstawie wyznaczyć zapas modułu i fazy oraz określić dopuszczalne opóźnienie
%
%GP - przekształtnik tyrystorowy
tau0=3.3e-3
Lp=[0 Kp]
Mp=[tau0 1]
Gp=tf(Lp,Mp)
%GRi - regulator prądu
Lgri=Kri*[Tri+1]   
Mgri=[Tri]
GRi=tf(Lgri,Mgri)
%GRw - regulator prędkości typu PI
Tr=4*beta
Kw=J/(2*Kt*kz*beta*psien)
Lgrw=Kw*[Tr 1]
Mgrw=[Tr 0]
GRw=tf(Lgrw, Mgrw)
%Gs 
sigma=sum(T)
Mgs=[sigma 1]
Gs=tf(1,Mgs)

%Gj
Gj=tf(1,[J 0])

%sklejamy transmitancje ze sobą
Gskladowe=(GRi*Gp*Gs*psien*Gj)/(1+GRi*Gp*Gs*psien*Gj*Y)
%i dodaje trzy obok siebie
G1=series(GRw,Gskladowe)
G=series(G1,Kt)

margin(G)
nyquist(G)
opoznienie=allmargin(G)

%% Ad 4 - Dokonać dyskretyzacji regulatorów działania ciągłego i wykonać symulacje dla różnych Tp
%Model DyskretnyUklad.slx
Kr=(J*Y)/(4*Kt*tau0*psien)
Ti=T
%Tp=0.0006 % czas próbkowania
%Tp= 0.0009
Tp=0.000001
K1=Kr
K2=Kr*((Tp/Ti)-1)


%% Ad 5 - Wykorzystując kwantyzatory. Doprowadzić do powstania cyklu granicznego
%W modelu DyskretnyUkladKwant dodałem kwantyzatory
UZak = 10 %zakres pomiarowy
Nkwant= 12 %ilość bitów kwantyztora
qkwant = UZak/(2^(Nkwant)-1) %poziom kwantyzacji


%% Ad 6 - Dokonać symulacji rozruchu napędem z momentem obciążenia. Na rysunku z wynikiem nanieść linią prostą wartość dopuszczalną prądu twornika.
V=(Kp*Y/Rt)*((B*beta)/(sqrt(B*T)-beta))
Mu=Mn
dI=(psien*V*Mu)/(((psien.^2)*V)+(J*Kp*Y)) %Wartośc dopuszczalna prądu twornika według wzoru


%% Ad 7 - Dokonać obciążenia udarowego momentem nominalnym Mn napędu pracującego z prędkością omegan

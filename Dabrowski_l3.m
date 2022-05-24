%%% Jakub Dąbrowski wt TP 15:15
clear all;
close all;
% 1) Parametry wejściowe
m=12; % rozmiar populacji m -> liczba parzysta
l=12; % długość łańcucha binarnego chromosomu l
a=0; % przedział [0,2]
b=2; % przedział [0,2]

% 1) Dekodowanie reprezentacji binarnej na rzeczywistą
xb = population(m, l);
C = num2str(xb);
xd = bin2dec(C); % binary into decimal

l = size(xb,2); % przedział [0 2]

xr = (a + xd*(b-a))/(2^l-1);

% 2) Funckja przystosowania z R
f = xr+sin(3*cos(5*xr))+0.8;

% 3) Druga metoda selekcji metoda odwrotnej dystrybuanty
Pr = f/sum(f); % utworzenie wektora prawdopodobieństw
prs = cumsum(Pr); % utworzenie wektora sum cząstkowych
for i = 1:m 
    x(i,1) = rand(1); % losowanie dwóch z przedziału [1 m]
    A = find(prs<=x(i,1));
    Z(i,1) = A(end, 1); % do ktorego przedzialu wpadnie
end

% 4) Wizualizacja ewolucji

%%
% Populacja początkowa 
function P = population(m, l)
    P = randi([0,1],m,l);
end

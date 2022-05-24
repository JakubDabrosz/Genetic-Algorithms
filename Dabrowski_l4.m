%%% Jakub Dąbrowski wt TP 15:15
clear all;
close all;
% 1) Parametry wejściowe
m=40; % rozmiar populacji m -> liczba parzysta
l=20; % długość łańcucha binarnego chromosomu l
lg=15; % liczba generacji lg, przebiegów pętli głównej
n=2; % ilość podłańcuchów
li=1:n; % wektor podziału podłańcuchów
pc=0.25; % prawdopodobieństwo krzyżowania pc
pm=0.02; % prawdopodobieństwo mutacji pm
t=3; % liczba rodziców
i=0; % początek pętli

%%% PROCEDURA GA
% P = population(m, l); % m wierszy, l kolumn INICJALIZACJA POPULACJI 
% save('P', "P")
load("P.mat")

f = evaluate(P, n); % P populacja, n ilość podłańcuchów FUNCKJA PRZYSTOSOWANIA
while i < lg
    i = i + 1;
    P1 = parent_selection(P, f, m, t); % SELEKCJA RODZICÓW
    P2 = recombine(P1, m, l, pc);      % KRZYŻOWANIE
    P3 = mutation(P2, m, l, pm);       % MUTACJA
    P = P3; 
    f = evaluate(P,n);                 % FUNCKJA PRZYSTOSOWANIA NOWEJ POPULACJI
    f_max(i, 1) = max(f);              % WEKTOR F-MAX
    f_sr(i, 1) = mean(f);              % WEKTOR F-ŚR
end

% save('f_max_m10', "f_max")
% save('f_max_m20', "f_max")
% save('f_max_m40', "f_max")
% subplot(1,2,1)
% subplot(1,2,2)


subplot(3,1,1)
load('f_max_m10')
plot(1:1:lg, f_max);
grid on;
xticks([1:1:lg])
xlabel('Generacja lg')
title('fmax w kolejnych generacjach, m = 10') 

subplot(3,1,2)
load('f_max_m20')
plot(1:1:lg, f_max);
grid on;
xticks([1:1:lg])
xlabel('Generacja lg')
title('fmax w kolejnych generacjach, m = 20') 

subplot(3,1,3)
load('f_max_m40')
plot(1:1:lg, f_max);
grid on;
xticks([1:1:lg])
xlabel('Generacja lg')
title('fmax w kolejnych generacjach, m = 40') 


% plot(1:1:lg, f_sr)
% grid on;
% xticks([1:1:lg])
% xlabel('Generacja lg')
% title('fśr w kolejnych generacjach')





% sgtitle('wielkość populacji m = 40') 


%% LAB4
% 1) i 2) Podział podłańcuchów i wyliczenie funkcja przystosowania z wykorzystaniem funckji bin na dec (lab4)
function f = evaluate(P, n)
    m = size(P,1);
    l = size(P,2);
    for i=1:m
        X{i,1} = P(i,1:10); % pierwszy podział (różna lub taka sama długość łańcuchów)
        X{i,2} = P(i,11:20);
        x1 = evaluate_bd(X{i,1});
        x2 = evaluate_bd(X{i,2});
        f(i,1) = (25-(x1-5)^2)*cos(2*x1-5) + (25-(x2-5)^2)*cos(2*(x2-5))+50; % dwuwymiarowa funkcja przystosowania
    end
end




%% LAB2
% *% 2) Populacja początkowa* 

function P = population(m, l)
    P = randi([0,1],m,l);
end

%%
% *% 3) Wyliczenie funkcji przystosowania dla populacji* 

function f = evaluate_sum(P)
    f = sum(P, 2);
end

%%
% *% 4) Selekcja m rodziców turniejowa binarna (t=2) * 

% najpierw losujemy 2 z każdego wiersza o m elementach
function P1 = parent_selection(P, f, m, t)
    for i = 1:m 
        x = randperm(m,t); % losowanie dwóch z przedziału [1 m]
        if f(x(1)) > f(x(2))
            P1(i,:) = P(x(1),:);  
        elseif f(x(2)) > f(x(1))
            P1(i,:) = P(x(2),:);  
        else 
            P1(i,:) = P(x(1),:);   
        end
    end
end

%%
% *% 5) Krzyżowanie jednopunktowe z prawdopodobieństwem pc*=0.7* 
    function P2 = recombine(P1, m, l, pc)
        P2 = P1;
        for i = 1:2:m-1
            if rand(1) < pc
                k = randi([1,l-1]); % punkt krzyżowania k z przedziału [1, l-1]
                P2(i,k+1:l) = P1(i+1,k+1:l); % zamiana końcowek od pozycji k+1
                P2(i+1,k+1:l) = P1(i,k+1:l);
            end
        end
    end
    
%%
% *% 5) Mutacja z prawdopodobieństwem pm*=0.12* 
function P3 = mutation(P2, m, l, pm)
    P3 = P2;
    for i = 1:m
        if rand(1) < pm
                los = randi([1,l]); % losowanie bitu w i-wierszu
                P3(i,los) = ~P2(i,los);
        end
    end
end


%% LAB3
% *% 1) Dekodowanie reprezentacji binarnej na rzeczywistą

function xr = evaluate_bd(P)
    a=0; % przedział [0,10]
    b=10; % przedział [0,10]
    C = num2str(P);
    xd = bin2dec(C); % binary into decimal
    l = size(P,2); % przedział [0 10] dla lab4 
    xr = (a + xd*(b-a))/(2^l-1);
end 
%%
% *% 2) Funckja przystosowania z R
% f = xr+sin(3*cos(5*xr))+0.8;

%%
% *% 3) Druga metoda selekcji metoda odwrotnej dystrybuanty
% Pr = f/sum(f); % utworzenie wektora prawdopodobieństw
% prs = cumsum(Pr); % utworzenie wektora sum cząstkowych
% for i = 1:m 
%    x(i,1) = rand(1); % losowanie dwóch z przedziału [1 m]
%    A = find(prs<=x(i,1));
%    Z(i,1) = A(end, 1); % do ktorego przedzialu wpadnie
% end


%%% Jakub Dąbrowski wt TP 15:15
clear all;
close all;
% 1) Parametry wejściowe
m=10;   % rozmiar populacji m -> liczba parzysta
n=6;    % wymiarowość przestrzeni poszukiwań 
lg=15;  % liczba generacji lg, przebiegów pętli głównej
 
pc=0.7; % prawdopodobieństwo krzyżowania pc=0.7-0.9
pm=1 ;  % prawdopodobieństwo mutacji pm
t=3;    % liczba rodziców 
range = [0 10]; % dziedzina funkcji przystosowania
sigma = 1/100;  % odchylenie standardowe


%%% PROCEDURA GA
P = population(m, n, range);  % m wierszy, n kolumn INICJALIZACJA POPULACJI POCZĄTKOWEJ
f = evaluate(P, range);              % P populacja, n ilość podłańcuchów FUNCKJA PRZYSTOSOWANIA

% while i < lg
%     i = i + 1;
%     P1 = parent_selection(P, f, m, t); % SELEKCJA RODZICÓW
%     P2 = recombine(P1, m, n, pc);      % KRZYŻOWANIE
%     P3 = mutation(P2, m, n, pm);       % MUTACJA
%     P = P3; 
%     f = evaluate(P);                   % FUNCKJA PRZYSTOSOWANIA NOWEJ POPULACJI
%     f_max(i, 1) = max(f);              % WEKTOR F-MAX
%     f_sr(i, 1) = mean(f);              % WEKTOR F-ŚR
% end
% 
% 
% figure(1)
% plot(1:1:lg, f_max);
% grid on;
% xticks([1:1:lg])
% xlabel('Generacja lg')
% title('fmax w kolejnych generacjach, m = 10') 

% figure(2)
% plot(1:1:lg, f_sr)
% grid on;
% xticks([1:1:lg])
% xlabel('Generacja lg')
% title('fśr w kolejnych generacjach')

%% LAB5 %%

%% 1) Populacja - m osobników wektory rzeczywistoliczbowe o n współrzędnych
%% 2) Populacja pcozątkowa losowana z zakresu dziedziny funkcji
% przystosowania
function P = population(m, n, range)
    P = randi([range(1) range(2)],m,n);
end

%% 2)** Obliczenie funkcji przystosowania
% najpierw losujemy 2 z każdego wiersza o m elementach

% Dekodowanie z wykorzystaniem sum 
function xr = evaluate_sum(P, range)
    a=range(1);              % przedział [0,10]
    b=range(2);              % przedział [0,10]
    xd = sum(P,2);           % realnum into dec // sum
    n = size(P,2);           % przedział [0 10] dla lab4 
    xr = (a + xd*(b-a))/(2^n-1);
end 

% Funckja przystosowania wykorzystując laboratorium 4 (podział na łańcuchy)
function f = evaluate(P, range)
    m = size(P,1);
    n = size(P,2);
    for i=1:m
        X{i,1} = P(i,1:n/2);    % pierwszy podział 
        X{i,2} = P(i,n/2+1:n);
        x1 = evaluate_sum(X{i,1}, range);
        x2 = evaluate_sum(X{i,2}, range);
        f(i,1) = (25-(x1-5)^2)*cos(2*x1-5) + (25-(x2-5)^2)*cos(2*(x2-5))+50; % dwuwymiarowa funkcja przystosowania
    end
end

%% 3) Selekcja m rodziców turniejowa binarna (t=2) * 

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

%% 4) Rekombinacja arytmetyczna 

function P2 = recombine(P1, m, n, pc)
        P2 = P1;
        for i = 1:m
            for j = 1:n
                if rand(1) < pc
                    u = rand(1);
                    x(i) = u*P(i, j) + (1-u)*P(i+1, j)
                end
            end
        end
end

%% 5) Mutacja z rozkładem normalnym (gaussowska), z parametrem sigma

% u = 2;      % Wartość oczekiwana (o największym prawdopodobieństwie wystąpienia)
% s = 1;      % Odchylenie standardowe (będące miarą rozrzutu generowanych wartości 
%             % wokół wartości oczekiwanej)
%             
% y1 = s*randn(N,M) + u; 



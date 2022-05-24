%%% Jakub Dąbrowski wt TP 15:15
%%% laboratorium 6
%% Parametry wejściowe
clear all;
close all;
% 1) Parametry wejściowe
m=20;   % rozmiar populacji m -> liczba parzysta
n=2;    % wymiarowość przestrzeni poszukiwań 
lg=15;  % liczba generacji lg, przebiegów pętli głównej
pc=0.7; % prawdopodobieństwo krzyżowania pc=0.7-0.9
pm=1 ;  % prawdopodobieństwo mutacji pm
t=3;    % liczba rodziców 
range=[-500 500]; % dziedzina funkcji przystosowania
sigma=1/100 * 10;  % odchylenie standardowe 1/100 zakresu dziedziny funkcji przystosowania
i=0;    % początek pętli

%% PROCEDURA GA
P = population(m, n, range);  % m wierszy, n kolumn INICJALIZACJA POPULACJI POCZĄTKOWEJ
f = evaluate(P);              % P populacja, n FUNCKJA PRZYSTOSOWANIA
r = sum(f)/m;
R = rand*r;                    % rand z przedziału (0,r), daje numer pierwszego wylosowanego osobnika 

for i = 1:m
    fact(i,1) = sum(f(1:i));  % sumuje pojedyncze elementy i tworzy wektor np. dla elementu f(4) = f(1) + f(2) + f(3)
end

for i = 1:m
    A = find(R<=fact)
    A(1);
    P1(i,:) = P(A(1),:); 
    R = R + r;
end




% while i < lg
%     i = i + 1;
%     P1 = parent_selection(P, f, m, t); % SELEKCJA RODZICÓW
%     P2 = recombine(P1, m, n, pc);      % KRZYŻOWANIE
%     P3 = mutation(P2, m, n, sigma);       % MUTACJA
%     P = P3; 
%     f = evaluate(P);                     % FUNCKJA PRZYSTOSOWANIA NOWEJ POPULACJI
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
% 
% figure(2)
% plot(1:1:lg, f_sr)
% grid on;
% xticks([1:1:lg])
% xlabel('Generacja lg')
% title('fśr w kolejnych generacjach')

%% LAB5 %%

% 1) Populacja - m osobników wektory rzeczywistoliczbowe o n współrzędnych
% 2) Populacja pcozątkowa losowana z zakresu dziedziny funkcji przystosowania
function P = population(m, n, range)   % rzeczywistoliczbowe
    a = range(1);                      % -500 - dziedzina [-500 500]
    b = range(2);                      %  500 - dziedzina [-500 5000]
    P = a + (b-a).*rand(m,n);
end

% 2)** Obliczenie funkcji przystosowania

% Funckja przystosowania wykorzystując laboratorium 4
% wstawić dla n = 2, x1 x2
function f = evaluate(P)
    m = size(P,1);
    n = size(P,2);
    for i=1:m
        x1 = P(i,1);
        x2 = P(i,2);
        f(i,1) = 1000-x1*sin(sqrt(abs(x1))) - x2*sin(sqrt(abs(x2))); % dwuwymiarowa funkcja przystosowania
    end
end

% 3) Selekcja m rodziców turniejowa (t=2) * 
% najpierw losujemy 2 (t=2) z każdego wiersza o m elementach
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

% 4) Rekombinacja arytmetyczna 

function P2 = recombine(P1, m, n, pc)
        P2 = P1;
        for i = 1:2:m
            if rand(1) < pc     % prawdopodobieństwo krzyżowania DLA RODZICÓW!!!!
                u = rand(1);    % liczba losowa z rozkładu jednostajnego z przedziału [0 1]
                for j = 1:n
                    P2(i,j) = u*P1(i, j) + (1-u)*P1(i+1, j);  % REKOMBINACJA
                    P2(i+1,j) = u*P1(i+1, j) + (1-u)*P1(i, j); 
                end
            end
        end
end

% 5) Mutacja z rozkładem normalnym (gaussowska), z parametrem sigma

% trzeba sprawdzić czy po mutacji nie wyjdzie po za [0 10], jeśli wyjdzie
% zostawić starego osobnika lub wstawic dla [1 9] (bezpieczny brzeg)
%%
% *% 5) Mutacja z prawdopodobieństwem pm*=0.12* 
function P3 = mutation(P2, m, n, sigma)
    P3 = P2;
    mut = sigma*randn(m, n) + 10; % rozkład normalny (gauss)
    for i = 1:m
        for j = 1:n
            if mut(i, j) < 10         % jeśli wyjdzie zostawiam starego osobnika 
                P3(i, j) = mut(i, j);  % jeśli jest w dziedzinie to biore nowego osobnika
            end
        end
    end
end

% u = 2;      % Wartość oczekiwana (o największym prawdopodobieństwie wystąpienia)
% s = 1/100;      % Odchylenie standardowe (będące miarą rozrzutu generowanych wartości 
%             % wokół wartości oczekiwanej)
            
% y1 = s*randn(N,M) + u; 



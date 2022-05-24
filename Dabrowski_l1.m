%%% Jakub Dąbrowski wt TP 15:15
clear all;
close all;
% 1) Parametry wejściowe
m=26; % rozmiar populacji m -> liczba parzysta
l=26; % długość łańcucha binarnego chromosomu l
lg=15; % liczba generacji lg, przebiegów pętli głównej
pc=0.7; % prawdopodobieństwo krzyżowania pc
pm=0.12; % prawdopodobieństwo mutacji pm
t=2; % liczba rodziców
i=0; % początek pętli


P = population(m, l);
f = evaluate(P, t); 
while i < lg
    i = i + 1;
    P1 = parent_selection(P, f, m, t);
    P2 = recombine(P1, m, l, pc)
    %mutation(P2, m, l, pm);
    P = P2; 
    f = evaluate(P, t); % funkcja przystosowania nowej populacji
    f_max(i, 1) = max(f) 
    f_sr(i, 1) = mean(f)
end

figure(1)
plot(1:1:lg, f_max);
grid on;
xticks([1:1:lg])
xlabel('Generacja lg')
title('fmax w kolejnych generacjach ') 

figure(2)
plot(1:1:lg, f_sr)
grid on;
xticks([1:1:lg])
xlabel('Generacja lg')
title('fśr w kolejnych generacjach ')


%%
% *% 2) Populacja początkowa* 
function P = population(m, l)
    P = randi([0,1],m,l);
end

%%
% *% 3) Wyliczenie funkcji przystosowania dla populacji* 

function f = evaluate(P, t)
    f = sum(P, t);
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
function mutation(P2, m, l, pm)
    for i = 1:m
        if rand(1) < pm
                los = randi([1,l]); % losowanie bitu w i-wierszu
                P2(i,los) = ~P2(i,los);
        end
    end
end
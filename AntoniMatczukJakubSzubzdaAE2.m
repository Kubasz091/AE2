% Algorytmy Ewolucyjne - Projekt 2: Problem Plecakowy

clear;
clc;
close all;

numerAlbumu = 331129;
rng(numerAlbumu);

N = 32;

item_weights = round(0.1 + 0.9 * rand(N, 1), 1);
item_values = round(1 + 99 * rand(N, 1));

items = table(item_weights, item_values, 'VariableNames', {'Waga', 'Wartosc'});

disp('--- Lista wszystkich przedmiotów ---');
disp(items);

generateLatexTable(items, 'c:\Users\Kubas\OneDrive\Pulpit\AE2\items_table.tex');

W_max = 0.30 * sum(items.Waga);
fprintf('\nMaksymalna pojemność plecaka (W): %.2f\n', W_max);

fitnessFunction = @(x) -sum(x' .* items.Wartosc);

Aineq = items.Waga';
bineq = W_max;

lb = zeros(1, N);
ub = ones(1, N);
IntCon = 1:N;

populationSize = 100;
eliteCount = ceil(0.05 * populationSize);
crossoverFraction = 0.8;
mutationRate = 0.04;

maxGenerations = 200;
stallGenerations = 50;
functionTolerance = 1e-7;

options = optimoptions('ga', ...
    'PopulationSize', populationSize, ...
    'EliteCount', eliteCount, ...
    'CrossoverFraction', crossoverFraction, ...
    'MutationFcn', {@mutationuniform, mutationRate}, ...
    'MaxGenerations', maxGenerations, ...
    'MaxStallGenerations', stallGenerations, ...
    'FunctionTolerance', functionTolerance, ...
    'Display', 'off', ...
    'PlotFcn', {}, ...
    'OutputFcn', @knapsackOutputFcn);

assignin('base', 'ga_stats_history', []);
[solution_x, fval, exitflag, output, final_population, final_scores] = ga(fitnessFunction, N, Aineq, bineq, [], [], lb, ub, [], IntCon, options);

fprintf('\n--- Wyniki Algorytmu Genetycznego ---\n');

if exitflag > 0 || exitflag == 0
    fprintf('Wektor binarny rozwiązania (x):\n');
    disp(solution_x');

    selected_items_indices = find(solution_x == 1);
    selected_weights = items.Waga(selected_items_indices);
    selected_values = items.Wartosc(selected_items_indices);

    total_weight_solution = sum(selected_weights);
    total_value_solution = sum(selected_values);

    fprintf('\nSumaryczna waga przedmiotów w plecaku: %.2f (Maksymalna dozwolona: %.2f)\n', total_weight_solution, W_max);
    fprintf('Sumaryczna wartość przedmiotów w plecaku: %d\n', total_value_solution);

    fprintf('\nWybrane przedmioty:\n');
    if ~isempty(selected_items_indices)
        disp(items(selected_items_indices,:));
        % Generowanie tabeli LaTeX dla wybranych przedmiotów
        selected_items_table = items(selected_items_indices,:);
        generateLatexTable(selected_items_table, 'c:\Users\Kubas\OneDrive\Pulpit\AE2\selected_items_table.tex');
    else
        disp('Żaden przedmiot nie został wybrany.');
    end

    fprintf('\n=== ANALIZA ROZWIĄZANIA PROBLEMU PLECAKOWEGO ===\n');
    fprintf('Problem: Maksymalizacja wartości przedmiotów w plecaku przy ograniczeniu wagowym\n');
    fprintf('Funkcja celu: max Σ(pi * xi), gdzie pi - wartość przedmiotu i, xi ∈ {0,1}\n');
    fprintf('Ograniczenie: Σ(wi * xi) ≤ W, gdzie wi - waga przedmiotu i\n\n');

    fprintf('PARAMETRY PROBLEMU:\n');
    fprintf('- Liczba przedmiotów (n): %d\n', N);
    fprintf('- Maksymalna waga plecaka (W): %.2f\n', W_max);
    fprintf('- Suma wag wszystkich przedmiotów: %.2f\n', sum(items.Waga));
    fprintf('- Współczynnik ograniczenia wagowego: 30%% wszystkich przedmiotów\n\n');

    fprintf('ROZWIĄZANIE ZNALEZIONE PRZEZ ALGORYTM GENETYCZNY:\n');
    fprintf('- Liczba wybranych przedmiotów: %d z %d\n', length(selected_items_indices), N);
    fprintf('- Wykorzystanie pojemności plecaka: %.1f%% (%.2f/%.2f)\n', ...
            100*total_weight_solution/W_max, total_weight_solution, W_max);
    fprintf('- Łączna wartość wybranych przedmiotów: %d\n', total_value_solution);

    if ~isempty(selected_items_indices)
        value_to_weight_ratio = selected_values ./ selected_weights;
        fprintf('- Średni stosunek wartość/waga wybranych przedmiotów: %.2f\n', mean(value_to_weight_ratio));
        fprintf('- Najlepszy stosunek wartość/waga w rozwiązaniu: %.2f\n', max(value_to_weight_ratio));
        fprintf('- Najgorszy stosunek wartość/waga w rozwiązaniu: %.2f\n', min(value_to_weight_ratio));
    end

    % Analiza efektywności rozwiązania
    all_ratios = items.Wartosc ./ items.Waga;
    fprintf('\n- Średni stosunek wartość/waga wszystkich przedmiotów: %.2f\n', mean(all_ratios));
    fprintf('- Algorytm wybrał przedmioty o średnim stosunku: %.2f vs %.2f (wszystkie)\n', ...
            mean(value_to_weight_ratio), mean(all_ratios));

    fprintf('\nWEKTOR BINARNY ROZWIĄZANIA:\n');
    fprintf('x = [');
    for i = 1:N
        fprintf('%d', solution_x(i));
        if i < N
            fprintf(' ');
        end
    end
    fprintf(']\n\n');

    fprintf('INTERPRETACJA WEKTORA BINARNEGO:\n');
    fprintf('- 1 oznacza: przedmiot został wybrany do plecaka\n');
    fprintf('- 0 oznacza: przedmiot nie został wybrany\n');
    fprintf('- Pozycja w wektorze odpowiada numerowi przedmiotu\n\n');

    fprintf('\n--- Parametry Algorytmu Genetycznego i Wyniki Uruchomienia ---\n');

    fprintf('USTAWIENIA ALGORYTMU GENETYCZNEGO:\n');
    fprintf('  Rozmiar populacji (PopulationSize): %d\n', populationSize);
    fprintf('  Liczba osobników elitarnych (EliteCount): %d\n', eliteCount);
    fprintf('  Współczynnik krzyżowania (CrossoverFraction): %.2f\n', crossoverFraction);
    fprintf('  Współczynnik mutacji (MutationRate): %.3f\n', mutationRate); % .3f for typical mutation rates like 0.015
    if iscell(options.MutationFcn) && ~isempty(options.MutationFcn) && isa(options.MutationFcn{1}, 'function_handle')
        fprintf('  Funkcja mutacji: @%s\n', func2str(options.MutationFcn{1}));
    elseif isa(options.MutationFcn, 'function_handle')
        fprintf('  Funkcja mutacji: @%s\n', func2str(options.MutationFcn));
    else
        fprintf('  Funkcja mutacji: (nieznana lub domyślna)\n');
    end

    fprintf('\nKRYTERIA ZATRZYMANIA ALGORYTMU:\n');
    fprintf('  Maksymalna liczba pokoleń (MaxGenerations): %d\n', maxGenerations);
    fprintf('  Liczba pokoleń bez poprawy (MaxStallGenerations): %d\n', stallGenerations);
    fprintf('  Tolerancja funkcji celu (FunctionTolerance): %e\n', functionTolerance);

    fprintf('\nSTATYSTYKI WYKONANIA ALGORYTMU:\n');
    fprintf('  Liczba wykonanych pokoleń: %d\n', output.generations);
    fprintf('  Liczba obliczeń funkcji celu: %d\n', output.funccount);
    fprintf('  Komunikat zakończenia GA: %s\n', output.message);

    num_non_elite_offspring = populationSize - eliteCount;
    num_non_elite_offspring = max(0, num_non_elite_offspring); % Ensure non-negative

    num_crossover_children = round(crossoverFraction * num_non_elite_offspring);
    num_mutation_children = num_non_elite_offspring - num_crossover_children;

    fprintf('\nLICZNOŚCI POTOMKÓW W KAŻDEJ GENERACJI:\n');
    fprintf('  Elitarni: %d\n', eliteCount);
    fprintf('  Skrzyżowani: %d\n', num_crossover_children);
    fprintf('  Zmutowani: %d\n', num_mutation_children);

else
    fprintf('Algorytm genetyczny nie znalazł rozwiązania lub został przerwany.\n');
    fprintf('Kod zakończenia (exitflag): %d\n', exitflag);
    fprintf('Komunikat: %s\n', output.message);
end

ga_stats_history_data = evalin('base', 'ga_stats_history');

if ~isempty(ga_stats_history_data)
    generations = ga_stats_history_data(:,1);
    max_actual_values_pop = -ga_stats_history_data(:,2);
    mean_actual_values_pop = -ga_stats_history_data(:,3);
    min_actual_values_pop = -ga_stats_history_data(:,4);
    variance_actual_values_pop = ga_stats_history_data(:,5);

    fig1 = figure('Name', 'Statystyki Wartości Funkcji Celu', 'Position', [100, 100, 800, 600]);
    hold on;
    plot(generations, max_actual_values_pop, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Max Wartość w Populacji');
    plot(generations, mean_actual_values_pop, 'g-', 'LineWidth', 1.5, 'DisplayName', 'Średnia Wartość w Populacji');
    plot(generations, min_actual_values_pop, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Min Wartość w Populacji');
    hold off;
    xlabel('Pokolenie');
    ylabel('Wartość Plecaka');
    title('Statystyki Wartości Plecaka w Populacji vs Pokolenie');
    legend show;
    grid on;

    saveas(fig1, 'c:\Users\Kubas\OneDrive\Pulpit\AE2\statistics_plot.png');
    fprintf('Wykres statystyk zapisany jako: statistics_plot.png\n');

    fig2 = figure('Name', 'Wariancja Wartości Funkcji Celu', 'Position', [200, 200, 800, 600]);
    plot(generations, variance_actual_values_pop, 'm-', 'LineWidth', 1.5, 'DisplayName', 'Wariancja Wartości');
    xlabel('Pokolenie');
    ylabel('Wariancja Wartości Plecaka');
    title('Wariancja Wartości Plecaka w Populacji vs Pokolenie');
    legend show;
    grid on;

    saveas(fig2, 'c:\Users\Kubas\OneDrive\Pulpit\AE2\variance_plot.png');
    fprintf('Wykres wariancji zapisany jako: variance_plot.png\n');

else
    disp('Nie zebrano danych do wykresów statystyk.');
end

function [state, options, optchanged] = knapsackOutputFcn(options, state, flag)
    optchanged = false;
    persistent history_data;

    if strcmp(flag, 'init')
        fprintf('Inicjalizacja algorytmu genetycznego...\n');
        if isfield(options, 'MaxGenerations') && ~isempty(options.MaxGenerations)
             max_gen_val = options.MaxGenerations;
        else
            max_gen_val = 200;
        end
        history_data = NaN(max_gen_val, 5);
    end

    current_generation = state.Generation;

    if isempty(current_generation) || current_generation == 0
        return;
    end

    scores = state.Score;

    if mod(current_generation, 10) == 0
        best_actual_value = -min(scores);
        mean_actual_value = -mean(scores);
        fprintf('Pokolenie %d: Najlepsza wartość = %d, Średnia = %.1f\n', ...
                current_generation, best_actual_value, mean_actual_value);
    end

    if current_generation > size(history_data, 1)
        extension_size = current_generation - size(history_data,1) + 100;
        history_data = [history_data; NaN(extension_size, 5)];
    end

    history_data(current_generation, 1) = current_generation;
    history_data(current_generation, 2) = min(scores);
    history_data(current_generation, 3) = mean(scores);
    history_data(current_generation, 4) = max(scores);
    history_data(current_generation, 5) = var(scores);

    if strcmp(flag, 'done')
        fprintf('Algorytm genetyczny zakończony.\n');
        last_gen_idx = find(~isnan(history_data(:,1)), 1, 'last');
        if isempty(last_gen_idx)
            assignin('base', 'ga_stats_history', []);
        else
            assignin('base', 'ga_stats_history', history_data(1:last_gen_idx,:));
        end
        clear history_data;
    end
end

function generateLatexTable(itemsTable, filename)
    % Funkcja generująca tabelę LaTeX dla przedmiotów
    fid = fopen(filename, 'w');
    if fid == -1
        warning('Nie można utworzyć pliku: %s', filename);
        return;
    end

    fprintf(fid, '\\begin{table}[h!]\n');
    fprintf(fid, '\\centering\n');
    fprintf(fid, '\\begin{tabular}{|c|c|c|c|}\n');
    fprintf(fid, '\\hline\n');
    fprintf(fid, 'Nr & Waga & Wartość & Stosunek V/W \\\\\n');
    fprintf(fid, '\\hline\n');

    if contains(filename, 'selected')
        original_indices = evalin('caller', 'selected_items_indices');
        for i = 1:height(itemsTable)
            ratio = itemsTable.Wartosc(i) / itemsTable.Waga(i);
            fprintf(fid, '%d & %.1f & %d & %.2f \\\\\n', ...
                    original_indices(i), itemsTable.Waga(i), itemsTable.Wartosc(i), ratio);
        end
    else
        for i = 1:height(itemsTable)
            ratio = itemsTable.Wartosc(i) / itemsTable.Waga(i);
            fprintf(fid, '%d & %.1f & %d & %.2f \\\\\n', ...
                    i, itemsTable.Waga(i), itemsTable.Wartosc(i), ratio);
        end
    end

    fprintf(fid, '\\hline\n');
    fprintf(fid, '\\end{tabular}\n');
    if contains(filename, 'selected')
        fprintf(fid, '\\caption{Wybrane przedmioty w rozwiązaniu problemu plecakowego}\n');
        fprintf(fid, '\\label{tab:selected_items}\n');
    else
        fprintf(fid, '\\caption{Lista wszystkich przedmiotów w problemie plecakowym}\n');
        fprintf(fid, '\\label{tab:all_items}\n');
    end
    fprintf(fid, '\\end{table}\n');

    fclose(fid);

    if contains(filename, 'selected')
        fprintf('Tabela LaTeX wybranych przedmiotów zapisana do: %s\n', filename);
    else
        fprintf('Tabela LaTeX wszystkich przedmiotów zapisana do: %s\n', filename);
    end
end

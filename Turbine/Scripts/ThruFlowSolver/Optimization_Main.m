% Main Optimization Script:
clear;clc;close all;
% addpath ../
% ddpath ../FlowTools/OneDim_FlowTools

% Steps:
% RNG a generation.
% Check all of their fitness
% Do sex
% Repeat generation

% Controls:
evolution = 9; %Just a parameter for keeping track of different evolutions (plotting)
population_size = 100; %Number of turbines in a generation.
parent_count = 2;
mutation_potency = 0.20; %How strong is a genetic mutation
pwr_weight    = 10;
tt_eff_weight = 0.0;
ts_eff_weight = 0.2;
stress_weight = 0.0;

weights = [pwr_weight, tt_eff_weight, ts_eff_weight, stress_weight];

%Limits:
lim_reac        = [0, 1];
lim_height      = [1, 100];
lim_backpress   = [1, 3].*100.*1000;
lim_rotorH2C    = [0.1, 0.9];
lim_rotorT2C    = [0.1, 0.9];
lim_statorH2C   = [0.1, 0.9];
lim_statorT2C   = [0.1, 0.9];

% if length(dir(fullfile('Optimization/', 'gen*'))) == 0
%     %% Create the first generation
%     gen1 = struct('inp', {}, ...
%         'res',{}, ...
%         'pwr',{}, ...
%         'tte',{}, ...
%         'tse',{}, ...
%         'strs', {}, ...
%         'fail', {}, ...
%         'fit', {});
%     for i=1:population_size
%         random_inp = [...
%             rand_interval(lim_reac), ...
%             rand_interval(lim_height), ...
%             rand_interval(lim_backpress), ...
%             rand_interval(lim_rotorH2C), ...
%             rand_interval(lim_rotorT2C), ...
%             rand_interval(lim_statorH2C), ...
%             rand_interval(lim_statorT2C)];
% 
%         res = Calc_Stage_Perf(random_inp(1), random_inp(2), random_inp(3), random_inp(4), random_inp(5), random_inp(6), random_inp(7)); 
% 
%         gen1(i).inp  = random_inp;
%         gen1(i).res  = res;
%         gen1(i).pwr  = res.power;
%         gen1(i).tte  = res.tt_eff;
%         gen1(i).tse  = res.ts_eff;
%         gen1(i).strs = res.stress_frac;
%         gen1(i).fail = res.fail;
%     end
% 
%     % Calc fitnesses
%     for i=1:population_size
%         indiv = gen1(i);
%         if indiv.fail == "success"
%             gen1(i).fit = calc_fitness(gen1, indiv.pwr, indiv.tte, indiv.tse, indiv.strs, weights);
%         else
%             gen1(i).fit = 0;
%         end
%     end
% 
%     save('Optimization\gen1.mat', "gen1");
% end


% Do turbine sex:
if ~isempty(dir(fullfile('Optimization', strcat('ev',char(string(evolution)),'gen*'))))
    [last_gen_num, last_gen] = get_latest_gen(evolution);
else
    last_gen_num = 0;
end

new_gen = struct('inp', {}, ...
    'res',{}, ...
    'pwr',{}, ...
    'tte',{}, ...
    'tse',{}, ...
    'strs', {}, ...
    'fail', {}, ...
    'fit', {}, ...
    'ind_name', {});
for i=1:population_size
    random_inp = [...
        rand_interval(lim_reac), ...
        rand_interval(lim_height), ...
        rand_interval(lim_backpress), ...
        rand_interval(lim_rotorH2C), ...
        rand_interval(lim_rotorT2C), ...
        rand_interval(lim_statorH2C), ...
        rand_interval(lim_statorT2C)];
   
    % Do mutations:
    if (last_gen_num == 0)
        new_inp = random_inp;
    else
        parents = pick_parents(parent_count, last_gen);

        new_inp = mean(vertcat(parents.genes), 1);

        % Do mutation:
        new_inp = new_inp + random_inp.*mutation_potency;
    end

    res = Calc_Stage_Perf(new_inp(1), new_inp(2), new_inp(3), new_inp(4), new_inp(5), ...
        new_inp(6), new_inp(7));

    new_gen(i).inp  = new_inp;
    new_gen(i).res  = res;
    new_gen(i).pwr  = res.power;
    new_gen(i).tte  = res.tt_eff;
    new_gen(i).tse  = res.ts_eff;
    new_gen(i).strs = res.stress_frac;
    new_gen(i).fail = res.fail;
    new_gen(i).ind_name = "turbine_g"+(last_gen_num+1)+"n"+i;
end

% Calc fitnesses
for i=1:population_size
    indiv = new_gen(i);
    if indiv.fail == "success"
        new_gen(i).fit = calc_fitness(new_gen, indiv.pwr, indiv.tte, indiv.tse, indiv.strs, weights);
    else
        new_gen(i).fit = 0;
    end
end

gen_name = "ev" + char(string(evolution)) + "gen" + (last_gen_num+1) + ".mat";
save('Optimization\'+gen_name, "new_gen");


function num = rand_interval(interval)
    num = (interval(2)-interval(1)).*rand()+interval(1);
end

function fitness = calc_fitness(generation, pwr, tte, tse, strs, weights)
    % max_pwr = max([generation.pwr].*([generation.fail] == "success")); 
    targ_pwr = 80.*1000;

    pwr_frac = abs(targ_pwr - pwr)./targ_pwr;

    pwr_frac = 1 - pwr_frac;
    
    norm_weights = weights./sum(weights);

    fitness = sum([pwr_frac, tte, tse, (1-strs)].*norm_weights);    
end

function [gen_num, gen_data] = get_latest_gen(evolution)
    gen_files = dir(fullfile('Optimization', strcat('ev',char(string(evolution)),'gen*')));
     
    gen_nums = zeros(size(gen_files));
    for i = 1:length(gen_files)
        gen_nums(i) = str2double(gen_files(i).name(7));
    end
    
    gen_num = max(gen_nums);
    fname = "Optimization/" + strcat('ev',char(string(evolution)),'gen', char(string(gen_num))) + ".mat";
    gen_data = importdata(fname);
end

function parents = pick_parents(count, gen)
    parents = struct('gen_idx', {}, 'indiv', {}, 'genes', {});
    
    for n=1:count
        weighted_probs = [gen.fit]./sum([gen.fit]);
        cum_weights = weighted_probs;
        for i = 1:length(weighted_probs)
            if (weighted_probs(i) > 0)
                cum_weights(i) = weighted_probs(i) + sum(weighted_probs(1:i-1));
            end
        end
    
        rand_num = rand;
    
        picked_idx = find(rand_num < cum_weights, 1, 'last');
        picked_parent = gen(picked_idx);
    
        parents(n).gen_idx = picked_idx;
        parents(n).indiv = picked_parent;
        parents(n).genes = picked_parent.inp;

        gen(picked_idx).fit = 0; %Remove already selected parents from consideration
    end
end
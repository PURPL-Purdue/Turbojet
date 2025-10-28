clear;clc

type = "stator"
evolution = 1

gen_files = dir(fullfile(type + '_optimization', type + strcat('_ev',char(string(evolution)),'gen*')));

gens = struct('gen', {});
for i = 1:length(gen_files)
    gens{i} = importdata(fullfile(type + '_optimization', gen_files(i).name));
end


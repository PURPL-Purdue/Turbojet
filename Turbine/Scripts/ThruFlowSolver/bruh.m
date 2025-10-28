clear;clc;

evolution = 17;
searchgen = 2;

leeway = 0.2;
target_alpha = 64.48;
target_beta2 = 23.67;
target_beta3 = 51;

gen_files = dir(fullfile('Optimization', strcat('ev',char(string(evolution)),'gen*')));

gens = struct('gen', {});
for i = 1:length(gen_files)
    gens{i} = importdata(fullfile('Optimization', gen_files(i).name));
end

alphas_gen = ones(1,300);
% beta2s_gen = ones(1,300);
% beta3s_gen = ones(1,300);
gen = gens{searchgen};

for i = 1:length(alphas_gen)
    alphas_gen(i) = real(gen(i).res.sol.alpha2_deg);
end
% for i = 1:length(beta2s_gen)
%     beta2s_gen(i) = real(gen(i).res.sol.beta2_deg);
% end
% for i = 1:length(beta3s_gen)
%     beta3s_gen(i) = real(gen(i).res.sol.beta3_deg);
% end


found = false;
for i = 1:length(alphas_gen)
    if abs(alphas_gen(i) - target_alpha) < leeway
        i
    end
end
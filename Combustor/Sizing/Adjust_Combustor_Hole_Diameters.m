%% outer liner dimensions
outerPrimaryDiameter = 6.6953; % outer primary hole Diameter (mm)
outerSecondaryDiameter = 8.415; % outer secondary hole Diameter (mm)
outerDillutionDiameter = 9.829; % outer Dillution hole Diameter (mm)

outerPrimaryCount = 24; % outer primary hole count
outerSecondaryCount = 16; % outer secondary hole count 
outerDillutionCount = 16; % outer Dillution hole count

outerPrimaryLayers = 2; % number of layers of holes in primary
outerSecondaryLayers = 2; % number of layers of holes in secondary
outerDillutionLayers = 2; % number of layers of holes in dillution

outerLinerDiameterAvg = (84.8 + 88.8) / 2; % avg diameter of the outer liner (mm)

%% inner liner dimensions
innerPrimaryDiameter = 4.3274; % inner primary hole Diameter (mm)
innerSecondaryDiameter = 5.374; % inner primary hole Diameter (mm)
innerDillutionDiameter = 6.0811; % inner primary hole Diameter (mm)

innerPrimaryCount = 24; % inner primary hole count
innerSecondaryCount = 16; % inner primary hole count
innerDillutionCount = 16; % inner primary hole count

innerPrimaryLayers = 2; % number of layers of holes in primary
innerSecondaryLayers = 2; % number of layers of holes in secondary
innerDillutionLayers = 2; % number of layers of holes in dillution

innerLinerDiameterAvg = (38.1 + 41.1) / 2; % avg diameter of the inner liner (mm)

innerPrimaryRequestedDiameter = input("Enter inner primary hole diameter (mm): ");
innerSecondaryRequestedDiameter = input("Enter inner secondary hole diameter (mm): ");
innerDillutionRequestedDiameter = input("Enter inner dillution hole diameter (mm): ");

fprintf("\n");

%% calculate areas
outerPrimaryArea = outerPrimaryCount * pi * (outerPrimaryDiameter / 2) ^ 2; % total outer primary hole area (mm^2)
outerSecondaryArea = outerSecondaryCount * pi * (outerSecondaryDiameter / 2) ^ 2; % total outer secondary hole area (mm^2)
outerDillutionArea = outerDillutionCount * pi * (outerDillutionDiameter / 2) ^ 2; % total outer dillution hole area (mm^2)

innerPrimaryArea = innerPrimaryCount * pi * (innerPrimaryDiameter / 2) ^ 2; % total inner primary hole area (mm^2)
innerSecondaryArea = innerSecondaryCount * pi * (innerSecondaryDiameter / 2) ^ 2; % total inner secondary hole area (mm^2)
innerDillutionArea = innerDillutionCount * pi * (innerDillutionDiameter / 2) ^ 2; % total inner dillution hole area (mm^2)

totalPrimaryArea = outerPrimaryArea + innerPrimaryArea; % total primary hole area (mm^2)
totalSecondaryArea = outerSecondaryArea + innerSecondaryArea; % total secondary hole area (mm^2)
totalDillutionArea = outerDillutionArea + innerDillutionArea; % total dillution hole area (mm^2)
 
%% calculated height of holes
outerPrimaryAdjustedDiameter = sqrt((totalPrimaryArea / pi - innerPrimaryCount * (innerPrimaryRequestedDiameter / 2) ^ 2) / outerPrimaryCount) * 2; % adjusted outer primary hole diamter (mm)
outerSecondaryAdjustedDiameter = sqrt((totalSecondaryArea / pi - innerSecondaryCount * (innerSecondaryRequestedDiameter / 2) ^ 2) / outerSecondaryCount) * 2; % adjusted outer secondary hole diamter (mm)
outerDillutionAdjustedDiameter = sqrt((totalDillutionArea / pi - innerDillutionCount * (innerDillutionRequestedDiameter / 2) ^ 2) / outerDillutionCount) * 2; % adjusted outer dillution hole diamter (mm)

%% print values
fprintf("Outer Primary Hole Height (mm): %.4f\n", outerPrimaryAdjustedDiameter);
fprintf("Outer Secondary Hole Height (mm): %.4f\n", outerSecondaryAdjustedDiameter);
fprintf("Outer Dillution Hole Height (mm): %.4f\n", outerDillutionAdjustedDiameter);

fprintf("\n");


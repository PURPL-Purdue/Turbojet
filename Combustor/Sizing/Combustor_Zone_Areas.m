%% outer liner dimensions
outerPrimaryDiameter = 6.6953; % outer primary hole diamter (mm)
outerSecondaryDiameter = 8.415; % outer secondary hole diamter (mm)
outerDilutionDiameter = 9.829; % outer Dilution hole diamter (mm)

outerPrimaryCount = 24; % outer primary hole count
outerSecondaryCount = 16; % outer secondary hole count 
outerDilutionCount = 16; % outer Dilution hole count

outerPrimaryLayers = 2; % number of layers of holes in primary
outerSecondaryLayers = 2; % number of layers of holes in secondary
outerDillutionLayers = 2; % number of layers of holes in dillution

outerLinerDiamterAvg = (98.9 + 94.9) / 2; % avg diameter of the outer liner (mm)

%% inner liner dimensions
innerPrimaryDiameter = 4.3274; % inner primary hole diamter (mm)
innerSecondaryDiameter = 5.374; % inner primary hole diamter (mm)
innerDilutionDiameter = 6.0811; % inner primary hole diamter (mm)

innerPrimaryCount = 24; % inner primary hole count
innerSecondaryCount = 16; % inner primary hole count
innerDilutionCount = 16; % inner primary hole count

innerPrimaryLayers = 2; % number of layers of holes in primary
innerSecondaryLayers = 2; % number of layers of holes in secondary
innerDillutionLayers = 2; % number of layers of holes in dillution

innerLinerDiamterAvg = (42.1 + 38.1) / 2; % avg diameter of the inner liner (mm)

%% calculate areas
outerPrimaryArea = outerPrimaryCount * pi * (outerPrimaryDiameter / 2) ^ 2 ; % outer total primary hole area
outerSecondaryArea = outerSecondaryCount * pi * (outerSecondaryDiameter / 2) ^ 2 ; % outer total secondary hole area
outerDilutionArea = outerDilutionCount * pi * (outerDilutionDiameter / 2) ^ 2 ; % outer total dilution hole area

innerPrimaryArea = innerPrimaryCount * pi * (innerPrimaryDiameter / 2) ^ 2 ; % outer total primary hole area
innerSecondaryArea = innerSecondaryCount * pi * (innerSecondaryDiameter / 2) ^ 2 ; % outer total secondary hole area
innerDilutionArea = innerDilutionCount * pi * (innerDilutionDiameter / 2) ^ 2 ; % outer total dilution hole area

%% calculated height of holes
outerPrimaryHeight = outerPrimaryArea / (pi * outerLinerDiamterAvg) / outerPrimaryLayers; % outer primary zone hole height for axisymmetric model (mm)
outerSecondaryHeight = outerSecondaryArea / (pi * outerLinerDiamterAvg) / outerSecondaryLayers; % outer secondary zone hole height for axisymmetric model (mm)
outerDillutionHeight = outerDilutionArea / (pi * outerLinerDiamterAvg) / outerDillutionLayers; % outer dillution zone hole height for axisymmetric model (mm)

innerPrimaryHeight = innerPrimaryArea / (pi * innerLinerDiamterAvg) / innerPrimaryLayers; % outer primary zone hole height for axisymmetric model (mm)
innerSecondaryHeight = innerSecondaryArea / (pi * innerLinerDiamterAvg) / innerSecondaryLayers; % outer secondary zone hole height for axisymmetric model (mm)
innerDillutionHeight = innerDilutionArea / (pi * innerLinerDiamterAvg) / innerDillutionLayers; % outer dillution zone hole height for axisymmetric model (mm)

%% print values
fprintf("Outer Primary Hole Height (m): %f\n", outerPrimaryHeight / 100);
fprintf("Outer Secondary Hole Height (m): %f\n", outerSecondaryHeight / 100);
fprintf("Outer Dillution Hole Height (m): %f\n", outerDillutionHeight / 100);

fprintf("\n");

fprintf("Inner Primary Hole Height (m): %f\n", innerPrimaryHeight / 100);
fprintf("Inner Secondary Hole Height (m): %f\n", innerSecondaryHeight / 100);
fprintf("Inner Dillution Hole Height (m): %f\n", innerDillutionHeight / 100);
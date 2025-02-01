folder = fileparts(fileparts(mfilename('fullpath')))
a = [2 3 4 5]
writematrix(a, fullfile(folder, 'Curve_Files\asdf.txt'))
function [error, n, method, symbols, equations, initialConditions, max_iter, epsilon] = readfile()
%{
file format:
3                       ->  number of equations (also number of variables)
LU Decomposition        ->  root finding method
a b c                   ->  space separated variable symbols
10*a - 7*b              ->  function 1                           \
-3*a + 2.099*b + 6*c    ->  function 2                            |--> n equations in n lines
5*a - b + 5*c           ->  function 3                           /
0 0 0                   ->  space separated initial conditions   \
100                     ->  maximum number of iterations          |--> only in Gauss Seidel method
1e-05                   ->  minimum allowable relative error     /
%}

% open file
filter = {'.txt'};
[name, path] = uigetfile(filter);
directory = [path name];
if length(directory) == 2
    fid = -1;
else
    fid = fopen(directory);
end
% check if file exists
if (fid < 0)
    error = 'You did not choose a file';
    n = ''
    method = '';
    symbols = '';
    equations = '';
    initialConditions = '';
    max_iter = '';
    epsilon = '';
    return
end
% read number of equations
n = fgetl(fid);
% read method name
method = fgetl(fid);
% read symbols
symbols = fgetl(fid);
% read equations
equations = fgetl(fid);
for i = 2 : n
    equations = equations + "\n" + fgetl(fid);
end
compose(equations);     % convert \n to actual newlines
% read special inputs for Gauss-Seidel Method
if strcmp(method, 'Gauss-Seidel Method')
    initialConditions = fgetl(fid);
    max_iter = fgetl(fid);
    epsilon = fgetl(fid);
end
% close file
fclose(fid);

% validate epsilon
if epsilon == -1
    epsilon = '0.00001';
else
    es = str2num(epsilon);
    if length(es) ~= 1
        error = 'Invalid accuracy format';
        return
    end
    if (es < 0) || (es > 100)
        error = 'Invalid accuracy format';
        return
    end
end
% validate max_iter
if max_iter == -1
    max_iter = '50';
else
    i = str2num(epsilon);
    if length(i) ~= 1
        error = 'Invalid number of maximum iterations';
        return
    end
    if i < 1
        error = 'Invalid number of maximum iterations';
        return
    end
end

% No errors
error = 0;

end
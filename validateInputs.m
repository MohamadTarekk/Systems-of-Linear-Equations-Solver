function [error] = validateInputs(n, method, initialConditions, max_iter, epsilon, isIterative)
%{
n                   ->  number of equations (also number of variables)
method              ->  root finding method
isIterative         ->  if isIterative == true -> validate: initialConditions, max_iter, epsilon
initialConditions   ->  space separated initial conditions
max_iter            ->  maximum number of iterations
epsilon             ->  minimum allowable relative error
%}

    % initialize error as 0, ie no error
    error = 0;
    % validate n
    if n == -1
        error = 'Number of Equations is missing';
        return;
    end
    n = str2double(n);
    if(~isequaln(n, NaN) && length(n)==1)
        if n < 1
            error = 'Invalid format for Number of Equations';
            return;
        end
    end
    % validate method
    if method == -1
        error = 'Method is missing';
        return;
    end
    if ~strcmp(method, 'Gauss Elimination') && ~strcmp(method, 'Gauss-Jordan Elimination') && ~strcmp(method, 'LU Decomposition') && ~strcmp(method, 'Gauss-Seidel Method') && ~strcmp(method, 'All')
        error = 'Unsupported method';
        return;
    end
    % validate iterative method parameters
    if(isIterative)
        % validate initialConditions
        if initialConditions == -1
            error = 'Initial Conditions are missing';
            return;
        end
        [initialConditions, success] = str2num(initialConditions); %#ok<ST2NM>
        if ~success
            error = 'Invalid format for Initial Conditions';
            return;
        end
        if ~(length(initialConditions) == n)
            error = 'Initial Conditions are missing';
            return;
        end
        % validate Maximum Number of Iterations
        if ~(max_iter == -1)
            max_iter = str2double(max_iter);
            if(~isequaln(max_iter, NaN) && length(max_iter)==1)
                if max_iter < 1
                    error = 'Invalid format for Maximum Number of Iterations';
                    return;
                end
            end
        end
        % validate Epsilon
        if ~(epsilon == -1)
            epsilon = str2double(epsilon);
            if(~isequaln(epsilon, NaN) && length(epsilon)==1)
                if epsilon < 0 || epsilon > 100
                    error = 'Invalid format for Epsilon';
                    return;
                end
            end
        end
    end
end
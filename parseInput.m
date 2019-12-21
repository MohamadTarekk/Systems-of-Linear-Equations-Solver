function [error, n, initialConditions, max_iter, epsilon, coeffsMatrix, resultsMatrix, symbols] = parseInput(n, isIterative, initialConditions, max_iter, epsilon, equations, symbols)
    % symbols: symbols fo variables in the input
    % equations: input functions as array of strings
    % n: number of equations (also number of variables)

    % parse n
     n = str2double(n);
    % parse iterative method parameters
    if isIterative
        % parse Initial Conditions
        [initialConditions] = str2num(initialConditions); %#ok<ST2NM>
        % parse Maximum Number of Iterations
        if max_iter == -1
            max_iter = 50;
        else
             max_iter = str2double(max_iter);
        end
        % parse Maximum Number of Iterations
        if epsilon == -1
            epsilon = 0.00001;
        else
             epsilon = str2double(epsilon);
        end
    end
    
    % initialize error to 0, ie no error
    error = 0;
    % split symbols
    symbols = strtrim(symbols);
    symbols = split(symbols, ' ');
    if n ~= length(symbols)
        error = 'Symbols are missing';
        [n, initialConditions, max_iter, epsilon, coeffsMatrix, resultsMatrix, symbols] = setDefaults();
        return;
    end
    % parse equations
    equations = splitlines(equations);
    % initialize matrices
    coeffsMatrix = zeros(n);
    resultsMatrix = zeros(n, 1);
    % fill the matrices
    vars = convertStringsToChars(symbols);
    vars{length(vars)+1} = '1';
    for row = 1 : n
        try
            % create symbolic function
            symbolic_expression = str2sym(string(equations(row)));
            % extract coeffiecients
            [coefficients,terms] = coeffs(symbolic_expression);
        catch
            error = sprintf('Invalid syntax at equation %d', row);
            [n, initialConditions, max_iter, epsilon, coeffsMatrix, resultsMatrix, symbols] = setDefaults();
            return;
        end
        % add coefficients to the matrices
        termsStrings = arrayfun(@char, terms, 'uniform', 0);
        if ~all(ismember(termsStrings, vars))
            error = 'Variables of Equations Field do not match Variables filed';
            [n, initialConditions, max_iter, epsilon, coeffsMatrix, resultsMatrix, symbols] = setDefaults();
            return;
        end
        [coeffsMatrix, resultsMatrix, error] = appendCoeffs(coeffsMatrix, resultsMatrix, symbols, coefficients, terms, row);
        if(error ~= 0)
            [n, initialConditions, max_iter, epsilon, coeffsMatrix, resultsMatrix, symbols] = setDefaults();
            return;
        end
    end    
end

function [coeffsMatrix, resultsMatrix, error] = appendCoeffs(coeffsMatrix, resultsMatrix, symbols, coefficients, terms, row)
    resSymbol = sym('1');
    for column = 1 : length(symbols)
        for j = 1 : length(terms)
            if isequaln(sym(symbols(column)), terms(j))
                coeffsMatrix(row, column) = coefficients(j);
                terms(j) = sym('0');    %mark as used
            elseif isequaln(resSymbol, terms(j))
                resultsMatrix(row) = -1 * coefficients(j);
                terms(j) = sym('0');    %mark as used
            end
        end
    end
    % check if there are unused terms, ie invalid symbols in the equations
    error = 0;
    for column = 1 : length(terms)
        if ~isequaln(sym('0'), terms(j))
            error = 'Unknown variables were found in the equations';
            return;
        end
    end
end

function [n, initialConditions, max_iter, epsilon, coeffsMatrix, resultsMatrix, symbols] = setDefaults()
    n = '';
    initialConditions = '';
    max_iter = '50';
    epsilon = '0.00001';
    coeffsMatrix = '';
    resultsMatrix = '';
    symbols = '';
end
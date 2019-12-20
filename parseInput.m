function [coeffsMatrix, resultsMatrix, symbols, error] = parseInput(symbols, equations, n)
    % symbols: symbols fo variables in the input
    % str: input functions as array of strings
    % n: number of equations (also number of variables)

    % split symbols
    symbols = split(symbols, ' ');
    if n ~= length(symbols)
        error = true;
        coeffsMatrix = [];
        resultsMatrix = [];
        return;
    end
    
    % initialize matrices
    coeffsMatrix = zeros(n);
    resultsMatrix = zeros(n, 1);
    
    % fill the matrices
    for row = 1 : n
        % create symbolic function
        symbolic_expression = str2sym(string(equations(row)));
        % extract coeffiecients
        [coefficients,terms] = coeffs(symbolic_expression);
        % add coefficients to the matrices
        [coeffsMatrix, resultsMatrix, error] = appendCoeffs(coeffsMatrix, resultsMatrix, symbols, coefficients, terms, row);
        if(error)
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
    error = false;
    for column = 1 : length(terms)
        if ~isequaln(sym('0'), terms(j))
            error = true;
            break;
        end
    end    
end
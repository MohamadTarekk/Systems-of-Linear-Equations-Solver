function [X, isSingular] = gaussElimination(coefficients, results, n, tolerance)
    % coefficients is the coeffecients square matrix
    % results is the results matrix
    % n is the number of equations
    % tolerance is the smallest allowable scaled pivot
    % X is the solution matrix
    % isSingular is a boolean indicating if the system is singular, ie not solvable
    scaledFactors = zeros(n, 1);    % an n-element array for storing scaling factors
    for i = 1 : n
        scaledFactors(i) = abs(coefficients(i, 1));
        for j = 2 : n
            if abs(coefficients(i, j)) > scaledFactors(i)
                scaledFactors(i) = abs(coefficients(i, j));
            end
        end
    end
    
    [coefficients, results, isSingular] = forwardElimination(coefficients, results, scaledFactors, n, tolerance);
    if ~isSingular
        [X] = backwardSubstitution(coefficients, results, n);
    end
    disp(isSingular);
end

function [coefficients, results, isSingular] = forwardElimination(coefficients, results, scaledFactors, n, tolerance)
    for row = 1 : n - 1
        [coefficients, results, scaledFactors] =  pivot(coefficients, results, scaledFactors, n, row);
        if abs(coefficients(row, row) / scaledFactors(row))  < tolerance
            isSingular = true;
            return;
        end
        % forward elimination
        for i = row + 1 : n
            factor = coefficients(i, row) / coefficients(row, row);
            for j = row + 1 : n
                coefficients(i, j) = coefficients(i, j) - factor * coefficients(row, j);
            end
            results(i) = results(i) - factor * results(row);
        end
    end
    if abs(coefficients(n, n) / scaledFactors(n)) < tolerance
        isSingular = true;
        return;
    end
    isSingular = false;
end

function [coefficients, results, scaledFactors] =  pivot(coefficients, results, scaledFactors, n, row)
    pivotRow = row;
    % Find the largest scaled coefficient in column k
    max = abs(coefficients(row, row) / scaledFactors(row));
    for i = row + 1 : n
        temp = abs(coefficients(i, row) / scaledFactors(i));
        if temp > max
            max = temp;
            pivotRow = i;
        end
    end
    
    % swapping the rows
    if pivotRow ~= row
        for j = row : n
            temp = coefficients(pivotRow, j);
            coefficients(pivotRow, j) = coefficients(row, j);
            coefficients(row, j) = temp;
        end
        
        temp = results(pivotRow);
        results(pivotRow) = results(row);
        results(row) = temp;
        
        temp = scaledFactors(pivotRow);
        scaledFactors(pivotRow) = scaledFactors(row);
        scaledFactors(row) = temp;
    end
end

function [X] = backwardSubstitution(coefficients, results, n)
    X(n) = results(n) / coefficients(n, n);
    for i = n - 1 : -1 : 1
        sum = 0;
        for j = i + 1 : n
            sum = sum + coefficients(i, j) * X(j);
        end
        X(i) = (results(i) - sum) / coefficients(i, i);
    end
end

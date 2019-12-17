function [X, isSingular] = gaussJordanElimination(coefficients, results, n, tolerance)
    % coefficients is the coeffecients square matrix
    % results is the results matrix
    % n is the number of equations
    % tolerance is the smallest allowable scaled pivot
    % X is the solution matrix
    % isSingular is a boolean indicating if the system is singular, ie not solvable
    [coefficients, results, isSingular] = elimination(coefficients, results, n, tolerance);
    if ~isSingular
        [X] = substitution(coefficients, results, n);
    end
end

function [coefficients, results, isSingular] = elimination(coefficients, results, n, tolerance)
    for row = 1 : n
        [coefficients, results] =  pivot(coefficients, results, n, row);
        if abs(coefficients(row, row))  < tolerance
            isSingular = true;
            return;
        end
        % forward elimination
        range = 1 : n;
        for i = range(range ~= row)
            factor = coefficients(i, row) / coefficients(row, row);
            for j = row + 1 : n
                coefficients(i, j) = coefficients(i, j) - factor * coefficients(row, j);
            end
            results(i) = results(i) - factor * results(row);
        end
    end
    isSingular = false;
end

function [coefficients, results] =  pivot(coefficients, results, n, row)
    pivotRow = row;
    % Find the largest coefficient in column k
    max = abs(coefficients(row, row));
    for i = row + 1 : n
        temp = abs(coefficients(i, row));
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
    end
end

function [X] = substitution(coefficients, results, n)
    X = zeros(1, n);
    for i = 1 : n
        X(i) = results(i) / coefficients(i, i);
    end
end

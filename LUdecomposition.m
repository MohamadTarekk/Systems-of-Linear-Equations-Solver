function [X, isSingular] = LUdecomposition(A, B, n, tolerance)
    % Assume: AX = LUX = B
    % A: 2-D (Sqaure) matrix of Coefficients
    % B: 1-D vector that contains RHS of the equations
    % n: Dimension of the system of equations
    % X: 1-D vector to store the results
    % tolerance: the smallest allowable scaled pivot
    % isSingular: returns true if matrix is singular, ie not solvable
    [scalingFactors] = getScalingFactors(A, n);
    [L, U, B, isSingular] = decompose(A, B, n, scalingFactors, tolerance);
    if(isSingular)
        X = [];
    else
        [X] = substitute(L, U, B, n);
    end
end

function [scalingFactors] = getScalingFactors(a, n)
    scalingFactors = zeros(n, 1);
    for i = 1 : n
        scalingFactors(i) = abs(a(i, 1));
        for j = 2 : n
            if abs(a(i, j)) > scalingFactors(i)
                scalingFactors(i) = abs(a(i, j));
            end
        end
    end
end

function [L, A, B, isSingular] = decompose(A, B, n, scalingFactors, tolerance)
    % use A to store new values of U

    % initialize L matrix
    L = zeros(n);
    for i = 1 : n
        L(i,i) = 1;
    end
    % decompose A into U and L
    for row = 1 : n - 1
        % apply partial pivoting
        [L, A, B, scalingFactors] =  pivot(L, A, B, scalingFactors, n, row);
        % check for sigularity
        if abs(A(row, row) / scalingFactors(row))  < tolerance
            isSingular = true;
            return;
        end
        % forward elimination for A to get U
        % the multiplied factors to get U from A are used to populate L
        for i = row + 1 : n
            factor = A(i, row) / A(row, row);
            for j = row + 1 : n
                A(i, j) = A(i, j) - factor * A(row, j);
            end
            B(i) = B(i) - factor * B(row);
            L(i, row) = factor;
        end
    end
    % check for sigularity
    if abs(coefficients(n, n) / scalingFactors(n)) < tolerance
        isSingular = true;
        return;
    end
    isSingular = false;
end

function [L, U, B, scalingFactors] =  pivot(L, U, B, scalingFactors, n, row)
    pivotRow = row;
    % Find the largest scaled coefficient in column k
    max = abs(U(row, row) / scalingFactors(row));
    for i = row + 1 : n
        temp = abs(U(i, row) / scalingFactors(i));
        if temp > max
            max = temp;
            pivotRow = i;
        end
    end
    % swapping the rows
    if pivotRow ~= row
        % swap U
        for j = row : n
            temp = U(pivotRow, j);
            U(pivotRow, j) = U(row, j);
            U(row, j) = temp;
        end
        % swap L
        for j = 1 : row-1
            temp = L(pivotRow, j);
            L(pivotRow, j) = L(row, j);
            L(row, j) = temp;
        end
        % swap B
        temp = B(pivotRow);
        B(pivotRow) = B(row);
        B(row) = temp;
        % swap scaling factors
        temp = scalingFactors(pivotRow);
        scalingFactors(pivotRow) = scalingFactors(row);
        scalingFactors(row) = temp;
    end
end

function [X] = substitute(L, U, B, n)
    [Y] = forwardSubstitute(L, B, n);
    [X] = backwardSubstitute(U, Y, n);
end

function [Y] = forwardSubstitute(L, B, n)
    Y = zeros(n, 1);
    Y(1) = B(1) / L(1,1);
    for i = 2 : n
        sum = 0;
        for j = 1 : i - 1
            sum = sum + L(i, j) * Y(j);
        end
        Y(i) = (B(i) - sum) / L(i,i);
    end
end

function [X] = backwardSubstitute(U, Y, n)
    X(n) = Y(n) / U(n, n);
    for i = n - 1 : -1 : 1
        sum = 0;
        for j = i + 1 : n
            sum = sum + U(i, j) * X(j);
        end
        X(i) = (Y(i) - sum) / U(i, i);
    end
end
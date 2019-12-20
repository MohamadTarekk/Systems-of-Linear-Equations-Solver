function [X, iterations, data] = gaussSeidel(coefficients, results, initialGuesses, n, maxIterations, tolerance)
    data = zeros(maxIterations + 1, n);
    for i = 1 : n
        data(1, i) = initialGuesses(i);
    end
    index = ones(1, n);
    X = zeros(1, n);
    precision = zeros(maxIterations, n);
    tolerance = tolerance * 100;
    
    iterations = maxIterations;
    for i = 1 : maxIterations
        stop = true;
        for j = 1 : n
            sum = 0;
            range = 1 : n;
            for k = range(range ~= j)
                sum = sum + coefficients(j, k) * data(index(k), k);
            end
            index(j) = index(j) + 1;
            data(index(j), j) = (results(j) - sum) / coefficients(j, j);
            precision(i, j) = abs( ( data(index(j), j) - data(index(j) - 1, j) ) / data(index(j), j) ) * 100;
            stop = stop && (precision(i, j) < tolerance);
        end
                
        if stop
            iterations = i;
            break;
        end
    end
    
    for i = 1 : n
        X(i) = data(iterations + 1, i);
    end
    
    % remove initial guess row (first row)
    data(1,:) = [];
    % concatinate the precision
    result = zeros(iterations, 2 * n);
    for i = 1 : n
        result(1:iterations,2 * i - 1) = data(1:iterations,i);
        result(1:iterations,2 * i) = precision(1:iterations,i);
    end
    data = result;
end

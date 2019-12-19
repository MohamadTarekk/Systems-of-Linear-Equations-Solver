function [X, iterations, data, precision] = gaussSeidel(coefficients, results, initialGuesses, n, maxIterations, tolerance)
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
end

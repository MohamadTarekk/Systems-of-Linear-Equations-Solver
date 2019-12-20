% sample run
[error, n, method, symbols, equations, initialConditions, max_iter, epsilon, isIterative] = readfile();
if error == 0
    [error] = validateInputs(n, method, initialConditions, max_iter, epsilon, isIterative);
    if error == 0
        [error, n, initialConditions, max_iter, epsilon, coeffsMatrix, resultsMatrix, symbols] = parseInput(n, isIterative, initialConditions, max_iter, epsilon, equations, symbols);
        if error == 0
            [X, isSingular] = LUdecomposition(coeffsMatrix, resultsMatrix, n, 0.000001);
        else
            errordlg(error);
        end
    else
        errordlg(error);
    end
else
    errordlg(error);
end
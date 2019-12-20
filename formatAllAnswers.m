function [result] = formatAllAnswers(Xgauss, XgaussJordan, XluDecomp, XgaussSeidel, variables)
    result = "";
    
    n = length(variables);

    result = result + sprintf("Gauss Elimination:\n");
    for i = 1 : n
        result = result + sprintf("%s = %s\n", string(variables(i)), string(Xgauss(i)));
    end

    result = result + sprintf("Gauss-Jordan Elimination:\n");
    for i = 1 : n
        result = result + sprintf("%s = %s\n", string(variables(i)), string(XgaussJordan(i)));
    end

    result = result + sprintf("LU Decomposition:\n");
    for i = 1 : n
        result = result + sprintf("%s = %s\n", string(variables(i)), string(XluDecomp(i)));
    end

    result = result + sprintf("Gauss-Seidel Method:\n");
    for i = 1 : n
        result = result + sprintf("%s = %s\n", string(variables(i)), string(XgaussSeidel(i)));
    end
end


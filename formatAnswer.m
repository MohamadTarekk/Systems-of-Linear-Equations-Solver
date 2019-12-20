function [result] = formatAnswer(X, variables)
    result = "";
    for i = 1 : length(X)
        result = result + sprintf("%s = %s\n", string(variables(i)), string(X(i)));
    end
end

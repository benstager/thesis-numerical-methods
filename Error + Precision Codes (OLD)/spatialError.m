function [s_error] = spatialError(approx,reference)

if ~ isValidGridSize(length(reference))
    error('Invalid reference size');
end

if ~isValidGridSize(length(approx))
    error('Invalid approx size');
end

nr = log(length(reference)+1)/log(2);
na = log(length(approx)+1)/log(2);
gamma = nr - na;

restrictedRef = reference(2^gamma:2^gamma:end);

s_error = norm(restrictedRef-approx,'inf');
end

function flag = isValidGridSize(n)

exponent = round(log(n+1)/log(2));
flag = 2^exponent == n+1;

end
function [P,sgP] = LinApproximator(y,r,funcs)
N = size(y, 2);
M  = size(funcs, 1);

g = zeros(N, M);
for i = 1 : N
    for j = 1 : M
        f = cell2mat(funcs(j));
        vec = num2cell(r(:, i));
        g(i, j) =  f(vec{:});
    end
end

P = g\(y');
Err_0 = sqrsum(y - g*P)/N;

Gk = trace(G);
sgP = sqrt(Err_0/(2*Gk));

end

function[sum] = sqrsum(x)
    len = length(x);
    sum = 0;
    for k = 1 : len
        sum = sum + x(k)^2;
    end
end
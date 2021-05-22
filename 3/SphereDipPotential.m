function [F, X, Y, P] = SphereDipPotential(XYZ, Q, D, R, r0, a, b, Dx, Dy, Nxy)

N = length(R);
F = zeros(Nxy(1),Nxy(2));
b = b - a .* (a(1) * b(1) + a(2) * b(2) + a(3) * b(3)) ./ (a(1) * a(1) + a(2) * a(2) + a(3) * a(3));
P = [a(1) b(1); a(2) b(2); a(3) b(3)];
X = [Dx(1):(Dx(2) - Dx(1)) / (Nxy(2) -1):Dx(2)];
Y = [Dy(1):(Dy(2) - Dy(1)) / (Nxy(1) -1):Dy(2)];

for i = 1:Nxy(2)
    for j = 1:Nxy(1)
        x0 = r0(1) + a(1) * X(i) + b(1) * Y(j);
        y0 = r0(2) + a(2) * X(i) + b(2) * Y(j);
        z0 = r0(3) + a(3) * X(i) + b(3) * Y(j);
        oldr = [x0; y0; z0];
        for k = 1:N
            newr = oldr - XYZ(:, k);
            if(norm(newr) > R(k))
                F(j,i) = F(j,i) + dot(newr, D(k,:)') / norm(newr)^3;
            else
                F(j,i) = F(j, i) + Q(k) / R(k) + dot(newr, D(k,:)') / R(k)^3;
            end
        end
    end
end      

X = repmat(X, Nxy(1), 1);
Y = Y';
Y = repmat(Y, 1, Nxy(2));
end



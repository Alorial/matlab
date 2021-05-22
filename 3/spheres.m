R = [5 2]; 
F = [1; 2];
XYZ = [5 20; -10 1; 10 -5];
Nxy = [100 200];
a = [0; 1; 0];
b = [1; 1; 0];
Dx = [-10 10];
Dy = [-10 10];
r0 = [1;1;1];

Q = ElectroStaticBalls(XYZ, R, F);
disp('Q = ');
disp( Q);
[F,X,Y,P] = SpherePotential(XYZ,Q,R,r0,a,b,Dx,Dy,Nxy)

figure; hold on; grid on; mesh(X,Y,F); 

function [Q] = ElectroStaticBalls(XYZ, R, F)
N = length(R);
dist = zeros(N,N);

for i = 1 : N 
    for j = 1 : N 
        dist(i,j) = sqrt((XYZ(1,i) - XYZ(1,j))^2+((XYZ(2,i) - XYZ(2,j))^2+(XYZ(3,i) - XYZ(3,j))^2));
    end
end
    function [L] = check(XYZ, R)
        N = length(R);
        L = 1;
        for i = 1:(N - 1)
            for j = (i + 1):N
                if norm(XYZ(:,i)-XYZ(:,j)) <= R(i) + R(j)
                    error('WRONG COORDS');
                    L = 0;
                end
            end
        end
    end

if check(XYZ, R) == 1
    disp('rel_dist: ');
    disp(dist);
    M = dist .^ (-1);
    for i = 1:N
        M(i,i) = 1/R(i);
    end
    disp('M: ');
    disp(M);
    Q = F\M;
end
end


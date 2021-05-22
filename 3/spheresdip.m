global k
k = 9*10^9;
R = [5 2];
F = [k*1; k*20];
XYZ = [10 5; -7 -1; 0 0];
Nxy = [100 200];
a = [0; 1; 0];
b = [1; 1; 0];
Dx = [-10 10];
Dy = [-10 10];
r0 = [1;1;1];

[Q,D] = ElectroStaticDipoles(XYZ,R,F);
disp('Q = ');
disp( Q);
disp('D = ');
disp( D);

[F,X,Y,P] = SphereDipPotential(XYZ,Q,D,R,r0,a,b,Dx,Dy,Nxy);
figure; hold on; grid on; mesh(X,Y,F); 

function [Q,D] = ElectroStaticDipoles(XYZ,R,F)
global k
N = length(R);
dist = zeros(N,N);
for i = 1 : N
    for j = 1 : N 
        dist(i,j) =   ((XYZ(1,i) - XYZ(1,j))^2+((XYZ(2,i) - XYZ(2,j))^2+(XYZ(3,i) - XYZ(3,j))^2  )) .^ 0.5;
    end
end
    function [L] = check(XYZ, R)
        N = length(R);
        L = 1;
        for i = 1:(N - 1)
            for j = (i + 1):N
                if norm(XYZ(:,i)-XYZ(:,j)) <= R(i) + R(j)
                    error('WRONG COORDS');
                end
            end
        end
    end

Rx = zeros(N, N);
Ry = zeros(N, N);
Rz = zeros(N, N);

for i = 1 : N 
    for j = 1 : N 
        Rx(i,j) =   abs((XYZ(1,i) - XYZ(1,j)));
    end
end

for i = 1 : N 
    for j = 1 : N 
        Ry(i,j) =   abs((XYZ(2,i) - XYZ(2,j)));
    end
end

for i = 1 : N 
    for j = 1 : N
        Rz(i,j) =   abs((XYZ(3,i) - XYZ(3,j)));
    end
end

Res = zeros(4*N, 1);
Res(1:N) = F;

bigMatrix = zeros(4*N, 4*N);
pointPot = 1 ./ (dist);
for i = 1: N
    pointPot(i, i) = 0;
end
dipPotx =  Rx ./ (dist .^ 2);
for i = 1: N
    dipPotx(i, i) = 0;
end
dipPoty =  Ry ./ (dist .^ 2);
for i = 1: N
    dipPoty(i, i) = 0;
end
dipPotz =  Rz ./ (dist .^ 2);
for i = 1: N
    dipPotz(i, i) = 0;
end
pointEx = Rx ./ (dist .^ 3);
for i = 1: N
    pointEx(i, i) = 0;
end

pointEy = Ry ./ (dist .^ 3);

for i = 1: N
    pointEy(i, i) = 0;
end
pointEz = Rz ./ (dist .^ 3);
for i = 1: N
    pointEz(i, i) = 0;
end
dipExx = (3* Rx .^ 2) ./ (dist .^ 5) - 1 ./ (dist .^ 3);
for i = 1: N
    dipExx(i, i) = 0;
end
dipExy = (3* Rx .* Ry)  ./ (dist .^ 5);
for i = 1: N
    dipExy(i, i) = 0;
end
dipEyz = (3* Ry .* Rz) ./ (dist .^ 5);
for i = 1: N
    dipEyz(i, i) = 0;
end
dipEzx = (3* Rz .* Rx)./ (dist .^ 5);
for i = 1: N
    dipEzx(i, i) = 0;
end
dipEyy = (3* Ry .^ 2 )./ (dist .^ 5) - 1 ./ (dist .^ 3);
for i = 1: N
    dipEyy(i, i) = 0;
end
dipEzz = (3* Rz .^2 )./ (dist .^ 5) - 1 ./ (dist .^ 3);
for i = 1: N
    dipEzz(i, i) = 0;
end

bigMatrix(1:N, 1:N) = pointPot;

bigMatrix(1:N, N+1:2*N) = dipPotx;
bigMatrix(1:N, 2*N+1:3*N) = dipPoty;
bigMatrix(1:N, 3*N+1:4*N) = dipPotz;

bigMatrix(N+1: 2*N, 1:N) = pointEx;
bigMatrix(2*N+1: 3*N, 1:N) = pointEy;
bigMatrix(3*N+1: 4*N, 1:N) = pointEz;

bigMatrix(N+1:2*N, N+1:2*N) = dipExx;
bigMatrix(N+1:2*N, 2*N+1:3*N) = dipExy;
bigMatrix(N+1:2*N, 3*N+1:4*N) = dipEzx;

bigMatrix(2*N+1:3*N, N+1:2*N) = dipExy;
bigMatrix(2*N+1:3*N, 2*N+1:3*N) = dipEyy;
bigMatrix(2*N+1:3*N, 3*N+1:4*N) = dipEyz;

bigMatrix(3*N+1:4*N, N+1:2*N) = dipEzx;
bigMatrix(3*N+1:4*N, 2*N+1:3*N) = dipExy;
bigMatrix(3*N+1:4*N, 3*N+1:4*N) = dipEzz;
disp(bigMatrix);

A = zeros(4*N, 1);
A = Res \ bigMatrix;
A = A * k;

disp(A);
Q = A(1:N);
D = zeros(3, N);
for i =1 : N
    D(1, i) = A(N+i);
    D(2, i) = A(2*N+i);
    D(3, i) = A(3*N+i);
end
D = D';

end
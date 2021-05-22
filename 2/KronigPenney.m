function [E] = KronigPenney(k, m0, a0, b0, U00, Emax)
global m;
global a;
global b;
global U0;
global int_k;
int_k = k(1);
eV = 1.602176634e-19;
%переводим в нужные ед.измерения
a = a0*10^(-9); 
b = b0*10^(-9); 
U0 = U00*eV;
m = m0;
Eees = U00+0.01:0.001:Emax;
ff = Eees;

for i = 1:(length(Eees))
    Ee = Eees(i)*eV;
    ff(i) = funcf(m, a, b, U0, Ee);
end

figure, hold on, grid on, plot(ff, Eees);

extr = U00+0.01;
dif = diff(ff);
for j = 1:(length(dif)-1)
    if (sign(dif(j)) * sign(dif(j+1))) <= 0
        extr = [extr, Eees(j)];
    end
end

fk = @fun;
k_len= length(k);
for q = 1:(length(extr)-1)
    E_int = [extr(q), extr(q+1)];
    for qk = 1:k_len
        int_k = k(qk);
        if (sign(fun(E_int(1))*sign(fun(E_int(end))))) < 0
            E(q, qk) = fzero(fk, E_int);
        end
    end
end

end

function f = funcf(m, a, b, U0, E)
h = 1.054571817e-34;
mu_q = 2*m*E/(h*h);
lambda_q = 2*m*(E-U0)/(h*h);
mu = sqrt(mu_q);
lambda = sqrt(lambda_q);
f = cos(mu*a)*cos(lambda*b) - (lambda_q + mu_q)/(2*mu*lambda)*sin(mu*a)*sin(lambda*b);
end

function fk = fun(Ee)
eV = 1.602176634e-19;
global int_k;
global m;
global a;
global b;
global U0;
h = 1.054571817e-34;
E = Ee*eV;

mu_q = 2*m*E/(h*h);
lambda_q = 2*m*(E-U0)/(h*h);
mu = sqrt(mu_q);
lambda = sqrt(lambda_q);

fk = cos(mu*a)*cos(lambda*b) - (lambda_q + mu_q)/(2*mu*lambda)*sin(mu*a)*sin(lambda*b) - cos(int_k*(a+b));
end
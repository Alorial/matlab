U_0 = -4;
Emax = 25;
m = 9.109383701528 * 10^(-31);
a = 0.5;
b = 2;
k = -0.6e+9: 0.01e+9: 0.6e+9;
E = KronigPenney(k, m, a, b, U_0, Emax);
figure; hold on; grid on; 
for q = 1:length(k)
    plot(k(q), E(:, q), 'linestyle', 'none', 'marker', '*');
end
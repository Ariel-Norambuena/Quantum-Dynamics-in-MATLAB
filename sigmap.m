function x = sigmap(i,L,Is,Icav)
up = [1 0]';
down = [0 1]';
sigma_p = up*down';
sigma_p = kron(sigma_p,Icav);
Op_total = cell(1,L);
for site = 1:L
    Op_total{site} = Is+double(eq(i,site))*(sigma_p-Is);
end
x = Op_total{1};
for site = 2:L
    x = kron(x,Op_total{site});
end
end
function x = acav(i,L,Nphoton,Is,Iatom)
a = diag(sqrt(1:Nphoton)',1);
a = kron(Iatom,a);
Op_total = cell(1,L);
for site = 1:L
    Op_total{site} = Is+double(eq(i,site))*(a-Is);
end
x = Op_total{1};
for site = 2:L
    x = kron(x,Op_total{site});
end
end
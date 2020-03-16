function [R_sort,L_sort,lambda_sort] = sortingEigenvalues(dim,TOL,L)

[R,DR] = eig(L);       % Right eigenvectors and eigenvalues
[L,DL] = eig(L');      % Left eigenvectors and eigenvalues
eig_R = diag(DR);      % Right eigenvalues written as a vector
eig_L = diag(DL);      % Left eigenvalues written as a vector
ind_RL = zeros(dim*dim,2);
count = 1;
for n=1:dim*dim                 % Sorting of eigenvalues
    an = eig_R(n);
    for m=1:dim*dim
        bm = eig_L(m);
        if(abs(real(an)-real(bm))<TOL && abs(imag(an)-imag(bm))<TOL && count<=dim*dim)
            ind_RL(count,1) = n;
            ind_RL(count,2) = m;
            count = count + 1;
        end
    end
end
eig_L = eig_L(ind_RL(:,2)');    % Final sorting
eig_R = eig_R(ind_RL(:,1)');
L = L(:,ind_RL(:,2)');
R = R(:,ind_RL(:,1)');
lambda = eig_R;
[~,ind] = sort(lambda);         % Sorting of eigenvalues
lambda_sort = lambda(ind);           % \lambda_k eigenvalues
L = L(:,ind);
R = R(:,ind);

% Final R_sort and L_sort matrices
R_sort = cell(1,length(lambda_sort));
L_sort = cell(1,length(lambda_sort));
for k=1:length(lambda_sort)
    R_sort{k} = reshape(R(:,k),dim,dim);
    L_sort{k} = reshape(L(:,k),dim,dim);
    Rk = R_sort{k};
    Lk = L_sort{k};
    Ck = trace(Lk*Rk);
    Lk = Lk/sqrt(Ck);          % Normalized left eigenmatrices
    Rk = Rk/sqrt(Ck);          % Normalized right eigenmatrices
    R_sort{k} = Rk;
    L_sort{k} = Lk;
end
end
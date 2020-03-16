dim = 4;                  % Dimension of the total Hilbert space
Is = eye(dim);            % Identity matrix for the total Hilbert space
J = 1;                    % Value of the coupling term J
B = 0.1*J;                % Value of the magnetic field B
Sx = [0 1;1 0];    	  % S_x operator for one spin
Sz = [1 0;0 -1];          % S_z operator for one spin
I = eye(2);               % Identity matrix for one spin 1/2
H = -J*kron(Sx,Sx)-B*(kron(Sx,I)+kron(I,Sx)); % Hamiltonian of the system
L_H = -1i*kron(Is,H)+1i*kron(H.',Is);         % Lindblad operator L_H
gamma_1 = 0.1*B;       % Decay rate gamma_1
gamma_2 = 0.5*B;       % Decay rate gamma_2
S_minus =[0 0; 1 0];   % Lowering operator of the particle 1
L1 = kron(S_minus,I);  % Lowering operator of the particle 1 in the total Hilbert space
DL_1 = gamma_1*(kron(conj(L1),L1)-0.5*kron(Is,L1'*L1)-0.5*kron(L1.'*conj(L1),Is));
L2 = kron(I,S_minus);  % Lowering operator of the particle 2 in the total Hilbert space
DL_2 = gamma_2*(kron(conj(L2),L2)-0.5*kron(Is,L2'*L2)-0.5*kron(L2.'*conj(L2),Is));
L_diss = DL_1 + DL_2;  % Lindbladian L_diss
L = L_H + L_diss;      % Total Lindblad operator
TOL = 1e-10;
[R_sort,L_sort,lambda_sort] = sortingEigenvalues(dim,TOL,L);

down = [0 1]';                    % Quantum state down = [0 1]
Psi_0 = kron(down,down);          % Initial wavefunction
rho_0 = Psi_0*Psi_0';             % Initial density matrix
Nt = 3000;     	                  % Number of steps for time
T = 2*pi/(2*B) ;                  % Period of time
ti = 0;                           % Intial time
tf = 4*T;                         % Final time
dt = (tf-ti)/(Nt-1);              % Step for time
t = ti:dt:tf;                     % Time vector
Mz = zeros(size(t));              % Initial average magnetization
SSz = (kron(Sz,I)+kron(I,Sz))/2;  % Operator S1^z+S2^z
for n=1:length(t)        % Iteration to find general sulution of rho(t)
    rho = zeros(dim,dim);
    for k=1:length(lambda_sort)
        Lk = L_sort{k};
        Rk = R_sort{k};
        ck = trace(rho_0*Lk);
        rho = rho + ck*exp(lambda_sort(k)*t(n))*Rk; % General sulution for rho(t)
    end
    Mz(n) = trace(SSz*rho);   % General sulution for Mz(t)
end
figure()
hold on
plot(t/T,real(Mz),'r-','LineWidth',3)
plot(t/T, exp(real(lambda_sort(end))*t),'k--','LineWidth',2)
plot(t/T,-exp(real(lambda_sort(end))*t),'k--','LineWidth',2)
hold off
xlabel('$Bt$','Interpreter','LaTex','Fontsize', 30)
ylabel('$\langle M_z \rangle$','Interpreter','LaTex','Fontsize', 30)
set(gca,'fontsize',21)
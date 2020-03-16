Omega = 1;                      % Raby frequency from |1> to |2>
gamma0 = 0.2*Omega;             % Decay rate
dim = 2;                        % Dimension Hilbert space two-level system
Is = eye(dim);                  % Identity matrix Hilbert space
v1 = Is(:,1);                   % Excited state for the atom
v2 = Is(:,2);                   % Ground state fro the atom
s11 = v1*v1';                   % Atom operator sigma_{11}
s22 = v2*v2';                   % Atom operator sigma_{22}
sp = v1*v2';                    % Atom operator sigma_{+}
sm = v2*v1';                    % Atom operator sigma_{-}
HL = -0.5*Omega*(sp+sm);        % Hamiltonian of the two-level system
N = 0;                          % Mean number of photons at zero temperature
Lrad = gamma0*(N+1)*(kron(conj(sm),sm)-0.5*kron(Is,sm'*sm)-0.5*kron(sm.'*conj(sm),Is)) + ...
    gamma0*N*(kron(conj(sp),sp)-0.5*kron(Is,sp'*sp)-0.5*kron(sp.'*conj(sp),Is));
L = -1i*kron(Is,HL)+1i*kron(HL.',Is)+Lrad;  % Lindblad operator
TOL = 1e-6;
[R_sort,L_sort,lambda_sort] = sortingEigenvalues(dim,TOL,L);
psi_0 = v2;                     % Ground state
rho_0 = psi_0*psi_0';           % Initial density matrix
Nt = 100000;                    % Number of steps for time
ti = 0;                         % Initial time
tf = 50/Omega;                  % Final time
dt = (tf-ti)/(Nt-1);            % Step time dt
t = ti:dt:tf;                   % Time vector
pe = zeros(size(t));         % Matrix elements \rho_{ij} = <i|\rho|j>
sigmap = zeros(size(t));
for n=1:Nt                      % General solution
    rho = zeros(dim,dim);
    for k=1:length(lambda_sort)
        Lk = L_sort{k};
        Rk = R_sort{k};
        ck = trace(rho_0*Lk);
        rho = rho + ck*exp(lambda_sort(k)*t(n))*Rk;
    end
    pe(n) = rho(1,1);
    sigmap(n) = trace(rho*sp);
end
% Exact Solution at zero temperature
% p_e: population excited state, sp = <\sigma_x>
mu = 1i*sqrt((gamma0/4)^2-Omega^2);
pe_exact = Omega^2/(gamma0^2+2*Omega^2)*(1-exp(-3*gamma0*t/4).*(cos(mu*t)+3*gamma0/4/mu*sin(mu*t)));
if gamma0~=0
    sp_exact = -1i*Omega*gamma0/(gamma0^2+2*Omega^2)*(1-exp(-3*gamma0*t/4).*(cos(mu*t)+(gamma0/4/mu)-Omega^2/gamma0/mu*sin(mu*t)));
else
    sp_exact = -1i/2*sin(Omega*t);
end
figure
box on
hold on
plot(t*Omega,real(pe),'r-','Linewidth',2)
plot(t*Omega,pe_exact,'b--','Linewidth',2)
hold off
xlabel('$\Omega t$','Interpreter','LaTex','Fontsize', 30)
ylabel('$p_e(t)$','Interpreter','LaTex','Fontsize', 30)
legend({'$\mbox{Numerical}$','$\mbox{Exact}$'},'Interpreter','latex','Fontsize', 21,'Location','northeast')
set(gca,'fontsize',21)
xlim([0 50])

figure
box on
hold on
plot(t*Omega,imag(sigmap),'r-','Linewidth',2)
plot(t*Omega,imag(sp_exact),'b--','Linewidth',2)
hold off
xlabel('$\Omega t$','Interpreter','LaTex','Fontsize', 30)
ylabel('$\mbox{Im}\langle \hat{\sigma}_+ \rangle$','Interpreter','LaTex','Fontsize', 30)
legend({'$\mbox{Numerical}$','$\mbox{Exact}$'},'Interpreter','latex','Fontsize', 21,'Location','northeast')
set(gca,'fontsize',21)
xlim([0 50])
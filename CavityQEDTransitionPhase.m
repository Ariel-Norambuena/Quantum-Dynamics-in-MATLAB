Nsim = 25;                      % Number of simulations 
L = 2;                          % Number of cavities
wc = 1;                         % Cavity frequency
g = 1e-2*wc;                    % Atom-light coupling
J = 1e-4*wc;                    % Coupling between cavities
Nph = 2;                        % Number of photons per cavity
dimFock = Nph+1;                % Dimension Fock space for photons
dimT = 2*dimFock;               % Dimension atom+cavity system
Deltai = 10^(-2)*g;             % Initial detuning
Deltaf = 10^(+2)*g;             % Final detuning
xi = log10(Deltai/g);           % Initial detuning in log scale
xf = log10(Deltaf/g);           % Final detuning in log scale
dx = (xf-xi)/(Nsim-1);          % Step dx 
x = xi:dx:xf;                   % Vector x to plot transition phase
OP_JC = zeros(size(x));         % Order parameter Jaynes-Cummings model
OP_R = zeros(size(x));          % Order parameter Rabi model
A = cell(1,L);                  % Cell array to storage \hat{a}_i operators
Sp = cell(1,L);                 % Cell array to storage \hat{\sigma}_i^+ operators
N_ex = cell(1,L);               % Cell array to storage N_i operators
Iatom = eye(2);                 % Identity matrix atom system
Icav = eye(dimFock );           % Identity matrix cavity system
Is = eye(2*dimFock);            % Identity matrix atom+cavity system
for i=1:L
    A{i} = acav(i,L,Nph,Is,Iatom);     % Photon operator \hat{a}_i
    Sp{i} = sigmap(i,L,Is,Icav);       % Atom operator \hat{\sigma}_i 
    N_ex{i} = A{i}'*A{i}+Sp{i}*Sp{i}'; % number of excitations per cavity
end  
Ad = ones(L);                   % Adyacent matrix
Ad = triu(Ad)-eye(L);           % A(i,j)=1 for j>i
Hhopp = zeros(dimT^L,dimT^L);   % Interaction Hamiltonian
for i=1:L
    for j=1:L
        Hhopp = Hhopp - J*Ad(i,j)*A{i}'*A{j} - J*Ad(i,j)*A{i}*A{j}';
    end
end
Nt = 10000;                     % Length time vector 
ti = 0.01/J;                    % Initial time
tf = 1/J;                       % Final time
dt = (tf-ti)/(Nt-1);            % Step time dt
t = ti:dt:tf;                   % Time vector 
Lambda_R = zeros(Nsim,Nt);
Lambda_JC = zeros(Nsim,Nt);
parfor n=1:Nsim
    D = g*10^(x(n));            % Detuning in each simulation
    Model = 'Rabi';             
    [OP_R(n), P1m_R] = QuantumSimulationCavityArray(wc,D,g,Nph,L,A,Sp,N_ex,Hhopp,t,Model);
    Model = 'Jaynes-Cummings';  
    [OP_JC(n),P1m_JC] = QuantumSimulationCavityArray(wc,D,g,Nph,L,A,Sp,N_ex,Hhopp,t,Model);
    Lambda_R(n,:) = -1/L*log2(P1m_R);     % Rate function for the Rabi model
    Lambda_JC(n,:) = -1/L*log2(P1m_JC);   % Rate function for the Jaynes-Cummings model
end  
figure()
box on
hold on
plot(J*t,Lambda_R(end,:),'b-','Linewidth',3)
plot(J*t,Lambda_JC(end,:),'r-','Linewidth',3)
hold off
xlabel('$Jt$','Interpreter','LaTex','Fontsize', 30)
ylabel('$\Lambda(t)$','Interpreter','LaTex','Fontsize', 30)
set(gca,'fontsize',21)
legend({'$\mbox{RH}$','$\mbox{JCH}$'},'Interpreter','latex','Fontsize', 21,'Location','best')

figure()
box on
hold on
plot(x,real(OP_R),'.b','Markersize',30)
plot(x,real(OP_JC),'.r','Markersize',30)
hold off
xlabel('$\mbox{Log}_{10}(\Delta/g)$','Interpreter','LaTex','Fontsize', 30)
ylabel('$\mbox{OP}$','Interpreter','LaTex','Fontsize', 30)
set(gca,'fontsize',21)
legend({'$\mbox{RH}$','$\mbox{JCH}$'},'Interpreter','latex','Fontsize', 21,'Location','best')
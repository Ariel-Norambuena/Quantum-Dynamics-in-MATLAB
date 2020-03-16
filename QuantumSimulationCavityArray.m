function [OP,P1m]= QuantumSimulationCavityArray(wc,D,g,Nphoton,L,A,Sp,N_ex,Hhopp,t,Model)
HJC = zeros((Nphoton+1)^L*2^L,(Nphoton+1)^L*2^L);
for i=1:L
    switch Model
        case 'Jaynes-Cummings'
            HJC = HJC + wc*A{i}'*A{i} + (D+wc)*Sp{i}*Sp{i}' + g*(Sp{i}*A{i}+Sp{i}'*A{i}');
        case 'Rabi'
            HJC = HJC + wc*A{i}'*A{i} + (D+wc)*Sp{i}*Sp{i}' + g*(Sp{i}+Sp{i}')*(A{i}+A{i}');
    end
end
H = HJC + Hhopp;        % Total Hamiltonian
dt = t(2)-t(1);         % Step time dt
U = expm(-1i*H*dt);     % Time propagator with step dt
up = [1 0]';            % Excited state for the two-level system
down = [0 1]';          % Ground state for the two-level system
Fock = eye(Nphoton+1);  % Fock states
theta1 = 0.5*atan(2*g*sqrt(1)/D);
phi_1m = cos(theta1)*kron(down,Fock(:,2))...
    -sin(theta1)*kron(up,Fock(:,1));  % Initial state |1,->
PSI_0 = phi_1m;
for k=1:L-1
    PSI_0 = kron(PSI_0,phi_1m);   % Many body wavefunction |1,->...|1,-> (L terms)
end
dn_T = zeros(size(t));            % Standard deviation dn_i = <n_i^2>-<n_i>^2
P1m = zeros(size(t));             % Population |1,-> at time t
for n=1:length(t)
    if n==1
        PSI = PSI_0;              % Initial wavefunction
    else
        PSI = U*PSI;              % Wavefunction at time t_n
    end
    dn = 0;
    for i=1:L
        dn = dn + PSI'*N_ex{i}^2*PSI-(PSI'*N_ex{i}*PSI)^2;
    end
    P1m(n) = abs(PSI_0'*PSI)^2;   % Ground state probability
    dn_T(n) = dn;                 % Standard deviation of the number of excitations
end
OP = mean(dn_T);  % Order parameter
end
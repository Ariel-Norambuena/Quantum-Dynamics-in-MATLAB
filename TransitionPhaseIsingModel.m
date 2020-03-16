Nspins = 6;   % Number of spins
H0 = zeros(2^Nspins,2^Nspins);
H1 = zeros(2^Nspins,2^Nspins);
J = 0;
alpha = 0.2;  % Parameter \alpha
for j=1:Nspins
    for i=j+1:Nspins
        Jij = abs(i-j)^(-alpha);
        J = J + (Nspins-1)^(-1)*Jij;
    end
end
B = J/0.42;   % Magnetic field
Sx = [0 1;1 0];
Sz = [1 0;0 -1];
for i=1:Nspins
    Szi = getSci(Sz,i,Nspins);
    H0 = H0 - B*Szi;       % Hamiltonian for the magnetic field
    for j=1:Nspins
        if i~=j
            Sxi = getSci(Sx,i,Nspins);
            Sxj = getSci(Sx,j,Nspins);
            Vij = abs(i-j)^(-alpha)/J;
            H1 = H1 - Vij*Sxi*Sxj; % Interaction Hamiltonian
        end
    end
end
xr = [1 1]'/sqrt(2);        % Single-particle state |Psi_{-->}>
xl = [-1 1]'/sqrt(2);       % Single-particle state |Psi_{<--}>
Xr = xr;
Xl = xl;
for n=1:Nspins-1
    Xr = kron(Xr,xr);        % Many-body state |Psi_{-->}>
    Xl = kron(Xl,xl);        % Many-body state |Psi_{<--}>
end
PSI_0 = Xr;                 % Initial condition
Mx = 0;
for i=1:Nspins
    Sxi = getSci(Sx,i,Nspins);
    Mx = Mx + Sxi/Nspins; % Magnetization operator
end
H = H0 + H1;                % Total hamiltonian
ti = 0;                     % Initial time
tf = 22;                    % Final time
Nt = 10000;                 % Number of steps
dt = (tf-ti)/(Nt-1);        % Step time dt
t = ti:dt:tf;               % Time vector
U = expm(-1i*H*dt);         % Time propagator operator U(dt)
Lambda = zeros(size(t));    % Rate function \Lambda(t)
Av_Mx = zeros(size(t));     % Average Magnetization <M_x(t)>
for n=1:length(t)
    if n==1
        PSI = PSI_0;        % Initial wavefunction
    else
        PSI = U*PSI;        % Wavefunction at time t_n
    end
    Pr = abs(Xr'*PSI)^2;    % Probability state |Psi_{-->}>
    Pl = abs(Xl'*PSI)^2;    % Probability state |Psi_{<--}>
    Lambda(n) = min(-Nspins^(-1)*log(Pr),-Nspins^(-1)*log(Pl));
    Av_Mx(n) =  PSI'*Mx*PSI;  % Average magnetization along x-axis
end

figure()
plot(B*t,Lambda,'b-','Linewidth',3)
xlabel('$B t$','Interpreter','LaTex','Fontsize', 30)
ylabel('$\Lambda(t)$','Interpreter','LaTex','Fontsize', 30)
set(gca,'fontsize',21)
xlim([0 5])

figure()
box on
plot(t*B,real(Av_Mx),'r-','Linewidth',2)
xlabel('$B t$','Interpreter','LaTex','Fontsize', 30)
ylabel('$\langle M_x\rangle$','Interpreter','LaTex','Fontsize', 30)
set(gca,'fontsize',21)
xlim([0 100])
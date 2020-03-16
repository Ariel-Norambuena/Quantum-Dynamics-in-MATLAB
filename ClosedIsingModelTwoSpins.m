J = 1;              % Value of the coupling J
B = 0.1*J;          % Value of the magnetic field
Sx = [0 1;1 0];     % S_x operator
Sz = [1 0; 0 -1];   % S_z operator
I = eye(2);         % Identity matrix
Hspins = -J*kron(Sx,Sx)-B*(kron(Sx,I)+kron(I,Sx));  % Hamiltonian
down = [0 1]';                   % Quantum state down = [0 1]^T
Psi_0 = kron(down,down);         % Initial wavefunction
T = 2*pi/(2*B);                  % Period of time
Nt = 1000;                       % Number of steps to construct time vector
ti = 0;                          % Initial time
tf = 2*T;                        % Final time
dt = (tf-ti)/(Nt-1);             % Step time dt
t = ti:dt:tf;                    % Time vector
U = expm(-1i*Hspins*dt);         % Time propagator operator U(dt)
SSz = (kron(Sz,I)+kron(I,Sz))/2; % Operator S1^z + S2^z
Mz = zeros(size(t));             % Average magnetization
for n=1:length(t)                % Iteration to find Psi(t) and Mz(t)
    if n==1
        Psi = Psi_0;             % Initial wavefunction
    else
        Psi = U*Psi;             % Wavefuntion at time t_n
    end
    Mz(n) = Psi'*SSz*Psi;        % Average magnetization at time t_n
end
plot(t/T,real(Mz),'r-','LineWidth',3)
xlabel('$t/T$','Interpreter','LaTex','Fontsize', 30)
ylabel('$\langle M_z \rangle$','Interpreter','LaTex','Fontsize', 30)
set(gca,'fontsize',21)
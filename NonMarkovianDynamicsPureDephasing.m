Nw = 5000;                            % Number of points for \omega
wi = 0.01;                            % Initial frequency
wf = 5;                               % Final frequency
dw = (wf-wi)/(Nw-1);                  % Step frequency d\omega
w = wi:dw:wf;                         % Frequency vector
alpha = 5;                            % Strength spectral density function J1(w)
s = 2.5;                              % Ohmic-parameter s
wc = 0.1;                             % Cut-off frequency
J0 = 0.2;                             % Strength spectral density function J2(w)
w0 = 2;                               % Resonant frequency
Gamma = 0.1;                          % Width of the spectral density function J2(w)
J1 = alpha*wc^(1-s)*w.^s.*exp(-w/wc); % Spectral density funcion J1(w)
J2 = J0*(Gamma/2)./((w-w0).^2+...
    (Gamma/2)^2).*w.^s./(w/w0+1).^2; % Spectral density funcion J2(w)
box on
hold on
plot(w,J1,'r-','Linewidth',2)
plot(w,J2,'b-','Linewidth',2)
hold off
xlabel('$\omega$','Interpreter','LaTex','Fontsize', 30)
ylabel('$J(\omega)$','Interpreter','LaTex','Fontsize', 30)
legend({'$J_1(\omega)$','$J_2(\omega)$'},'Interpreter','latex','Fontsize', 21,'Location','northeast')
set(gca,'fontsize',21)
xlim([0 5])
Nt = 5000;                 % Number of points for time
ti = 0;                    % Initial time
tf = 100;                  % Final time
dt = (tf-ti)/(Nt-1);       % Step time dt
t = ti:dt:tf; t= t';       % Time vector
Psi_0 = [1 1]'/sqrt(2);    % Initial wavefunction
rho1 = Psi_0*Psi_0';       % Initial density matrix for J_1(w)
rho2 = rho1;               % Initial density matrix for J_2(w)
p11 = zeros(size(t));      % Matrix element rho_11
p22 = zeros(size(t));      % Matrix element rho_22
p12 = zeros(size(t));      % Matrix element rho_12
p21 = zeros(size(t));      % Matrix element rho_21
wa = ones(size(t))*w;      % Auxiliar frecuency vector
J1a = ones(size(t))*J1;    % Auxiliar J1 vector
J2a = ones(size(t))*J2;    % Auxiliar J2 vector
ta = t*ones(size(w));      % Auxiliar time vector
T = 0.001*w0;              % Temperature
gamma1 = sum(J1a./wa.*sin(wa.*ta).*coth(wa/T/2),2)*dw;  % Rate gamma_1(t)
gamma_1_teo = alpha*wc*gamma(s)*sin(s*atan(wc*t))./(1+(wc*t).^2).^(s/2);
gamma2 = sum(J2a./wa.*sin(wa.*ta).*coth(wa/T/2),2)*dw;  % Rate gamma_2(t)
figure()
box on
hold on
plot(t,gamma1,'r-','Linewidth',2)
plot(t,gamma_1_teo,'k--','Linewidth',2)
plot(t,gamma2,'b-','Linewidth',2)
hold off
xlabel('$t$','Interpreter','LaTex','Fontsize', 30)
ylabel('$\gamma(t)$','Interpreter','LaTex','Fontsize', 30)
legend({'$\gamma_1(t)$','$\gamma_1^{\rm teo}(t)$','$\gamma_2(t)$'},'Interpreter','latex','Fontsize', 21,'Location','northeast')
set(gca,'fontsize',21)
xlim([0 100])
nk = 15;                  % nk steps per interval dt
C1 = zeros(size(t));
C2 = zeros(size(t));
N1 = zeros(size(t));
N2 = zeros(size(t));
sz = [1 0;                 % Pauli matrix s_z
    0 -1];
for n = 1:length(t)
    for k = 1:nk
        L1 = gamma1(n)*(sz*rho1*sz -rho1);      % Lindbladian for gamma_1
        L2 = gamma2(n)*(sz*rho2*sz -rho2);      % Lindbladian for gamma_2
        rho1_pred = rho1 + L1*dt/nk;            % Predictor \rho_1
        rho2_pred = rho2 + L2*dt/nk;            % Predictor \rho_2
        rho1_m = 0.5*(rho1+rho1_pred);          % Mean \rho_1
        rho2_m = 0.5*(rho1+rho2_pred);          % Mean \rho_1
        L1 = gamma1(n)*(sz*rho1_m*sz -rho1_m);  % Lindbladian using mean \rho_1
        L2 = gamma2(n)*(sz*rho2_m*sz -rho2_m);  % Lindbladian using mean \rho_2
        rho1 = rho1 + L1*dt/nk;                 % Density matrix \rho_1
        rho2 = rho2 + L2*dt/nk;                 % Density matrix \rho_2
    end
    C1(n) = 2*abs(rho1(1,2));                   % Coherence C1(t)
    C2(n) = 2*abs(rho2(1,2));                   % Coherence C2(t)
    f1 = (abs(gamma1(1:n))-gamma1(1:n));
    N1(n) = sum(f1)*dt;                         % Degree of non-Markovianity for gamma_1
    f2 = (abs(gamma2(1:n))-gamma2(1:n));
    N2(n) = sum(f2)*dt;                         % Degree of non-Markovianity for gamma_2
end
figure()
box on
hold on
plot(t,C1,'r-','Linewidth',2)
plot(t,C2,'b-','Linewidth',2)
hold off
xlabel('$t$','Interpreter','LaTex','Fontsize', 30)
ylabel('$C(t)$','Interpreter','LaTex','Fontsize', 30)
legend({'$C_1(t)$','$C_2(\omega)$'},'Interpreter','latex','Fontsize', 21,'Location','northeast')
set(gca,'fontsize',21)
xlim([0 50])

figure()
box on
hold on
plot(t,N1,'r-','Linewidth',2)
plot(t,N2,'b-','Linewidth',2)
hold off
xlabel('$t$','Interpreter','LaTex','Fontsize', 30)
ylabel('$\mathcal{N}_{\gamma}(t)$','Interpreter','LaTex','Fontsize', 30)
legend({'$\mathcal{N}_{\gamma_1}(t)$','$\mathcal{N}_{\gamma_2}(t)$'},'Interpreter','latex','Fontsize', 21,'Location','best')
set(gca,'fontsize',21)
xlim([0 100])
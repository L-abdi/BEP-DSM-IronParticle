function [t,x,Tp2,Tg,Cg,mdot_max,mdot_R,dT,Gd,Tgjsum,Tgj,Cgj,dmFeO_dt,Ap] = ParticleModel(rhoFe,rhoFeO,MW_N2,MW_O2,MW_Fe,MW_FeO, YO20, YN20 ...
,cp_g, rho_g,lambda_g,D,gamma_FeO_Fe, gamma_FeO_O2,hf_FeO,qFeO,k0FeO,TaFeO, ...
p0,Ru,k0,Ea,Sh,Nu,Sc,Pr, delta0, rp0, Tp0, Tg0,tEnd,h)

%% Forward Euler: yn+1 = yn + hf(yn)

% Initial Conditions
X0 = rp0*delta0; % Initial wustite layer thickness
    
rFe0 = rp0*(1.0-delta0); % Initial iron core radius
rFeO0 = rp0;             % Initial particle radius
mFe0 = 4/3*pi* (rFe0^3)*rhoFe; % Initial iron core mass
mFeO0 = 4/3*pi*(X0^3)*rhoFeO; % Initial wustite layer mass
ep0 = energyParticle(mFe0, mFeO0, Tp0); % Initial energy of particle
mFe_e = mFe0 * (0.01/100); % 0.01% of initial iron core mass
Cg0 = rho_g*YO20; %Initial oxidizer concentration in bulk gas

%% Governing ODEs
% ODE of the form dx = f(x), where dx1 = dmFeO/dt; dx2 = dmFe/dt; dx3 = dep/dt
f = cell(3,1); 
   
%% Forward Euler scheme

t = 0:h:tEnd;
N = length(t);

% Pre-allocating
mFe = zeros(1,N); mFeO = mFe; ep = mFeO; 
mdot_max = mFe; mdot_R = mdot_max; Tp2 = ep; Tg = Tp2; Ap=ep; dmFeO_dt = ep;
Tgj = Tg; Tgjsum = Tgj; Cgj = Tgj; Cgjsum = Tgj; Cg = Tg; Gd=Cg; dT= Tp2; 
xi = [0,0,0]; % Particle position
R = 3; % Dimension, for Green's function
% dx = zeros(length(xcoor)); % Distance from point to particle
dx = 0; % Equations solved locally, so distance to particle is zero

x = [mFe; mFeO; ep]; % ODE of the form dx = f(x)
x(1,1) = mFe0; %IC
x(2,1) = mFeO0; %IC
x(3,1) = ep0; %IC
Tg(1) = Tg0;
Cg(1) = Cg0;
Tp2(1) = Tp0;
mdot_R(1) = kineticRate(mFe0, mFeO0, ep0); % Kinetic rate
mdot_max(1) = 4 * pi * rFeO0 * D * Cg0; % Maximum diffusion rate

for n = 1:N

    mdot_max(n) = 4 * pi * radiusFeO(x(1,n), x(2,n), x(3,n))* D * Cg(n);
    if mdot_max(n) < kineticRate(x(1,n), x(2,n), x(3,n)) %Switch-type model
        dmFeO_dt(n) = mdot_max(n);
    else     
        dmFeO_dt(n) = kineticRate(x(1,n), x(2,n), x(3,n));
    end
    if x(1,n) <= mFe_e
        f{1} = @(mFe, mFeO, ep) 0;
        f{2} = @(mFe, mFeO, ep) 0;
        f{3} = @(mFe, mFeO, ep)  - areaParticle(mFe, mFeO, ep) .* hp(mFe, mFeO, ep) .* (temperatureParticle(mFe, mFeO, ep)-Tg(n));
    else
        f{1} = @(mFe, mFeO, ep) -dmFeO_dt(n) ./ gamma_FeO_Fe;
        f{2} = @(mFe, mFeO, ep) dmFeO_dt(n);
        f{3} = @(mFe, mFeO, ep) qFeO .* dmFeO_dt(n) - areaParticle(mFe, mFeO, ep) .* hp(mFe, mFeO, ep) .* (temperatureParticle(mFe, mFeO, ep)-Tg(n));
    end

    % Explicit forward Euler
    x(:,n+1) = x(:,n) + h .* [f{1}(x(1,n),x(2,n),x(3,n)); 
        f{2}(x(1,n),x(2,n),x(3,n)); f{3}(x(1,n),x(2,n),x(3,n))];
   
    Tp2(n) = real(temperatureParticle(x(1,n), x(2,n), x(3,n))); %Particle temperature
    h_Fe(n) = hFe(Tp2(n));
    h_FeO(n) = hFeO(Tp2(n));
    mdot_R(n) = kineticRate(x(1,n), x(2,n), x(3,n));
    Ap(n) = areaParticle(x(1,n), x(2,n), x(3,n));
     
    for i = 1:(n-1) % Approximating the convolution integral

        dt = h;
        if t(n-i) == 0
            Gd(n-i) = 0;
        else
            Gd(n-i) = (1./(4*pi*D*t(n-i)).^(R/2)).*...  % Green's function
                        exp(-abs(dx).^2./(4*D*t(n-i)));
        end

        dT(i) = Tp2(i)-Tg(n); % Difference in temperature, locally! 

        Tg(2) = Tg0 +  dt .* (1./(rho_g .* cp_g)).*(Gd(1) .* hp(x(1,1), x(2,1), x(3,1)) .* ...
            areaParticle(x(1,1), x(2,1), x(3,1)).* (dT(1)));
        Tgj(i+1) = Tgj(i) + (dt./ (rho_g .* cp_g) ) .* ( Gd(n-i) .* hp(x(1,i), x(2,i), x(3,i)) .* ...
            Ap(i).* dT(i) ) ;
        Tg(n+1) = Tg0 + Tgj(i+1);

 
        Cgj(i+1) = Cgj(i) + Gd(n-i) .* dmFeO_dt(i) .* dt;
        Cg(2) = Cg0 - dt* (Gd(1) .* mdot_R(1));
        Cg(n+1) = Cg0 - Cgj(i+1);

    end
end

return

% %% Visualisation
% close all
% 
% figure(1)
% subplot(2,2,1)
% plot(t,x(1,1:end-1))
% ylabel('m_{Fe} [kg]')
% subplot(2,2,2)
% plot(t,x(2,1:end-1))
% ylabel('m_{FeO} [kg]')
% subplot(2,2,3)
% plot(t,x(3,1:end-1))
% ylabel('e_p [J]')
% 
% figure(2)
% plot(t*1e3,Tp2)
% grid on
% xlabel('t [ms]')
% ylabel('T_p [K]')
% ylim([500 3500])
% xlim([0 50])
% hold on
% plot(t*1e3,Tg,'r-.')
% lgd = legend('T_p','T_g');
% lgd.Location = 'northwest';
% saveas(gcf,sprintf('TempHistories/Tp%g_Tg%g_COUPLED.png',Tp0,Tg0))

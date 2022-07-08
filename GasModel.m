%% Gasmodel
clc;
clear all;

tic
%% Input constant parameters %%

rhoFe = 7874.0;                        % Solid-phase iron(Fe) density, [kg/m^3]
rhoFeO  = 5745.0;                      % Solid-phase FeO density, [kg/m^3]
MW_N2 = 28.0134e-3;                    % Nitrogen molar weight, [kg/mol]
MW_O2 = 31.9988e-3;                    % Oxygen molar weight, [kg/mol]
MW_Fe = 55.845e-3;                     % Pure iron molar weight, [kg/mol]
MW_FeO = 71.844e-3;                    % FeO molar weight, [kg/mol]
YO20 = 0.21 ;                          % 02 fraction in air [-]                     
YN20 = 0.79 ;                          % N2 fraction in air [-]

% Fe + 1/2 O2 ---> FeO
n_O2 = 1/2 ;                           % The number of molecules of oxygen(O2), [-]
n_FeO = 1 ;                            % The number of moelcules of FeO, [-]
gamma_FeO_Fe = (n_FeO * MW_FeO/MW_Fe) ; % Fuel(Fe) - Oxidizer(O2) mass stoichiometric coefficient: [-]
gamma_FeO_O2 = (n_FeO/n_O2 * MW_FeO/MW_O2) ; % Product(FeO) - Oxidizer(O2) mass stoichiometric coefficient: [-]
hf_FeO = -272.04*1000 ;                % Standard enthalpy of Solid phase Wustite(FeO_s), [J/mol]
qFeO = -hf_FeO/MW_FeO ;                % Specific heat release(=heating value), [J/kg]
k0FeO = 2.6698e-4;                     % Parabolic growth rate [m2/s]
TaFeO = 20319;                         % Activation temperature [K]

p0 = 101325 ;                          % Standard Pressure, [Pa] (= 1[atm])
Ru = 8.3145 ;                          % Universal(Ideal) gas constant, [J/(mol*K)]
k0 = 5.34e-4 ;                         % (Calibrated) Parabolic Kinetic Rate for low temp oxidation process, [m^2/s]
Ea = 168.9*1000 ;                      % (Calibrated) Arnhenius Activation Energy for low temp oxidation process, [J/mol]
Sh = 2.0 ;                             % Sherwood number, [-], (the ratio of the convective mass transfer to diffusive mass transport)
Nu = 2.0 ;                             % Nusselt number, [-], (the ratio of convective to conductive heat transfer)
Sc = 1.0 ;                             % Schmidt number, mu/(rho*D), (the ratio of mass diffusivity)
Pr = 1.0 ;                             % Prantle number, cp*mu/lambda, (the ratio of thermal diffusivity)
Le = Sc/Pr ;                           % Lewis number, Sc/Pr, (the ratio of thermal diffusivity to mass diffusivity)
R_g = Ru * ( YO20 / MW_O2 + YN20 / MW_N2 ); % Specific gas constant

cp_g = 677;                            % Specific heat capacity gas [J/kgK]
rho_g = 1.225;                         % Density of gas [kg/m3]
lambda_O2 = 0.02659;                   % Thermal conductivity of O2 [W/mK]
lambda_N2 = 0.02614;                   % Thermal conductivity of N2 [W/mK]
lambda_g = 0.5*(YO20*lambda_O2+YN20*lambda_N2+1/(YO20/lambda_O2+YN20/lambda_N2)); % Thermal conductivity of gas [W/mK]
D = lambda_g / (cp_g * rho_g *Le);     % Diffusivity of the oxidizer [m^2/s]

%% Initial conditions %%

delta0 = 1e-3; % delta = rp/X
rp0 = 10.0e-6; % Initial particle radius
    
Tp0 = 1500; % Initial particle temperature 
Tg0 = 1150; % Initial gas temperature   

Cg0 = rho_g*YO20; % Initial oxidizer concentration in bulk gas

tEnd = 80e-3; % Time duration in seconds
h = 1e-4; % Time-step forward Euler

% Solving governing equations of the particle model
[t,x,Tp2,Tg,Cg,mdot_max,mdot_R,dT,Gd,Tgjsum,Tgj,Cgj,dmFeO_dt,Ap] = ParticleModel(rhoFe,rhoFeO,MW_N2,MW_O2,MW_Fe,MW_FeO, YO20, YN20 ...
,cp_g, rho_g,lambda_g,D,gamma_FeO_Fe, gamma_FeO_O2,hf_FeO,qFeO,k0FeO,TaFeO, ...
p0,Ru,k0,Ea,Sh,Nu,Sc,Pr, delta0, rp0, Tp0, Tg0,tEnd,h);
hp2 = hp(x(1,:), x(2,:), x(3,:));
Ap2 = areaParticle(x(1,:), x(2,:), x(3,:));

%% Bulk gas model %%
% Grids are solely defined for visualization of results!

% Pre-allocating
N = length(t);
M = N;
xi =[0,0]; % Particle position
R = 2; % Dimension
rEnd = 1e-3; % Right boundary of grid
r0 = xi(1);
dxx = (rEnd-r0)/(M-1); % Spatial step-size
xcoor=(r0:dxx:rEnd); 
ycoor=(r0:dxx:rEnd);
dt = h;
nxi = find(xcoor==xi(1)); 
dx2=sqrt((xcoor(:)-xi(1)).^2+(ycoor(:)-xi(2)).^2); % Distance from point in space to particle
Cgj = zeros(M,N);
Cg2 = zeros(M,N); % Gaseous oxidizer concentration
Gd = Cgj; Gh = Gd; % Greens function 
Tgj = Gd; Cgjsum = Cgj; Tgjsum = Tgj;
Tg2 = Cg2; % Bulk gas temperature
Tg2(:,1) = Tg0; % IC
Cg2(:,1) = Cg0; % IC

%% Solving gas equations

for n = 1:N % Approximating the convolution integral
    for i = 1:(n-1) % Time

        dt = h;
        if t(n-i) == 0
            Gd(:,n-i) = 0;
        else
            Gd(:,n-i) = (1./(4*pi*D.*t(n-i)).^(R/2)).*...    % Green's function
                        exp(-abs(dx2(:)).^2./(4.*D.*t(n-i)));
        end

        dT(i) = Tp2(i)-Tg2(1,n); % Difference in temperature, locally! 
        
        Tg2(:,2) = Tg0 +  dt .*  Gd(:,1) .* hp2(1) .* ...
            Ap(1) .* (dT(1));

        Tgj(:,i+1) = Tgj(:,i) + Gd(:,n-i) .* hp2(i) .* ...
            Ap(i).* dT(i) .* dt ;

        Tg2(:,n+1) = Tg0 + Tgj(:,i+1);
        

        Cg2(:,2) = Cg0 - dt .* (Gd(:,1) .* dmFeO_dt(1));

        Cgj(:,i+1) = Cgj(:,i) + Gd(:,n-i) .* dmFeO_dt(i) .* dt;
 
        Cg2(:,n+1) = Cg0 - Cgj(:,i+1);
    end
end

Tg2 = real(Tg2);
Cg2 = real(Cg2);
r = dx2;

%% Visualization
close all

figure(12)
plot(t(3:end)*1e3, real(mdot_max(3:end))*1e9)
hold on
plot(t(3:end)*1e3, real(mdot_R(3:end))*1e9,'r-.')
ylim([0 10])
legend('\textbf{$\dot{m}_{D,max}$}','\textbf{$\dot{m}_R$}', 'Interpreter','latex')
% title('Reaction rates over time')
xlabel('t [ms]')
ylabel('\textbf{$\dot{m}_{D,max}$} ,\textbf{$\dot{m}_R$  $\cdot 10^9$ [kg/s]}', 'Interpreter','latex')
txt = {sprintf('T_{p,0} = %g K',Tp0) ,sprintf('T_{g,0} = %g K',Tg0) ,...
    sprintf('\\delta_0 = %g m',delta0)};
text(78, 5, txt, 'Horiz','right', 'Vert','bottom')
grid on
set(gcf,'position',[0 0 1000 300])
% % saveas(gcf,sprintf('ReactionRates/RR_Tp%g_Tg%g.png',Tp0,Tg0))


figure(11)
plot(t*1e3,Tp2)
hold on
plot(t*1e3,real(Tg2(1,1:end-1)),'g-.')
grid on
xlabel('t [ms]')
ylabel('T [K]')
ylim([1100 2200])
legend('T_p','T_g')

figure(23)
subplot(2,1,1)
plot(t*1e3,real(Tg2(1,1:end-1)))
grid on
xlabel('t [ms]')
ylabel('T [K]')
ylim([1150 1250])
subplot(2,1,2)
plot(t*1e3,Cg2(1,1:end-1))
grid on
ylim([min(min(Cg2)) Cg0])
xlabel('t [ms]')
ylabel('C_{O_2} [kg/m^3]')


toc
return










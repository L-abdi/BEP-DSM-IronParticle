%% 

clear all
clc

%% Input constant parameters
rhoFe = 7874.0;                                     % Solid-phase iron(Fe) density, [kg/m^3]
rhoFeO  = 5745.0;                                     % Solid-phase FeO density, [kg/m^3]
MW_N2 = 28.0134e-3;                                    % Nitrogen molar weight, [kg/mol]
MW_O2 = 31.9988e-3;                                       % Oxygen molar weight, [kg/mol]
MW_Fe = 55.845e-3;                                      % Pure iron molar weight, [kg/mol]
MW_FeO = 71.844e-3;                                      % FeO molar weight, [kg/mol]
% MW_Fe3O4 = 231.533e-3;                                      % Fe3O4 molar weight, [kg/mol]
YO20 = 0.21 ;
YN20 = 0.79 ;
cp_g = 677; %[J/kgK]
rho_g = 1.225; %[kg/m3]
%NASA polynomials
lambda_O2 = 0.02659; %[W/mK]
lambda_N2 = 0.02614; %[W/mK]
lambda_g = 0.5*(YO20*lambda_O2+YN20*lambda_N2+1/(YO20/lambda_O2+YN20/lambda_N2));

%FeO wustite
%Fe3O4 magnetite
%Fe2O3 hematite


% Fe + 1/2 O2 ---> FeO
n_O2 = 1/2 ;                                           % The number of molecules of oxygen(O2), [-]
n_FeO = 1 ;                                            % The number of moelcules of FeO, [-]
gamma_FeO_Fe = (n_FeO * MW_FeO/MW_Fe) ;                % Fuel(Fe) - Oxidizer(O2) mass stoichiometric coefficient: [-]
gamma_FeO_O2 = (n_FeO/n_O2 * MW_FeO/MW_O2) ;           % Product(FeO) - Oxidizer(O2) mass stoichiometric coefficient: [-]
hf_FeO = -272.04*1000 ;                                % Standard enthalpy of Solid phase Wustite(FeO_s), [J/mol]
qFeO = -hf_FeO/MW_FeO ;                                % Specific heat release(=heating value), [J/kg]
k0FeO = 2.6698e-4;                                     % Parabolic growth rate [m2/s]
TaFeO = 20319;                                         % [K]

p0 = 101325 ;                                        % Standard Pressure, [Pa] (= 1[atm])
Ru = 8.3145 ;                                        % Universal(Ideal) gas constant, [J/(mol*K)]
k0 = 5.34e-4 ;                                       % (Calibrated) Parabolic Kinetic Rate for low temp oxidation process, [m^2/s]
Ea = 168.9*1000 ;                                    % (Calibrated) Arnhenius Activation Energy for low temp oxidation process, [J/mol]
Sh = 2.0 ;                                           % Sherwood number, [-], (the ratio of the convective mass transfer to diffusive mass transport)
Nu = 2.0 ;                                           % Nusselt number, [-], (the ratio of convective to conductive heat transfer)
Sc = 1.0 ;                                           % Schmidt number, mu/(rho*D), (the ratio of mass diffusivity)
Pr = 1.0 ;                                           % Prantle number, cp*mu/lambda, (the ratio of thermal diffusivity)
Le = Sc/Pr ;                                         % Lewis number, Sc/Pr, (the ratio of thermal diffusivity to mass diffusivity)
R_g = Ru * ( YO20 / MW_O2 + YN20 / MW_N2 );          % Specific gas constant
D = lambda_g / (cp_g * rho_g *Le); %Diffusivity of the oxidizer
%% Initial conditions

delta0 = 1e-3;
rp0 = 10.0e-6;
    
Tp0 = 1500;
Tg0 = 1150;   

Cg0 = rho_g*YO20; %Initial oxidizer concentration in bulk gas
%% Forward Euler: yn+1 = yn + hf(yn)
                    
Le = Sc/Pr ;                                         % Lewis number, Sc/Pr, (the ratio of thermal diffusivity to mass diffusivity)
R_g = Ru * ( YO20 / MW_O2 + YN20 / MW_N2 );          % Specific gas constant

% Initial Conditions
X0 = rp0*delta0;
    
rFe0 = rp0*(1.0-delta0);
rFeO0 = rp0;
mFe0 = 4/3*pi* (rFe0^3)*rhoFe;
mFeO0 = 4/3*pi*(X0^3)*rhoFeO;
ep0 = energyParticle(mFe0, mFeO0, Tp0);

mFe_e = mFe0 * (0.01/100); % 0.01% of initial iron core mass

%% Governing ODEs
% ODE of the form dx = f(x), where dx1 = dmFeO/dt; dx2 = dmFe/dt; dx3 = dep/dt

f = cell(3,1); 
   

%% Forward Euler scheme

r0 = 0;
rEnd = 5e-3;

% N = 8e2;
M = 6000
tEnd = 55e-3;
% dr = 1e-5;
% dt = tEnd/(N-1);
% dt = 2e-12;
% dt = tEnd/(N-1)
dr = (rEnd-r0)/(M-1)
rn = 5e-6;
dt = (0.40*dr^2)/D

r = r0:dr:rEnd;
t = 0:dt:tEnd;
s1 = dt/dr^2;
s2 = dt/(2*dr);
% D*s1
N = length(t)


mFe = zeros(1,N); mFeO = mFe; ep = mFeO; 
mdot_max = mFe; mdot_R = mdot_max; Tp2 = ep; Ap = ep; Tg = Tp2; Tf = Tp2;
Tgj = Tg; Tgjsum = Tgj; Cgj = Tgj; Cgjsum = Tgj; Cg = Tg; Gd=Cg; dT= Tp2;%pre-allocating
xi =[0,0,0]; %Particle position
R = 3; %Dimension
dxx = 1e-3/(N-1);
xcoor=(0:dxx:1e-3)+xi(1);
ycoor=(0:dxx:1e-3)+xi(2);
% dx = zeros(length(xcoor)); %Distance from point to particle


x = [mFe; mFeO; ep]; %ODE of the form dx = f(x)
x(1,1) = mFe0; %IC
x(2,1) = mFeO0; %IC
x(3,1) = ep0; %IC
Tg(1) = Tg0;
Tw = Tg0;
Cw = Cg0;
Cg(1) = Cg0;
Tp2(1) = Tp0;
mdot_R(1) = kineticRate(mFe0, mFeO0, ep0);
mdot_max(1) = 4 * pi * rFeO0 * D * Cg0;

C = zeros(M,N);
C(:,1) = Cg0; %IC
C(:,2) = Cg0; %IC
C(end,:) = Cw; %BC
C(end-1,:) = Cw; %BC

T = zeros(M,N);
T(:,1) = Tg0; %IC
T(:,2) = Tg0; %IC
T(end,:) = Tw; %BC
T(end-1,:) = Tw; %BC

a = lambda_g/(rho_g*cp_g);
%%

for n = 1:N
    h = dt;
    mdot_max(n) = 4 * pi * radiusFeO(x(1,n), x(2,n), x(3,n))* D * C(3,n);
    if mdot_max(n) < kineticRate(x(1,n), x(2,n), x(3,n))
        dmFeO_dt(n) = mdot_max(n);
    else     
        dmFeO_dt(n) = kineticRate(x(1,n), x(2,n), x(3,n));
    end
    if x(1,n) <= mFe_e
        f{1} = @(mFe, mFeO, ep) 0;
        f{2} = @(mFe, mFeO, ep) 0;
        f{3} = @(mFe, mFeO, ep)  - areaParticle(mFe, mFeO, ep) .* hp(mFe, mFeO, ep) .* (temperatureParticle(mFe, mFeO, ep)-T(2,n));
    else
        f{1} = @(mFe, mFeO, ep) -dmFeO_dt(n) ./ gamma_FeO_Fe;
        f{2} = @(mFe, mFeO, ep) dmFeO_dt(n);
        f{3} = @(mFe, mFeO, ep) qFeO .* dmFeO_dt(n) - areaParticle(mFe, mFeO, ep) .* hp(mFe, mFeO, ep) .* (temperatureParticle(mFe, mFeO, ep)-T(2,n));
    end   
    x(:,n+1) = x(:,n) + h .* [f{1}(x(1,n),x(2,n),x(3,n)); 
        f{2}(x(1,n),x(2,n),x(3,n)); f{3}(x(1,n),x(2,n),x(3,n))];
   
    Tp2(n) = temperatureParticle(x(1,n), x(2,n), x(3,n)); %Particle temperature
    mdot_R(n) = kineticRate(x(1,n), x(2,n), x(3,n));
    Ap(n) = areaParticle(x(1,n), x(2,n), x(3,n));
    
    %FDM
    for i = 2:M-1 %position
        
        T(1,n) = Tp2(n); %
        C(end,n) = Cg0;

        if abs(r(i)) < rn
            np(i) = 1/( (4/3)*pi*rn^3 );
        else
            np(i) = 0;
        end

        C(i,n+1) = C(i,n) + D*((s1 - (2*s2/r(i)))*C(i-1,n) - 2*s1*C(i,n) ... 
             + (s1 + (2*s2/r(i)))*C(i+1,n)) - np(i)*dmFeO_dt(n)*dt;

        dT(n) = Tp2(n)-T(2,n);

        Tf(n+1) = Tf(n) + (dt/(cp_g*rho_g)) * np(i) * hp(x(1,n), x(2,n), x(3,n)) * Ap(n) * dT(n);

        T(i,n+1) =  T(i,n) + a*( (s1 - (2*s2/r(i)))*T(i-1,n) - 2*s1*T(i,n) ... 
             + (s1 + (2*s2/r(i)))*T(i+1,n) ) + Tf(n+1);

    end 

end

Tp2 = real(Tp2);
C = real(C);
T = real(T);
%% Visualization

figure(12)
plot(t(2:end)*1e3, real(mdot_max(2:end))*1e9)
hold on
plot(t(2:end)*1e3, real(mdot_R(2:end))*1e9,'r-.')
ylim([0 10])
legend('\textbf{$\dot{m}_{D,max}$}','\textbf{$\dot{m}_R$}', 'Interpreter','latex')
% title('Reaction rates over time')
xlabel('t [ms]')
ylabel('\textbf{$\dot{m}_{D,max}$} ,\textbf{$\dot{m}_R$  $\cdot 10^9$ [kg/s]}', 'Interpreter','latex')
txt = {sprintf('T_{p,0} = %g K',Tp0) ,sprintf('T_{g,0} = %g K',Tg0) ,...
    sprintf('\\delta_0 = %g m',delta0)};
text(99, 5, txt, 'Horiz','right', 'Vert','bottom')
grid on
set(gcf,'position',[0 0 1000 300])

% % figure(1)
% % clf
% % plot(t,Tp2)
% % % hold on
% % % plot(t,T(2,1:end-1))
% 
% figure(2)
% % contourf(real(T),200,'linecolor','none')
% surf(t,r,T(:,1:end-1))
% title('FDM')
% xlabel('t')
% ylabel('r')
% cb = colorbar;
% ylim([0 1e-3])
% % zlim([1150 1300])
% % caxis([Tg0 max(max(T))+100])
% colormap('hot')  
% shading interp
% view(0,90)


return





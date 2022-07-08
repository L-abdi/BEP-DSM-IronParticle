%% Solid iron enthaply (J/mol)
function cp_FeO = cpFeO(Tp)
%% Coefficients for [cp & del_h] of Solid Iron Oxide(FeO)
AA0 = 45.75120 ;
AA1 = 18.78553 ;
AA2 = -5.952201 ;
AA3 = 0.852779 ;
AA4 = -0.081265 ;
AA5 = -286.7429 ;
AA6 = -272.0441 ;

%% Calculate cp [J/mol/K]
% if (Tp <= 1650)
    cp_FeO = AA0 + AA1*(Tp/1000) ...              % [J/(mol*K)]
        + AA2*(Tp/1000)^2 ...
        + AA3*(Tp/1000)^3 ...
        + AA4/(Tp/1000)^2 ;
% end

return
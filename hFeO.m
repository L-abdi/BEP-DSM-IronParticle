%% Solid FeO enthaply (J/mol)
function h_FeO = hFeO(Tp)
%% Coefficients for [cp & del_h] of Solid Iron Oxide(FeO)
AA0 = 45.75120 ;
AA1 = 18.78553 ;
AA2 = -5.952201 ;
AA3 = 0.852779 ;
AA4 = -0.081265 ;
AA5 = -286.7429 ;
AA6 = -272.0441 ;

%% Calculate enthaply per mol
% if (Tp <= 1650)
    h_FeO = ( AA0*(Tp/1000) ...                   % [J/mol]
        + AA1*((Tp/1000)^2)/2 ...
        + AA2*((Tp/1000)^3)/3 ...
        + AA3*((Tp/1000)^4)/4 ...
        - AA4/(Tp/1000) ...
        + AA5 - AA6 )*1000 ;
% end
return
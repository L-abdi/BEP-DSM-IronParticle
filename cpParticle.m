%% Calculate iron particle's internal energy
function cp = cpParticle(mFe, mFeO, Tp)
MW_Fe = 55.845e-3;                                      % Pure iron molar weight, [kg/mol]
MW_FeO = 71.844e-3;                                     % FeO molar weight, [kg/mol]
cp = mFe / MW_Fe * cpFe(Tp) + mFeO / MW_FeO * cpFeO(Tp);
return
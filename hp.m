function hp2 = hp(mFe,mFeO,ep)
Nu = 2;
lambda_g = 0.026233866029137;

hp2 = Nu * lambda_g ./ (2 * radiusFeO(mFe,mFeO,ep)); % Convective heat transfer coefficient

return
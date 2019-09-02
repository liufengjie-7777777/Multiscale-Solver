clear all; clc

syms lr lt lz ufs positive;
syms thetaSMC gammar deltaM AS2 NCF NCU EAMp EAM positive;
syms nAMp nAM LMr Lfor LMz Lfoz I4SMCe positive;
syms eS2xr eS2yr eS2xz eS2yz positive;

MSMC = [sin(thetaSMC), cos(thetaSMC), 0];
MSMCr = [cos(thetaSMC), -sin(thetaSMC), 0];
MSMCz = [0, 0, 1];

%Deformation Gradients
F = diag([lr lt lz]);
Ffs = ufs*(MSMC'*MSMC)+(eye(3)-(MSMC'*MSMC))/sqrt(ufs);
Fe = simplify(F/Ffs); %=F*inv(Ffs)

%Right Cauchy-Green Tensors
C = transpose(F)*F;
Ce = simplify(transpose(Fe)*Fe);

%I4
I4SMCr = simplify(MSMCr*C*transpose(MSMCr));
I4SMCz = simplify(MSMCz*C*transpose(MSMCz));
%I4SMCe = simplify(MSMC*Ce*transpose(MSMC));

%sSMC
psiSMC = ( (nAMp*EAMp+nAM*EAM)/(4*deltaM) )*AS2*NCF * (gammar*Lfor*LMr+(1-gammar)*Lfoz*LMz) * power(sqrt(I4SMCe)-1,2);
dpsiSMC = simplify( diff(psiSMC,I4SMCe) );
sSMC = simplify( 2*Fe*(dpsiSMC*(MSMC'*MSMC))*transpose(Fe) );

PCU = simplify( MSMC* (2*Fe*(dpsiSMC*(MSMC'*MSMC)/Ffs'))* MSMC' );

%PMMy & SMMy
C1 = ((nAMp*EAMp+nAM*EAM)/deltaM)*AS2*NCU;
Cr = gammar*Lfor*LMr*eS2yr*(MSMCr'*MSMCr);
Cz = (1-gammar)*Lfoz*LMz*eS2yz*(MSMCz'*MSMCz);
PMMy = simplify( C1*(Cr+Cz) );

SMMy = simplify( PMMy*transpose(F) );

%sECM
syms Cp p a1 a4 positive;
c = sym('c',[2 4]);
alphaj = [a1 -a1 -a4 a4];

sECM = Cp/2*eye(3);
for j=1:1 %length(alphaj)
    MECMj = [0, cos(alphaj(j)), sin(alphaj(j))];
    I4fj = MECMj*(C*transpose(MECMj));
    sECM = sECM + 1*c(1,j)*(I4fj-1)*exp(c(2,j)*power(I4fj-1,2))*(MECMj'*MECMj);
end
sECM = simplify( 2*F*sECM*transpose(F) );

%Pisom Calc
Pisom = simplify( 2*Fe*(dpsiSMC*(MSMC'*MSMC)/Ffs') + PMMy );
Pisom = simplify(Pisom(2,2));


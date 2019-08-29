close all; clear all; clc

%Media
cM = 3*1e-3; %kPa to MPa
k1M = 2.3632*1e-3; %kPa to MPa
k2M = 0.8393;

HM = 0.26; %mm
betaM = 29*pi/180; %deg to rad

%Adventitia
cA = 0.3*1e-3; %kPa to MPa
k1A = 0.5620*1e-3; %kPa to MPa
k2A = 0.7112;

HA = 0.13; %mm
betaA = 62*pi/180; %deg to rad

%Reference Conf.
Ri = 0.71; alpha = 0*pi/180; %mm ; deg to rad
%Ri = 1.43; alpha = 160*pi/180; %mm ; deg to rad

Rm = Ri + HM;
Ro = Rm + HA;

k = 2*pi/(2*pi-alpha);
r = @(R,Ri,ri,lz) sqrt( (R.^2-Ri.^2)./(k.*lz)+ri.^2 );
R = @(r,Ri,ri,lz) sqrt( (r.^2-ri.^2).*(k.*lz)+Ri.^2 );

lr = @(R,Ri,ri,lz) R./(r(R,Ri,ri,lz)*k*lz);
lt = @(R,Ri,ri,lz) k.*r(R,Ri,ri,lz)./R;

syms ri;

Pin = linspace(0,25,50)*1e-3; %kPa to MPa
lz = 1.5:0.1:1.9;

for k=1:length(lz)
    strLegend{k} = num2str(lz(k));
    for i=1:length(Pin)
        rm = r(Rm,Ri,ri,lz(k));
        ro = r(Ro,Rm,rm,lz(k));

        [xM,wM] = lgwt(3,ri,rm);
        [xA,wA] = lgwt(3,rm,ro);

        for j=1:length(xM)
            Rj = R(xM(j),Ri,ri,lz(k));
            sECMM(:,j) = sECMCalc([-betaM,betaM],k1M,k2M,cM,lr(Rj,Ri,ri,lz(k)),lt(Rj,Ri,ri,lz(k)),lz(k));
        end
        for j=1:length(xA)
            Rj = R(xA(j),Rm,rm,lz(k));
            sECMA(:,j) = sECMCalc([-betaA,betaA],k1A,k2A,cA,lr(Rj,Rm,rm,lz(k)),lt(Rj,Rm,rm,lz(k)),lz(k));
        end

        S = sum( (sECMM(2,:)-sECMM(1,:))./xM.*wM ) - sum( (sECMA(2,:)-sECMA(1,:))./xA.*wA );

        riNew = vpasolve(S == Pin(i),ri,[0 2]);

        if ~isempty(riNew)
            riVec(i,k) = riNew;
            FVec(i,k) = vpa(subs(pi*sum( (2*sECMM(3,:)-sECMM(2,:)-sECMM(1,:)).*xM.*wM ) + sum( (2*sECMA(3,:)-sECMA(2,:)-sECMA(1,:)).*xA.*wA ),ri,riNew));
            fprintf('(%d-%d) ri=%.3f mm, F_T=%.3f N\n',k,i,riNew,FVec(i,k));
        end
    end
end

figure();

subplot(1,2,1);
plot(riVec,Pin*1e3);
ylabel('Pin (kPa)'); xlabel('ri (mm)');
legend(strLegend);
ylim([0 25]); xlim([0.3 1.5]);

subplot(1,2,2);
plot(riVec,FVec);
ylabel('F_T (N)'); xlabel('ri (mm)');
legend(strLegend);
xlim([0.3 1.5]);
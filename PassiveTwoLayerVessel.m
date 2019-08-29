clear variables; clc;

%Media
cM = 3e-3; %kPa to MPa
k1M = 2.3632e-3; %kPa to MPa
k2M = 0.8393;

HM = 0.26; %mm
betaM = 29*pi/180;

a01M = [0,cos(betaM), sin(betaM)]; a02M = [0, cos(betaM), -sin(betaM)];

%Adventitia
cA = 0.3e-3; %kPa to MPa
k1A = 0.5620e-3; %kPa to MPa
k2A = 0.7112;

HA = 0.13; %mm
betaA = 62*pi/180;

Ri = 1.43; %mm %0.71 for alpha=0
RoM = Ri+HM;
RoA = Ri+HM+HA;

a01A = [0,cos(betaA), sin(betaA)]; a02A = [0, cos(betaA), -sin(betaA)];

%Additional Parameters
k = 2*pi/(2*pi-160*pi/180); %

syms r ri positive;

syms lz positive;

roM = sqrt( ((RoM).^2 - Ri.^2)./(k*lz) + ri.^2 );
roA = sqrt( ((RoA).^2 - (RoM).^2)./(k*lz) + roM.^2 );

RM = sqrt( (r.^2 - ri.^2).*(k*lz) + Ri.^2 );
RA = sqrt( (r.^2 - roM.^2).*(k*lz) + RoM.^2 );

lrM = RM./(r.*k.*lz);
ltM = k.*r./(RM);

lrA = RA./(r.*k.*lz);
ltA = k.*r./(RA);

FM = diag([lrM ltM lz]);
CM = FM'*FM;

FA = diag([lrA ltA lz]);
CA = FA'*FA;

IiJ = @(a0iJ,C) a0iJ*C*a0iJ'; %trace((a0iJ'*a0iJ)*C); 

I4M = IiJ(a01M,CM); I6M = IiJ(a02M,CM);
I4A = IiJ(a01A,CA); I6A = IiJ(a02A,CA);

CauchyPassiveM = 2*FM*( cM/2*eye(3) +...
    k1M.*(I4M-1).*exp(k2M*(I4M-1).^2)*(a01M'*a01M)+k1M.*(I6M-1).*exp(k2M*(I6M-1).^2)*(a02M'*a02M) )*FM';
CauchyPassiveA = 2*FA*( cA/2*eye(3) +...
    k1A.*(I4A-1).*exp(k2A*(I4A-1).^2)*(a01A'*a01A)+k1A.*(I6A-1).*exp(k2A*(I6A-1).^2)*(a02A'*a02A) )*FA';

[xM,wM] = lgwt(3,ri,roM);
[xA,wA] = lgwt(3,roM,roA);

SrM = subs(CauchyPassiveM(1,1),r,xM);
SrA = subs(CauchyPassiveA(1,1),r,xA);

StM = subs(CauchyPassiveM(2,2),r,xM);
StA = subs(CauchyPassiveA(2,2),r,xA);

SzM = subs(CauchyPassiveM(3,3),r,xM);
SzA = subs(CauchyPassiveA(3,3),r,xA);

sM = sum(((StM-SrM)./xM).*wM);
sA = sum(((StA-SrA)./xA).*wA);

lzVal = linspace(1.5,1.9,5);
Pin = linspace(0,23*1e-3,20); %MPa
for j=1:length(lzVal)
    for i=1:length(Pin)
        y(i,j) = vpasolve(Pin(i) == subs(sM+sA,lz,lzVal(j)),[0 3]);
    end
end
figure(1); hold on;
for j=1:length(lzVal)
    plot(y(:,j),Pin*1e3);
end
hold off;
ylim([0 27]); xlim([0.3 1.5]); grid on;
ylabel('Pin (kPa)'); xlabel('ri (mm)');
legend('lz=1.5','','','','lz=1.9');

FM = subs(pi*sum(((2*SzM-StM-SrM).*xM).*wM),lz,lzVal);
FA = subs(pi*sum(((2*SzA-StA-SrA).*xA).*wA),lz,lzVal);

figure(2); hold on;
for j=1:length(lzVal)
    F(:,j) = double(subs(FA(j)+FM(j),ri,y(:,j)));
    plot(y(:,j),F(:,j));
end
hold off;
ylim([-0.06 0.3]); xlim([0.3 1.5]); grid on;
ylabel('F (N)'); xlabel('ri (mm)');
legend('lz=1.5','','','','lz=1.9');

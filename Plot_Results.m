%Plot Results

time = a.V.time;

ufs = a.V.ufs;
ri = a.V.ri;

b = ArteryVessel;

SamplePoints = length(a.V.time);
for i=1:SamplePoints
    b.cs.riNum = ri(i);
    b.xwNum;
    b.cs.r = b.cs.xNum;
    b.cs.ri = b.cs.riNum;
    b.cs.R = b.R;
    
    b.cs.lr = b.lr;
    b.cs.lt = b.lt;
    b.cs.lz = b.lz;
end

figure(1);
plot(obj.V.time,obj.V.Do);
grid on; ylim([600 1300]); xlim([0 obj.TotalTime]); %ylim([0 (ceil(max(obj.Do))+100)]);
ylabel('Do (um)'); xlabel('time (min)');
title(['lz=' num2str(obj.lz) ' Pin=' num2str(obj.cs.Pin/133.322387415*1e6) ' mmHg']);

figure(2);
plot(obj.V.time,obj.V.FT);
minFT = min(obj.V.FT);
maxFT = max(obj.V.FT);
if minFT>2
    minFT = ceil(minFT-2);
else
   minFT = 0;
end
if maxFT>0
    maxFT = ceil(maxFT+2);
else
    maxFT = 0;
end
grid on; ylim([minFT maxFT]); xlim([0 obj.TotalTime]);
ylabel('F_T (mN)'); xlabel('time (min)');

figure(3);
plot(obj.V.time,obj.V.stretch);
grid on; ylim([0 2]); xlim([0 obj.TotalTime]);
ylabel('Stretch Ratio'); xlabel('time (min)');
legend('\lambda_r','\lambda_\theta','\lambda_z','det(F)');

figure(4);
plot(obj.V.time,obj.V.PMMCU*1e3);
ylabel('Stress (kPa)'); ylim([0 ceil(max(obj.V.PMMCU*1e3,[],'all')+10)]);
hold on; yyaxis right;
plot(obj.V.time,obj.V.ufs);
h = legend('$P_{MM}$','$P_{CU}$','$\bar{u}_{fs}$ (right axis)');
ylabel('ufs'); xlabel('time (min)'); grid on; hold off;
set(h,'Interpreter','latex','fontsize',12);   

figure(5);
plot(obj.V.time, (obj.V.sECM+obj.V.sSMC+obj.V.sMMy)*1e3);
h = legend('$\bar{\sigma}_{r}$','$\bar{\sigma}_{\theta}$','$\bar{\sigma}_{z}$');
ylabel('Cauchy Stress (kPa)'); xlabel('time (min)'); grid on;
set(h,'Interpreter','latex','fontsize',12);

N=0;
if N
    figure();
    plot(a.V.time,a.V.dMAiMA0); title('dMAi/dMA0 vs time');
    legend('dMAr/dMA0','dMAz/dMA0');
    ylabel('dMAi/dMA0'); xlabel('time (min)');
    
    figure();
    plot(a.V.time,a.V.LMi*1e3);  title('L_{M,i} vs time');
    legend('LMr','LMz');
    ylabel('LMi (\mum)'); xlabel('time (min)');
    
    figure();
    plot(a.V.time,a.V.Lfoi); title('L_{fo,i} vs time');
    legend('L_{fo,r}','L_{fo,z}');
    ylabel('L_{fo,i}'); xlabel('time (min)');    
    
    figure();
    plot(a.V.time,a.V.eS2); title('\epsilon_{S2,i} vs time');
    legend('\epsilon_{S2xr}','\epsilon_{S2yr}','\epsilon_{S2xz}','\epsilon_{S2yz}');
    ylabel('\epsilon_{S2,i}'); xlabel('time (min)');    
    %ylim([0.06 0.15]);
    figure();
    plot(a.V.time,[a.V.x2c a.V.y2c]*1e3); title('x2c,i and y2c,i vs time');
    legend('x2cr','x2cz','y2cr','y2cz');
    ylabel('x2c & y2c (\mum)'); xlabel('time (min)');
    ylim([0.03 0.065]);
    
    figure();
    plot(a.V.time,[a.V.sECM(:,1) a.V.sSMC(:,1) a.V.sMMy(:,1)]*1e3);
    legend('\sigma_{ECMr}','\sigma_{SMCr}','\sigma_{MMyr}');
    ylabel('Cauchy Stress (kPa)'); xlabel('time (min)');
    title('Cauchy Stress Components in r Direction vs Time');
    
    figure();
    plot(a.V.time,[a.V.sECM(:,2) a.V.sSMC(:,2) a.V.sMMy(:,2)]*1e3);
    legend('\sigma_{ECM\theta}','\sigma_{SMC\theta}','\sigma_{MMy\theta}');
    ylabel('Cauchy Stress (kPa)'); xlabel('time (min)');
    title('Cauchy Stress Components in \theta Direction vs Time');
    
    figure();
    plot(a.V.time,[a.V.sECM(:,3) a.V.sSMC(:,3) a.V.sMMy(:,3)]*1e3);
    legend('\sigma_{ECMz}','\sigma_{SMCz}','\sigma_{MMyz}');
    ylabel('Cauchy Stress (kPa)'); xlabel('time (min)');
    title('Cauchy Stress Components in z Direction vs Time');
end
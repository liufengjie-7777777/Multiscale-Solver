clear all; close all; clc; %#ok<CLALL>

%Define a as artrey and calc n
if input('Which experiment do you want to simulate?\nEnter 1 for biaxial and 0 for uniaxial: ')
    a = ArteryVessel;
    a.cs.Pin = 90*133.322387415*1e-6; %Inner pressure
    a.cs.lambda = 1.5; %Axial loaded stretch ratio 
else
    a = ArteryStrip;
    a.cs.ltG = 1.69; %Circumferential stretch ratio
end

a.LMmax = 0.5e-3;

if ~a.InitialParameters
    SamplePoints = length(a.V.time);
    %Active Simulation
    for i=1:SamplePoints
        fprintf('(%d/%d) Time=%.1f min: ',i,SamplePoints,a.V.time(i));
        tic
        if a.stepCalc(i)
            fprintf('Error calculating current step\n');
            i = SamplePoints+1; %#ok<FXSET>
        else
            eTime = toc;
            fprintf('(Step time is %.2f sec. Elapsed time is %.1f min)\n',eTime,(SamplePoints-i)*eTime/60);
        end
    end
    if i==SamplePoints
        a.PlotResults;
    end
end

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
    
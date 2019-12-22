clear all; close all; clc; %#ok<CLALL>

%Define a as artrey and calc n
if input('Which experiment do you want to simulate?\nEnter 1 for biaxial and 0 for uniaxial: ')
    a = ArteryVessel;
    a.cs.Pin = 90*133.322387415*1e-6; %Inner pressure
    a.cs.lambda = 1.5; %Axial loaded stretch ratio
    
    Pin = linspace(0,90,10);
for i=1:length(Pin)
    a.cs.Pin = Pin(i)*133.322387415e-6;
    a.InitialParameters;
end

else
    a = ArteryStrip;
    a.cs.ltG = 1.69; %Circumferential stretch ratio
end

% a.LMmax = a.LMmax*0.7; %um to mm - same  units as delta_m
% a.EAMp = a.EAMp*1.43; %GPa to MPa
a.alphaPS = 30*pi/180;

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


if 0
save('SimulationPi90lz15-LMmax and EAMp decrease.mat','a')

clear all; close all;

load('SimulationPi90lz15.mat')
plot(a.V.time,a.V.Do);
hold on
load('SimulationPi90lz15-LMmax decrease.mat')
plot(a.V.time,a.V.Do);
load('SimulationPi90lz15-LMmax and EAMp decrease.mat')
plot(a.V.time(1:150:end),a.V.Do(1:150:end),'xr');

xlabel('Time (min)'); ylabel('Do (\mum)');
legend('Original values','0.7\timesL_{M,max}','0.7\timesL_{M,max} & 1.43\timesE_{AMp}');

hold off
end

clear all; close all; clc; %#ok<CLALL>

%Define a as artrey and calc n
if input('Which experiment do you want to simulate?\nEnter 1 for biaxial and 0 for uniaxial: ')
    a = ArteryVessel;
    a.cs.lambda = 1.5; %Axial loaded stretch ratio
    
    Pin = linspace(0,90,10); %Inner pressure (mmHg)
    for i=1:length(Pin)
        a.cs.Pin = Pin(i)*133.322387415e-6; %(mmHg to MPa)
        a.InitialParameters;
    end
else
    a = ArteryStrip;
    a.cs.ltG = 1.69; %Circumferential stretch ratio
end

%a.alphaPS = 40*pi/180; %rad
%a.cs.lambda = 1.58;

if ~a.InitialParameters
    SamplePoints = length(a.V.time);
    %Active Simulation
    for i=1:SamplePoints
        fprintf('(%d/%d) Time=%.1f min: ',i,SamplePoints,a.V.time(i));
        tic
        if a.stepCalc(i)
            fprintf('Error calculating current step\n');
            break;
        else
            eTime = toc;
            fprintf('(Step time is %.2f sec. Elapsed time is %.1f min)\n',eTime,(SamplePoints-i)*eTime/60);
        end
    end
    if i==SamplePoints
        a.PlotResults;
    end
end

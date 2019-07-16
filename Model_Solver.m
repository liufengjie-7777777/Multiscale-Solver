clear variables; close all; clc;

%Define a as artrey and calc n
if 1%input('Which experiment do you want to simulate?\nEnter 1 for biaxial and 0 for uniaxial: ')
    a = ArteryVessel;
else
    a = ArteryStrip;
end

err = a.InitialParameters;
if err
    fprintf('Error\n');
    return
end
SamplePoints = length(a.V.timeVec);

%Active Simulation
for i=1:SamplePoints
    fprintf('(%d/%d) Time=%.1f min: ',i,SamplePoints,a.V.timeVec(i));
    tic
    if a.stepCalc(i)
        fprintf('Error calculating current step\n');
        return
    end
    eTime = toc;
    fprintf('(Step time is %.2f sec. Elapsed time is %.1f min)\n',eTime,(SamplePoints-i)*eTime/60);
end

a.PlotResults;


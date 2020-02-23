%Simulation Program
clear all; close all; clc; %#ok<CLALL>

%Define a as artrey and calc n
if input('Which experiment do you want to simulate?\nEnter 1 for biaxial and 0 for uniaxial: ')
    a = ArteryVessel;
    a.cs.Pin = 90*133.322387415*1e-6; %Inner pressure
    a.cs.lambda = 1.5; %Axial loaded stretch ratio
    Sim = 'Biaxial';
else
    a = ArteryStrip;
    a.cs.ltG = 1.69; %Circumferential stretch ratio
    Sim = 'Biaxial';
end

varNames = {'lambda'};
saveDir = 'lz Simulations\';

for k=1:length(varNames)
    varName = varNames{k}; %'beta';
    mkdir(['Simulations\' saveDir varName ' Simulations\']);
    
    OriginalVarValue = a.cs.(varName);
    %Define values for sensitivity test
    varValues = linspace(0.5,1.5,6).*a.cs.(varName);
    
    a.PrintProgress = 0;
    
    if ~a.InitialParameters
        SamplePoints = length(a.V.time);
        if ~strcmp(Sim,'Biaxial')
            ufs = zeros(SamplePoints,31,length(varValues));
            ri = zeros(SamplePoints,length(varValues));
        end
    end
    
    for n=1:length(varValues)
        a.cs.(varName) = varValues(n);
        fprintf('(%d) Now simulating %s=%.2f\n',n,varName,varValues(n));
        if ~a.InitialParameters
            SamplePoints = length(a.V.time);
            %Active Simulation
            for i=1:SamplePoints
                if i==1 || (i>100 && mod(i,100)==0)
                    fprintf('(%d/%d) Time=%.1f min: ',i,SamplePoints,a.V.time(i));
                end
                tic
                if a.stepCalc(i)
                    fprintf('Error calculating current step\n');
                    i = SamplePoints+1; %#ok<FXSET>
                else
                    eTime = toc;
                    if i==1 || (i>100 && mod(i,100)==0)
                        fprintf('(Step time is %.2f sec. Elapsed time is %.1f min)\n',eTime,(SamplePoints-i)*eTime/60);
                    end
                end
            end
            save(['Simulations\' saveDir varName ' Simulations\' Sim 'Simulation-' varName '(' num2str(n) ').mat'],'a');
        end
    end
    a.cs.(varName) = OriginalVarValue;
end

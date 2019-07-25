%Simulation Program
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

varName = 'LMmax';
varUnits = 'nm'; varUnitConvertion = 1e6;

SimVariable = linspace(0.8,1.4,4)*a.(varName);
a.PrintProgress = 0;

if ~a.InitialParameters
	SamplePoints = length(a.V.time);
    ufs = zeros(SamplePoints,31,length(SimVariable));
    ri = zeros(SamplePoints,length(SimVariable));
end

for n=1:length(SimVariable)
    a.(varName) = SimVariable(n);
    fprintf('(%d) Now simulating %s = %.2f (%s)\n',n,varName,SimVariable(n)*varUnitConvertion,varUnits);
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
        
        ufs(:,:,n) = a.V.ufsN;
        ri(:,n) = a.V.ri;
        
        save(['Simulation-' varName '(' num2str(n) ').mat'],'a');
    end
end
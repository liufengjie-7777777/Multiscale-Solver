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
    Sim = 'Uniaxial';
end

paraName = 'lLA0';
paraVec = [7,10,15]*1e-6;

VarName = {'lambda','Pin'};
lambdaValues = linspace(1,1.6,5);
PinValues = linspace(80,110,4); %mmHg

saveDir = [paraName ' Simulation lzandPin\'];
mkdir(saveDir);

for k=1:length(paraVec)
    a.(paraName) = paraVec(k);
    
    a.PrintProgress = 0;
    
    SimCount = 1; 
    for n1=1:length(lambdaValues)
        for n2=1:length(PinValues)
            if ~a.InitialParameters
                SamplePoints = length(a.V.time);
            end

            a.cs.lambda = lambdaValues(n1);
            a.cs.Pin = PinValues(n2)*133.322387415*1e-6;

            fprintf('(%d) Now simulating %s_{(%d)}=%.3f %s=%.2f and %s=%.2f (mmHg)\n',SimCount,paraName,k,a.(paraName),'\lambda',lambdaValues(n1),'Pin',PinValues(n2));
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
                
                save([saveDir Sim paraName '(' num2str(k) ')' '_' 'Simulation' '(' num2str(SimCount) ').mat'],'a');
                SimCount = SimCount + 1;
            end
        end
    end

end

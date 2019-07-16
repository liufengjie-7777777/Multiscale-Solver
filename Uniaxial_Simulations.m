clear all; close all; clc; %#ok<CLALL>

lt = [1.4 1.6]; %0.5:0.1:2;

Pisom = zeros(length(lt),1);

for j=1:length(lt)
    %Define a as artrey and calc n
    a = ArteryStrip;
    a.cs.ltG = lt(j); %Circumferential stretch ratio
    a.cs.lrNum = 0;
    
    err = a.InitialParameters;
    if err
        fprintf('Error\n');
    else
        SamplePoints = length(a.V.timeVec);
        
        %Active Simulation
        for i=1:SamplePoints
            fprintf('(%d/%d) Time=%.1f min: ',i,SamplePoints,a.V.timeVec(i));
            tic
            if a.stepCalc(i)
                fprintf('Error calculating current step\n');
                break;
            end
            eTime = toc;
            fprintf('(Step time is %.2f sec. Elapsed time is %.1f min)\n',eTime,(SamplePoints-i)*eTime/60);
        end
        
        if i==SamplePoints
            Pisom(j) = double(a.V.PisomVec(end)*1e3); %kPa
        end
    end
end

plot(lt,Pisom);

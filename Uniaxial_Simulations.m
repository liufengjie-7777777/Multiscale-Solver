clear all; close all; clc; %#ok<CLALL>

lt = 1.2:0.01:2;
Pisom = zeros(1,length(lt));

%Define a as artrey
a = ArteryStrip;
%a.alphaPS = 25*pi/180;

for j=1:length(lt)
    a.cs.ltG = lt(j); %Circumferential stretch ratio
    a.cs.lrNum = 0; 
    
    err = a.InitialParameters;
    if err
        fprintf('Error\n');
    else
        SamplePoints = length(a.V.time);
        
        %Active Simulation
        for i=1:SamplePoints
            fprintf('(%d/%d) Time=%.1f min: ',i,SamplePoints,a.V.time(i));
            tic
            if a.stepCalc(i)
                fprintf('Error calculating current step\n');
                break;
            end
            eTime = toc;
            fprintf('(Step time is %.2f sec. Elapsed time is %.1f min)\n',eTime,(SamplePoints-i)*eTime/60);
        end
        
        if i==SamplePoints
            Pisom(j) = double(a.V.Pisom(end)*1e3); %kPa
        end
    end
end

N = 1;
%p = polyfit(lt(N:end),Pisom(N:end),2);
%xq = linspace(lt(N),lt(end),100);
%yq = polyval(p,xq);
%plot(lt(N:end),Pisom(N:end),'x',xq,yq);
plot(lt,Pisom);
ylabel('Steady State Pisom (kPa)'); xlabel('\lambda_\theta');
%ylim([0 25]);
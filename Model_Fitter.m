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

N = 3;
Cp = linspace(15,30,N)*1e-3;
c11 = linspace(2.5,4.5,N)*1e-3;
c21 = linspace(0.4,0.8,N);
c13 = linspace(8,10.5,N)*1e-3;
c23 = linspace(0.002,0.004,N);

for n1=1:N
    for n2=1:N
        for n3=1:N
            for n4=1:N
                for n5=1:N
                    a.Cp = Cp(n1);
                    a.c = [c11(n2),c11(n2),c13(n3),c13(n3) ; c21(n4),c21(n4),c23(n5),c23(n5)];
                    disp('asd');
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

                            if a.V.Pisom(i)*1e3 > 20
                                fprintf('\nTrying Next values (i=%d) %.4f kPa\n',i,a.V.Pisom(i)*1e3);
                                break;
                            end
                        end
                        if i==SamplePoints
                            a.PlotResults;
                        end
                    end
                end
            end
        end
    end
end

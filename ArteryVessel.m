classdef ArteryVessel < Artery
    properties
        dt = 1; %Seconds
        TotalTime = 15; %Minutes
        PrintProgress = 1;
    end
    methods
        %Length Dimensions Functions 
        function Ri = Ri(obj)
            Ri = obj.Ro-obj.H0; % inner radii
        end
        function hadv = hadv(obj)
            hadv = obj.H0/2; %Approximate diffusion distance through the Adventitia (same units as H0)
        end
        function r = r(obj)
        	r = sqrt( ((pi-obj.phi0/2)/(pi*obj.lz)) .* (obj.cs.R.^2-obj.Ri^2) + obj.cs.ri.^2 );
        end
        function R = R(obj)
        	R = sqrt( ((pi*obj.lz)/(pi-obj.phi0/2)) .* (obj.cs.r.^2-obj.cs.ri.^2) + obj.Ri^2  );
        end
        
        %Stretch Ratios
        function lr = lr(obj)
            lr = (pi-obj.phi0/2) .* obj.cs.R ./ (pi.*obj.cs.r.*obj.lz);
        end
        function lz = lz(obj)
            lz = obj.cs.lambda*obj.Deltaz;
        end
        function lt = lt(obj)
            lt = (pi/(pi-obj.phi0/2)).*(obj.cs.r./obj.cs.R);
        end
        
        %Some Shortcuts
        function obj = ro(obj,N)
            %N=1 for numeric ro
            if N %changing ri to its last known numeric value
                obj.cs.ri = obj.cs.riNum;
            end
            obj.cs.R = obj.Ro; %changing cs R to Ro
            if N %if numeric then save in cs roNum and return ri to sym
                obj.cs.roNum = obj.r;
                obj.cs.ri = obj.cs.riG;
            else %else save in cs ro
                obj.cs.ro = obj.r;
            end
            obj.cs.R = obj.R; %calc cs R again
        end
        
        %Calculate ri
        function err = riCalc(obj)
            %Calculates the new cs ri
            s = obj.cs.sECM + obj.cs.sSMC + obj.cs.sMMy;
            
            Sr = s(1,:); St = s(2,:); Sz = s(3,:);
            
            f = sum( ((St-Sr)./obj.cs.x).*obj.cs.w );
            
            if obj.cs.riNum %Use the previous riNum value for the numeric sol
            	riNew = double(vpasolve(f == obj.cs.Pin,[0.8*obj.cs.riNum 1.2*obj.cs.riNum]));
            else
            	riNew = double(vpasolve(f == obj.cs.Pin,[0 Inf]));
            end
            
            if isempty(riNew)
                err = 1;
            else
                obj.cs.riNum = riNew; %Updates riNum
                obj.xwNum; %Updates numeric x and w
                obj.cs.FT = obj.FTCalc(Sr,St,Sz); %Calcs F_T
                %Updates the new numeric values of the stretches
                obj.cs.r = obj.cs.xNum; %Temporarily change r and ri to numeric values
                obj.cs.ri = obj.cs.riNum;
                obj.cs.R = obj.R; %Updates R to numerics values
                
                obj.cs.lrNum = obj.lr;
                obj.cs.ltNum = obj.lt;
                
                %Update r, ri and R to sym values
                obj.cs.r = obj.cs.x;
                obj.cs.ri = obj.cs.riG;
                obj.cs.R = obj.R;
                
                err = 0;
            end
        end
        
        function err = stepCalc(obj,i)
            %Update current Step Myosin Fractions
            obj.cs.nAMp = obj.V.nAMp(i);
            obj.cs.nAM = obj.V.nAM(i);
            %Update ufs (also calculate and save numeric values of
            %dMAiMA0,LMi,Lfoi and eS2)
            obj.ufsUpdate;
            %Update the active stress components
            obj.sECM;
            obj.sSMC;
            obj.sMMy;
            
            if obj.riCalc %Calculate new cs riNum and FT
                fprintf('Error calculating ri\n');
                err = 1;
            else
                obj.V.UpdateVectors(i,obj.cs);
                obj.ufsVecUpdate(i);
                if obj.PrintProgress
                    fprintf('| Do=%.2f kPa | F_T=%.2f mN | ',obj.V.Do(i),obj.V.FT(i));
                end
                err = 0;
            end
        end
        
        function [err] = InitialParameters(obj)
            obj.nCalc;
            
            %2nd Step - Calculate Passive State
            obj.cs.ri = obj.cs.riG; %Update ri to be symbolic
            
            obj.xw; %Calculate Gaussian abscicess
            obj.cs.r = obj.cs.x; %Update current radii vector
            obj.cs.R = obj.R; %Update initial radii vector accordingly
            
            obj.cs.lr = obj.lr;
            obj.cs.lt = obj.lt;
            obj.cs.lz = obj.lz;
            
            obj.cs.lzNum = obj.lz;
            
            obj.sECM;
            obj.cs.sSMC = zeros(size(obj.cs.sECM));
            obj.cs.sMMy = zeros(size(obj.cs.sECM));
            
            if obj.riCalc
                fprintf('Error calculating ri\n');
                err = 1;
            else
                obj.ufs0; %calcs numeric ufs0
                obj.ufsVecUpdate(1);
                fprintf('Initial Passive Conditions: ');
                fprintf('Do=%.3f, lr=%.3f, lt=%.3f, lz=%.3f, detF=%.3f \n',obj.cs.roNum*2e3,obj.cs.lrNum(2),obj.cs.ltNum(2),obj.cs.lzNum,(obj.cs.lrNum(2)*obj.cs.ltNum(2)*obj.cs.lzNum));
                obj.V.InitialVectors(length(obj.V.time),1);
                err = 0;
            end
        end
        
        %Myosin Kinteics Functions
        function obj = nCalc(obj)
            function dydt = myode(t,y)
                k1 = obj.k1t(t);
                kM = [-k1, obj.k(2), 0, obj.k(7)
                    k1, -obj.k(2)-obj.k(3), obj.k(4), 0
                    0, obj.k(3), -obj.k(4)-obj.k(5), obj.k(6)
                    0, 0, obj.k(5), -obj.k(6)-obj.k(7)];
                dydt = kM*y;
            end
            [time, n] = ode15s(@(t,y) myode(t,y),0:obj.dt:obj.TotalTime*60,[1 0 0 0]);
            
            obj.V.time = time./60;
            obj.V.nAMp = n(:,3); obj.V.nAM = n(:,4);
        end
        function k1t = k1t(obj,t)
            k1t = obj.k(1)*erfc( obj.hadv/(2*sqrt(obj.DKCL*t)) );
        end
        
        %Calculate Axial Force for a Given Stress
        function FT = FTCalc(obj,Sr,St,Sz)
            f = subs(2*Sz-St-Sr,obj.cs.riNum);
            FT = double( pi*sum((f.*obj.cs.xNum).*obj.cs.wNum) );
        end
        
        %Calculate Legendre Points and Weights
        function obj = xw(obj)
            obj.ro(0);
            [obj.cs.x,obj.cs.w] = lgwt(3,obj.cs.ri,obj.cs.ro);
        end
        function obj = xwNum(obj)
            obj.ro(1);
            [obj.cs.xNum,obj.cs.wNum] = lgwt(3,obj.cs.riNum,obj.cs.roNum);
        end
        
        %Some Additional Shortcuts for Results Analysis       
        function obj = PlotResults(obj)
            figure(1);
            plot(obj.V.time,obj.V.Do);
            grid on; ylim([600 1300]); xlim([0 obj.TotalTime]); %ylim([0 (ceil(max(obj.Do))+100)]);
            ylabel('Do (um)'); xlabel('time (min)');
            title(['lz=' num2str(obj.lz) ' Pin=' num2str(obj.cs.Pin/133.322387415*1e6) ' mmHg']);
            
            figure(2);
            plot(obj.V.time,obj.V.FT);
            minFT = min(obj.V.FT);
            maxFT = max(obj.V.FT);
            if minFT>2
                minFT = ceil(minFT-2);
            else
               minFT = 0;
            end
            if maxFT>0
                maxFT = ceil(maxFT+2);
            else
                maxFT = 0;
                if minFT >= 0
                    minFT = floor(min(obj.V.FT));
                end
            end
            grid on; ylim([minFT maxFT]); xlim([0 obj.TotalTime]);
            ylabel('F_T (mN)'); xlabel('time (min)');
            
            figure(3);
            plot(obj.V.time,obj.V.stretch);
            grid on; ylim([0 2]); xlim([0 obj.TotalTime]);
            ylabel('Stretch Ratio'); xlabel('time (min)');
            legend('\lambda_r','\lambda_\theta','\lambda_z','det(F)');

            figure(4);
            plot(obj.V.time,obj.V.PMMCU*1e3);
            ylabel('Stress (kPa)'); ylim([0 ceil(max(obj.V.PMMCU*1e3,[],'all')+10)]);
            hold on; yyaxis right;
            plot(obj.V.time,obj.V.ufs);
            h = legend('$P_{MM}$','$P_{CU}$','$\bar{u}_{fs}$ (right axis)');
            ylabel('ufs'); xlabel('time (min)'); grid on; hold off;
            set(h,'Interpreter','latex','fontsize',12);   

            figure(5);
            plot(obj.V.time, (obj.V.sECM+obj.V.sSMC+obj.V.sMMy)*1e3);
            h = legend('$\bar{\sigma}_{r}$','$\bar{\sigma}_{\theta}$','$\bar{\sigma}_{z}$');
            ylabel('Cauchy Stress (kPa)'); xlabel('time (min)'); grid on;
            set(h,'Interpreter','latex','fontsize',12);
        end
        
        function ufsVecUpdate(obj,N)
            r = linspace(obj.cs.riNum,obj.cs.roNum,31);
            
            R = sqrt( ((pi*obj.lz)/(pi-obj.phi0/2)) .* (r.^2-obj.cs.riNum.^2) + obj.Ri^2  );
            lt = (pi/(pi-obj.phi0/2)).*(r./R);
            lr = 1./(lt.*obj.lz);
            
            if N<2
                obj.V.ufsN(1,:) = sqrt( (lr.^2)*power(sin(obj.thetaSMC),2) + (lt.^2)*power(cos(obj.thetaSMC),2) );
            else
                strain = vertcat(obj.cs.lr,obj.cs.lt);
                obj.cs.lr = lr;
                obj.cs.lt = lt;
                obj.cs.lz = obj.lz;
                ufs = obj.cs.ufs;
                obj.cs.ufs = obj.V.ufsN(N-1,:);
                
                obj.Lfoi;
                obj.eS2;
                
                ufsdot = double(obj.beta*(obj.PCU-obj.PMM));
                if ufsdot<0
                    obj.V.ufsN(N,:) = obj.V.ufsN(N-1,:) + ufsdot*obj.dt; %#ok<*MCNPN>
                end
                
                obj.cs.lr = obj.cs.lrNum;
                obj.cs.lt = obj.cs.ltNum;
                obj.cs.ufs = ufs;
                
                obj.Lfoi;
                obj.eS2;
                obj.PCU;
                obj.PMM;
                
                obj.cs.lr = strain(1,:);
                obj.cs.lt = strain(2,:);
            end
        end
    end
end

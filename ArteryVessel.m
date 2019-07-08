classdef ArteryVessel < Artery
    properties
        %Current State
        ri = sym('ri','positive'); %Inner Radius
        Pin = 90*133.322387415*1e-6; %Inner Pressure
        lambda = 1.25;
        
        riNum = 0;
        FT
        
        %Add current state class
        
        dt = 1; %Seconds
        TotalTime = 15; %Minutes
        
        %Vectors
        V = SimulationVectors;
    end
    methods
        %Radius and etc.
        function Ri = Ri(obj)
            Ri = obj.Ro-obj.H0; % inner radii
        end
        function hadv = hadv(obj)
            hadv = obj.H0/2; %Approximate diffusion distance through the Adventitia (same units as H0)
        end
        function r = r(obj,R)
        	r = sqrt( ((pi-obj.phi0/2)/(pi*obj.lz)) .* (R.^2-obj.Ri^2) + obj.ri.^2 );
        end
        function R = R(obj,r)
        	R = sqrt( ((pi*obj.lz)/(pi-obj.phi0/2)) .* (r.^2-obj.ri.^2) + obj.Ri^2  );
        end
        
        %Stretch Ratios
        function lr = lr(obj,r)
            lr = (pi-obj.phi0/2) .* obj.R(r) ./ (pi.*r.*obj.lz);
        end
        function lz = lz(obj)
            lz = obj.lambda*obj.Deltaz;
        end
        function lt = lt(obj,r)
            lt = (pi/(pi-obj.phi0/2)).*(r./obj.R(r));
        end
        
        %Some Shortcuts
        function ro = ro(obj)
        	ro = obj.r(obj.Ro);
        end
        function lto = lto(obj)
            lto = obj.lt(obj.ro);
        end
        function lti = lti(obj)
        	lti = obj.lt(obj.ri);
        end
        function lro = lro(obj)
            lro = obj.lr(obj.ro);
        end
        
        %Numerical Shortcuts
        function RNum = RNum(obj,r)
        	RNum = sqrt( ((pi*obj.lz)/(pi-obj.phi0/2)) .* (r.^2-obj.riNum.^2) + obj.Ri^2  );
        end
        function roNum = roNum(obj)
            roNum = sqrt( ((pi-obj.phi0/2)/(pi*obj.lz)) .* (obj.Ro.^2-obj.Ri^2) + obj.riNum.^2 );
        end
        function lrNum = lrNum(obj,r)
            lrNum = (pi-obj.phi0/2) .* obj.RNum(r) ./ (pi.*r.*obj.lz);
        end
        function ltNum = ltNum(obj,r)
            ltNum = (pi/(pi-obj.phi0/2)).*(r./obj.RNum(r));
        end
        
        %Calculate ri
        function err = riCalc(obj,sECM,sSMC,sMMy)
            [x,w] = obj.xw;
            if exist('sECM','var')
                s = sECM + sSMC + sMMy;
            else
                if obj.ufs ~= 0
                	[xNum,~] = obj.xwNum;
                    
                    [LMr,LMz] = obj.LMi(xNum);
                    [Lfor,Lfoz] = obj.Lfoi(xNum);
                    [eS2xr,eS2yr] = obj.eS2(1,xNum);
                    [eS2xz,eS2yz] = obj.eS2(2,xNum);
                    I4SMCe = subs(obj.I4SMCe(xNum),obj.riNum);
                    lr = obj.lrNum(xNum);
                    lt = obj.ltNum(xNum);
                    
                    obj.ufsUpdate(LMr,Lfor,LMz,Lfoz,eS2xr,eS2xz,I4SMCe,lr,lt,obj.ufsfit(xNum));
                    
                    I4SMCe = obj.I4SMCe(x);

                    s = obj.sECM(x) + obj.sSMC(I4SMCe,LMr,Lfor,LMz,Lfoz,lr,lt,obj.ufs) + obj.sMMy(LMr,Lfor,LMz,Lfoz,eS2yr,eS2yz,lr,lt);
                end
            end
            Sr = s(1,:); St = s(2,:); Sz = s(3,:);
            
            f = sum( ((St-Sr)./x).*w );
            
            if obj.riNum
            	riNew = double(vpasolve(f == obj.Pin,[0.9*obj.riNum 1.1*obj.riNum]));
            else
            	riNew = double(vpasolve(f == obj.Pin,[0 2]));
            end
            
            if isempty(riNew)
                fprintf('Error calculating ri\n');
                err = 1;
            else
                obj.riNum = riNew;
                obj.FT = obj.FTCalc(Sr,St,Sz);
                
                if obj.ufs == 0
                    [xNum,~] = obj.xwNum;
                    obj.ufs = sqrt( (obj.lrNum(xNum).^2)*power(sin(obj.thetaSMC),2) + (obj.ltNum(xNum).^2)*power(cos(obj.thetaSMC),2) );
                end
                err = 0;
            end
        end
        
        function err = stepCalc(obj,i,sECM)
            %Current Step Myosin Fractions
            obj.nAMp = obj.nAMpVec(i);
            obj.nAM = obj.nAMVec(i);
            
            [xNum,~] = obj.xwNum;

            [LMr,LMz] = obj.LMi(xNum);
            [Lfor,Lfoz] = obj.Lfoi(xNum);
            [eS2xr,eS2yr] = obj.eS2(1,xNum);
            [eS2xz,eS2yz] = obj.eS2(2,xNum);
            I4SMCe = subs(obj.I4SMCe(xNum),obj.riNum);
            lr = obj.lrNum(xNum);
            lt = obj.ltNum(xNum);
            
            obj.ufsUpdate(LMr,Lfor,LMz,Lfoz,eS2xr,eS2xz,I4SMCe,lr,lt,obj.ufsfit(xNum));
            [x,~] = obj.xw;
            
            I4SMCe = obj.I4SMCe(x);
%             [LMr,LMz] = obj.LMi(x);
%             [Lfor,Lfoz] = obj.Lfoi(x);
%             lr = obj.lr(x);
%             lt = obj.lt(x);
            
            if ~exist('sECM','var')
                sECM = zeros(3,length(x))*obj.ri;
                for j=1:length(x)
                    sECM(:,j) = obj.sECM(x(j));
                end
            end
            sSMC = obj.sSMC(I4SMCe,LMr,Lfor,LMz,Lfoz,lr,lt,obj.ufs);
            sMMy = obj.sMMy(LMr,Lfor,LMz,Lfoz,eS2yr,eS2yz,lr,lt);
            
            if obj.riCalc(sECM,sSMC,sMMy)
                fprintf('Error calculating ri\n');
                err = 1;
            else
                rm = (obj.riNum+obj.roNum)/2; %Middle radii

                PMMCU = subs([obj.PMM(LMr(2),Lfor(2),LMz,Lfoz(2),eS2xr(2),eS2xz) obj.PCU(LMr(2),Lfor(2),LMz,Lfoz(2),I4SMCe(2),lr(2),lt(2),obj.ufsfit(rm))],obj.riNum);
                
                obj.V.UpdateVectors(i,obj.ufsfit(rm),[obj.lrNum(rm),obj.ltNum(rm),obj.lz],subs(sECM(:,2),obj.riNum),subs(sSMC(:,2),obj.riNum),subs(sMMy(:,2),obj.riNum),PMMCU,0,obj.roNum*2e3,obj.FT*1e3);
                
                fprintf('| Do=%.2f kPa | F_T=%.2f mN | ',obj.V.DoVec(i),obj.V.FTVec(i));
                err = 0;
            end
        end
        
        function [err,sECM] = InitialParameters(obj)
            [obj.timeVec,n] = obj.nCalc;
            obj.nAMpVec = n(:,3); obj.nAMVec = n(:,4);

            [x,~] = obj.xw; sECM = zeros(3,length(x))*obj.ri;
            for j=1:length(x)
                sECM(:,j) = obj.sECM(x(j));
            end
            if obj.riCalc(sECM,0,0)
                fprintf('Error calculating ri\n');
                err = 1;
            else
                fprintf('Initial Passive Conditions: ');
                fprintf('Do=%.3f, lr=%.3f, lt=%.3f, lz=%.3f \n',obj.roNum*2e3,obj.lrNum(obj.riNum),obj.ltNum(obj.riNum),obj.lz);
                obj.V.InitialVectors(length(obj.timeVec),0);
                err = 0;
            end
        end
        
        function [time,n] = nCalc(obj)
            function dydt = myode(t,y)
                k1 = obj.k1t(t);
                kM = [-k1, obj.k(2), 0, obj.k(7)
                    k1, -obj.k(2)-obj.k(3), obj.k(4), 0
                    0, obj.k(3), -obj.k(4)-obj.k(5), obj.k(6)
                    0, 0, obj.k(5), -obj.k(6)-obj.k(7)];
                dydt = kM*y;
            end
            [time, n] = ode15s(@(t,y) myode(t,y),0:obj.dt:obj.TotalTime*60,[1 0 0 0]);
            
            obj.timeVec = time./60;
            obj.nAMpVec = n(:,3); obj.nAMVec = n(:,4);
        end
        function k1t = k1t(obj,t)
            k1t = obj.k(1)*erfc( obj.hadv/(2*sqrt(obj.DKCL*t)) );
        end
        
        %Calculate Axial Force for a Given Stress
        function FT = FTCalc(obj,Sr,St,Sz)
            [x,w] = obj.xwNum;
            if exist('Sr','var')
                f = subs(2*Sz-St-Sr,obj.riNum);
            end
            FT = double( pi*sum((f.*x).*w) );
        end
        
        %Calculate Legendre Points and Weights
        function [x,w] = xw(obj)
            [x,w] = lgwt(3,obj.ri,obj.ro);
        end
        function [x,w] = xwNum(obj)
            [x,w] = lgwt(3,obj.riNum,obj.roNum);
        end
        
        %Some Additional Shortcuts for Results Analysis
        function ufsfit = ufsfit(obj,r)
            p = polyfit(obj.xwNum,obj.ufs,1);
            ufsfit = p(1).*r + p(2);
        end
        
        function obj = PlotResults(obj)
            figure(1);
            plot(obj.timeVec./60,obj.V.DoVec);
            grid on; ylim([600 1300]); xlim([0 obj.TotalTime]); %ylim([0 (ceil(max(obj.DoVec))+100)]);
            ylabel('Do (um)'); xlabel('time (min)');
            title(['lz=' num2str(obj.lz) ' Pin=' num2str(obj.Pin/133.322387415*1e6) ' mmHg']);
            
            figure(2);
            plot(obj.timeVec./60,obj.V.FTVec);
            minFT = min(obj.V.FTVec);
            maxFT = max(obj.V.FTVec);
            if minFT<0
                minFT = ceil(minFT-2);
            else
               minFT = 0; 
            end
            if maxFT<0
                maxFT = ceil(maxFT+2);
            else
                maxFT = 0;
            end
            grid on; ylim([minFT maxFT]); xlim([0 obj.TotalTime]);
            ylabel('F_T (mN)'); xlabel('time (min)');
            
            figure(3);
            plot(obj.timeVec./60,obj.V.stretchVec);
            grid on; ylim([0 2]); xlim([0 obj.TotalTime]);
            ylabel('Stretch Ratio'); xlabel('time (min)');
            legend('\lambda_r','\lambda_\theta','\lambda_z','det(F)');

            figure(4);
            plot(obj.timeVec./60,obj.V.PMMCUVec*1e3);
            ylabel('Stress (kPa)'); ylim([0 ceil(max(obj.V.PMMCUVec*1e3,[],'all')+10)]);
            hold on; yyaxis right;
            plot(obj.timeVec./60,obj.V.ufsVec);
            h = legend('$P_{MM}$','$P_{CU}$','$\bar{u}_{fs}$ (right axis)');
            ylabel('ufs'); xlabel('time (min)'); grid on; hold off;
            set(h,'Interpreter','latex','fontsize',12);   

%             figure(4);
%             plot(obj.timeVec./60, (obj.sECMVec+obj.sSMCVec+obj.sMMyVec)*1e3);
%             h = legend('${\sigma}_{r}$','${\sigma}_{\theta}$','${\sigma}_{z}$');
%             ylabel('Cauchy Stress (kPa)'); xlabel('time (min)'); grid on;
%             set(h,'Interpreter','latex','fontsize',12);
        end
    end
end

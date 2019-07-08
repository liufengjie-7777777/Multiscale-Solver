classdef ArteryStrip < Artery
    properties
        %Current State
        lr = sym('lr','positive');
        lt = 1.69;
        lrNum
        
        dt = 1; %Seconds
        TotalTime = 5; %Minutes
        
        %Vectors
        V = SimulationVectors;
    end
    methods
        %Stretch Ratios
        function lz = lz(obj)
            lz = 1./(obj.lt.*obj.lr);
        end
        function lzNum = lzNum(obj)
            lzNum = 1./(obj.lt.*obj.lrNum);
        end
        
        function Pisom = Pisom(obj,LMr,LMz,Lfor,Lfoz,eS2yr)
            if ~exist('LMr','var')
                [LMr,LMz] = obj.LMi;
                [Lfor,Lfoz] = obj.Lfoi;
                [~,eS2yr] = obj.eS2(1);
            end
            PisomSMC = (obj.EAM*obj.nAM + obj.EAMp*obj.nAMp)*obj.AS2*obj.NCF*...
                (cos(obj.thetaSMC)^2).*obj.lt.*...
                (sqrt(obj.I4SMCe)-1).*(LMr.*Lfor*obj.gammar + LMz.*Lfoz*(1-obj.gammar))...    
                ./(2.*sqrt(obj.I4SMCe).*power(obj.ufs,2)*obj.deltam);
        
            PMMyt = (obj.EAM*obj.nAM + obj.EAMp*obj.nAMp)*obj.AS2*obj.NCU*...
                (sin(obj.thetaSMC)^2)*LMr.*Lfor.*eS2yr*obj.gammar/obj.deltam;
            Pisom = subs(PisomSMC + PMMyt,obj.lrNum);
        end
        
        %Calculate lr
        function err = lrCalc(obj,sECM,sSMC,sMMy)
            if ~exist('sECM','var')
                s = obj.sECM;
                if obj.ufs ~= 0
                    obj.ufsUpdate;
                    s = s + obj.sSMC + obj.sMMy;
                end
            else
                s = sECM + sSMC + sMMy;
            end
            Sr = s(1); Sz = s(3);
            lrNew = vpasolve(Sr == Sz,obj.lr,[0 2]);
            
            if isempty(lrNew)
                fprintf('Error calculating lr\n');
                err = 1;
            else
                err = 0;
                obj.lrNum = double(lrNew);
                if obj.ufs == 0
                    obj.ufs = sqrt( (lrNew.^2)*power(sin(obj.thetaSMC),2) + (obj.lt.^2)*power(cos(obj.thetaSMC),2) );
                end
            end
        end
        
        function err = stepCalc(obj,i,sECM)
            %Current Step Myosin Fractions
            obj.nAMp = obj.nAMpVec(i);
            obj.nAM = obj.nAMVec(i);
            
            [LMr,LMz] = obj.LMi;
            [Lfor,Lfoz] = obj.Lfoi;
            [eS2xr,eS2yr] = obj.eS2(1);
            [eS2xz,eS2yz] = obj.eS2(2);
            I4SMCe = subs(obj.I4SMCe,obj.lrNum);
            obj.ufsUpdate(LMr,Lfor,LMz,Lfoz,eS2xr,eS2xz,I4SMCe,obj.lrNum,obj.lt,obj.ufs);
            
            if ~exist('sECM','var')
                sECM = obj.sECM;
            end
            sSMC = obj.sSMC(obj.I4SMCe,LMr,Lfor,LMz,Lfoz,obj.lrNum,obj.lt,obj.ufs);
            sMMy = obj.sMMy(LMr,Lfor,LMz,Lfoz,eS2yr,eS2yz,obj.lrNum,obj.lt);
            
            if obj.lrCalc(sECM,sSMC,sMMy)
                fprintf('Error calculating lr\n');
                err = 1;
            else
                obj.V.UpdateVectors(i,obj.ufs,[obj.lrNum,obj.lt,obj.lzNum],subs(sECM,obj.lrNum),subs(sSMC,obj.lrNum),subs(sMMy,obj.lrNum),subs([obj.PMM, obj.PCU],obj.lrNum),obj.Pisom(LMr,LMz,Lfor,Lfoz,eS2yr));
                
                fprintf('| Pisom=%.2f kPa | ',obj.V.PisomVec(i)*1e3);
                fprintf('lr=%.3f, lt=%.3f, lz=%.3f | ',obj.V.stretchVec(i,1:3));
                err = 0;
            end    
        end
        
        function [err,sECM] = InitialParameters(obj)
            [obj.timeVec,n] = obj.nCalc;
            obj.nAMpVec = n(:,3); obj.nAMVec = n(:,4);
            %2nd Step - Calculate Passive State
            sECM = obj.sECM;
            if obj.lrCalc(sECM,0,0)
                fprintf('Error calculating lr\n');
                err = 1;
            else
                fprintf('Initial Passive Conditions: ');
                fprintf('lr=%.3f, lt=%.3f, lz=%.3f, det(F)=%.3f \n',obj.lrNum,obj.lt,obj.lzNum,(obj.lrNum*obj.lt*obj.lzNum));
                obj.V.InitialVectors(length(obj.timeVec),1);
                err = 0;
            end
        end
        
        function [time,n] = nCalc(obj)
            function dydt = myode(~,y)
                kM = [-obj.k(1) obj.k(2) 0 obj.k(7)
                    obj.k(1) -obj.k(2)-obj.k(3) obj.k(4) 0
                    0 obj.k(3) -obj.k(4)-obj.k(5) obj.k(6)
                    0 0 obj.k(5) -obj.k(6)-obj.k(7)];
                dydt = kM*y;
            end
            [time, n] = ode15s(@(t,y) myode(t,y),0:obj.dt:obj.TotalTime*60,[1 0 0 0]);
            
            obj.timeVec = time./60;
            obj.nAMpVec = n(:,3); obj.nAMVec = n(:,4);
        end
        
        function obj = PlotResults(obj)
%             for i=1:length(obj.sECMVec)
%                 obj.sECMVec(i,:) = double(subs(obj.sECMVec(i,:),obj.stretchVec(i,1)));
%                 obj.sSMCVec(i,:) = double(subs(obj.sSMCVec(i,:),obj.stretchVec(i,1)));
%                 obj.sMMyVec(i,:) = double(subs(obj.sMMyVec(i,:),obj.stretchVec(i,1)));
%             end
            pVec = obj.V.sECMVec(:,1) + obj.V.sSMCVec(:,1) + obj.V.sMMyVec(:,1);
            
            figure(1);
            plot(obj.timeVec./60,obj.V.PisomVec*1e3);
            grid on; ylim([0 (ceil(max(obj.V.PisomVec*1e3))+10)]); xlim([0 obj.TotalTime]);
            ylabel('Pisom (kPa)'); xlabel('time (min)');
            title(['lt=' num2str(obj.lt)]);

            figure(2);
            plot(obj.timeVec./60,obj.V.stretchVec);
            grid on; ylim([0 2]); xlim([0 obj.TotalTime]);
            ylabel('Stretch Ratio'); xlabel('time (min)');
            legend('\lambda_r','\lambda_\theta','\lambda_z','det(F)');

            figure(3);
            plot(obj.timeVec./60,obj.V.PMMCUVec*1e3);
            ylabel('Stress (kPa)'); ylim([0 50]);
            hold on; yyaxis right;
            plot(obj.timeVec./60,obj.V.ufsVec);
            h = legend('$P_{MM}$','$P_{CU}$','$\bar{u}_{fs}$ (right axis)');
            ylabel('ufs'); xlabel('time (min)'); grid on; hold off;
            set(h,'Interpreter','latex','fontsize',12);   

            figure(4);
            plot(obj.timeVec./60, (obj.V.sECMVec+obj.V.sSMCVec+obj.V.sMMyVec-pVec)*1e3);
            h = legend('${\sigma}_{r}$','${\sigma}_{\theta}$','${\sigma}_{z}$');
            ylabel('Cauchy Stress (kPa)'); xlabel('time (min)'); grid on;
            set(h,'Interpreter','latex','fontsize',12);
        end
    end
end
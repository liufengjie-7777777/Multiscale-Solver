classdef ArteryStrip < Artery
    properties
        dt = 1; %Seconds
        TotalTime = 5; %Minutes 
    end
    methods
        %Stretch Ratios
        function lz = lz(obj)
            lz = 1./(obj.cs.lt.*obj.cs.lr);
        end
        function lzNum = lzNum(obj)
            lzNum = 1./(obj.cs.ltNum.*obj.cs.lrNum);
        end
        
        function obj = Pisom(obj)
            PisomSMC = (obj.EAM*obj.cs.nAM + obj.EAMp*obj.cs.nAMp)*obj.AS2*obj.NCF*...
                (cos(obj.thetaSMC)^2).*obj.cs.ltNum.*...
                (sqrt(obj.cs.I4SMCeNum)-1).*(obj.cs.LMr.*obj.cs.Lfor*obj.gammar + obj.cs.LMz.*obj.cs.Lfoz*(1-obj.gammar))...    
                ./(2.*sqrt(obj.cs.I4SMCeNum).*power(obj.cs.ufs,2)*obj.deltam);
            
            PMMyt = (obj.EAM*obj.cs.nAM + obj.EAMp*obj.cs.nAMp)*obj.AS2*obj.NCU*...
                (sin(obj.thetaSMC)^2)*obj.cs.LMr.*obj.cs.Lfor.*obj.cs.eS2yr*obj.gammar/obj.deltam;
            obj.cs.Pisom = subs(PisomSMC + PMMyt,obj.cs.lrNum);
        end
        
        %Calculate lr
        function err = lrCalc(obj)
            s = obj.cs.sECM + obj.cs.sSMC + obj.cs.sMMy;
            
            Sr = s(1); Sz = s(3);
            if obj.cs.lrNum
                lrNew = vpasolve(Sr == Sz,obj.cs.lrG,[0.9*obj.cs.lrNum 1.1*obj.cs.lrNum]);
            else
                lrNew = vpasolve(Sr == Sz,obj.cs.lrG,[0 2]);
            end
            if isempty(lrNew)
                fprintf('Error calculating lr\n');
                err = 1;
            else
                obj.cs.lrNum = double(lrNew);
                obj.cs.lzNum = obj.lzNum;
                if obj.cs.ufs == 0
                    obj.cs.ufs = sqrt( (obj.cs.lrNum.^2)*power(sin(obj.thetaSMC),2) + (obj.cs.ltNum.^2)*power(cos(obj.thetaSMC),2) );
                end
                
                err = 0;
            end
        end
        
        function err = stepCalc(obj,i)
            %Current Step Myosin Fractions
            obj.cs.nAMp = obj.V.nAMpVec(i);
            obj.cs.nAM = obj.V.nAMVec(i);
            
            obj.ufsUpdate;
            
            obj.sSMC;
            obj.sMMy;
            
            if obj.lrCalc
                fprintf('Error calculating lr\n');
                err = 1;
            else
                obj.Pisom;
                
                obj.V.UpdateVectors(i,obj.cs);
                fprintf('| Pisom=%.2f kPa | ',obj.V.PisomVec(i)*1e3);
                fprintf('lr=%.3f, lt=%.3f, lz=%.3f | ',obj.V.stretchVec(i,1:3));
                err = 0;
            end    
        end
        
        function [err] = InitialParameters(obj)
            obj.nCalc;
            
            %2nd Step - Calculate Passive State
            obj.cs.lr = obj.cs.lrG;
            obj.cs.lt = obj.cs.ltG;
            obj.cs.lz = obj.lz;
            
            obj.cs.ltNum = obj.cs.lt;
            
            obj.sECM;
            obj.cs.sSMC = zeros(3,1);
            obj.cs.sMMy = zeros(3,1);
            
            %obj.cs.I4SMCe = 1;
            
            if obj.lrCalc
                fprintf('Error calculating lr\n');
                err = 1;
            else
                fprintf('Initial Passive Conditions: ');
                fprintf('lr=%.3f, lt=%.3f, lz=%.3f, det(F)=%.3f \n',obj.cs.lrNum,obj.cs.ltNum,obj.cs.lzNum,(obj.cs.lrNum*obj.cs.ltNum*obj.cs.lzNum));
                obj.V.InitialVectors(length(obj.V.timeVec),1);
                err = 0;
            end
        end
        
        function obj = nCalc(obj)
            function dydt = myode(~,y)
                kM = [-obj.k(1) obj.k(2) 0 obj.k(7)
                    obj.k(1) -obj.k(2)-obj.k(3) obj.k(4) 0
                    0 obj.k(3) -obj.k(4)-obj.k(5) obj.k(6)
                    0 0 obj.k(5) -obj.k(6)-obj.k(7)];
                dydt = kM*y;
            end
            [time, n] = ode15s(@(t,y) myode(t,y),0:obj.dt:obj.TotalTime*60,[1 0 0 0]);
            
            obj.V.timeVec = time./60;
            obj.V.nAMpVec = n(:,3); obj.V.nAMVec = n(:,4);
        end
        
        function obj = PlotResults(obj)
%             for i=1:length(obj.sECMVec)
%                 obj.sECMVec(i,:) = double(subs(obj.sECMVec(i,:),obj.stretchVec(i,1)));
%                 obj.sSMCVec(i,:) = double(subs(obj.sSMCVec(i,:),obj.stretchVec(i,1)));
%                 obj.sMMyVec(i,:) = double(subs(obj.sMMyVec(i,:),obj.stretchVec(i,1)));
%             end
            pVec = obj.V.sECMVec(:,1) + obj.V.sSMCVec(:,1) + obj.V.sMMyVec(:,1);
            
            figure(1);
            plot(obj.V.timeVec,obj.V.PisomVec*1e3);
            grid on; ylim([0 (ceil(max(obj.V.PisomVec*1e3))+10)]); xlim([0 obj.TotalTime]);
            ylabel('Pisom (kPa)'); xlabel('time (min)');
            title(['lt=' num2str(obj.cs.ltNum)]);
            
            figure(2);
            plot(obj.V.timeVec,obj.V.stretchVec);
            grid on; ylim([0 2]); xlim([0 obj.TotalTime]);
            ylabel('Stretch Ratio'); xlabel('time (min)');
            legend('\lambda_r','\lambda_\theta','\lambda_z','det(F)');
            
            figure(3);
            plot(obj.V.timeVec,obj.V.PMMCUVec*1e3);
            ylabel('Stress (kPa)'); ylim([0 50]);
            hold on; yyaxis right;
            plot(obj.V.timeVec,obj.V.ufsVec);
            h = legend('$P_{MM}$','$P_{CU}$','$\bar{u}_{fs}$ (right axis)');
            ylabel('ufs'); xlabel('time (min)'); grid on; hold off;
            set(h,'Interpreter','latex','fontsize',12);
            
            figure(4);
            plot(obj.V.timeVec, (obj.V.sECMVec+obj.V.sSMCVec+obj.V.sMMyVec-pVec)*1e3);
            h = legend('${\sigma}_{r}$','${\sigma}_{\theta}$','${\sigma}_{z}$');
            ylabel('Cauchy Stress (kPa)'); xlabel('time (min)'); grid on;
            set(h,'Interpreter','latex','fontsize',12);
        end
    end
end
classdef Artery < handle
    properties
        %Current State
        p = sym('p','real');
        
        ufs = 0;
        nAMp = 0;
        nAM = 0;
        
        %Vectors
        nAMpVec
        nAMVec
        timeVec
        
        %Stress Free Parameters
        L0 = 5.1; %mm
        Ro = 409*1e-3; %um to mm
        H0 = 102*1e-3; %um to mm
        phi0 = 12 * pi/180; %deg to rad
        Deltaz = 1; %My decision
        
        %Passive ECM and Collagen Parameters
        Cp = 27.6761 *1e-3; %kPa to MPa
        c = [ [3.8249 3.8249 9.1362 9.1362]*1e-3 ; 0.61 0.61 0.0034 0.0034]; %c1j is kPa to MPa and c2j is unitless
        alpha1 = 23.4 *pi/180; %deg to rad
        alpha4 = 81.1 *pi/180; %deg to rad
        
        %SMC Parameters
        alphaPS = 60*pi/180; %rad
        lMD = 5e-6; %nm to mm
        LLA0 = 10e-6; %nm to mm
        LS20 = 60e-6; %nm to mm
        AS2 = 1.9e-12; %nm^2 to mm^2 - same as N_CF
        deltam = 14.5e-6; %nm to mm - same units as L_Mmax
        LMmax = 1.1e-3; %um to mm - same  units as delta_m
        dMA0 = 53e-6; %nm to mm
        NCU = 167e6; %um^-2 to mm^-2 - (CUs per SMC in unit area)
        NCF = 244e6; %um^-2 to mm^-2 - same as A_S2
        ltOpt = 1.69; %lambda_theta_opt
        thetaSMC = 20 * pi/180; %deg to rad
        k = [0.1926 0.5 0.4 0.1 0.5 0.1926 0.01]; %s^-1
        
        %Parameters Obtained by Fitting Contraction Data (Table 2)
        beta = 0.14; %(s*MPa)^-1
        EAMp = 0.24e3; %GPa to MPa
        EAM = 0.072e3; %0.3*EAMp - According to the paper
        kMAi = [0.93 0.35]; %1-K_MAr, 2-K_MAz;
        cfo = 8.85;
        ufoOpt = 1.85;
        kqp = 0.91;
        gammar = 0.9;
        DKCL = 0.1e-4; %mm^2/s (Diffusion Rate) - needs to be same length units as h_adv
        lzOpt = 1.54; %in vivo axial stretch
    end
    methods
        %Defomation Tensor and etc.
        function F = F(obj,r)
            if exist('r','var')
                F = diag([obj.lr(r),obj.lt(r),obj.lz]);
            else
                F = diag([obj.lr,obj.lt,obj.lz]);
            end
        end
        function C = C(obj,r)
            if exist('r','var')
                F = obj.F(r);
            else
                F = obj.F;
            end
            C = transpose(F)*F;
        end
        function Ffs = Ffs(obj,r)
            if exist('r','var')
                Ufs = obj.ufsfit(r);
            else
                Ufs = obj.ufs;
            end
            Ffs = Ufs*(obj.MSMC'*obj.MSMC)+(eye(3)-(obj.MSMC'*obj.MSMC))/sqrt(Ufs);
        end
        function Fe = Fe(obj,r)
            if exist('r','var')
                Fe = zeros(3,3,length(r))*obj.ri;
                for j=1:length(r)
                    Fe(:,:,j) = obj.F(r(j))/obj.Ffs(r(j)); %Fe=F*inv(Ffs)
                end
            else
                Fe = obj.F/obj.Ffs; %Fe=F*inv(Ffs)
            end
        end
        function Ce = Ce(obj,r)
            if exist('r','var')
                Fe = obj.Fe(r);
            else
                Fe = obj.Fe;
            end
            Ce = zeros(3,3,length(Fe))*obj.ri; 
            for j=1:length(Fe)
                Ce(:,:,j) = transpose(Fe(:,:,j))*Fe(:,:,j);
            end
        end
        
        %Muscle Vectors
        function MSMC = MSMC(obj)
            MSMC = [sin(obj.thetaSMC), cos(obj.thetaSMC), 0];
        end
        function MSMCr = MSMCr(obj)
            MSMCr = [cos(obj.thetaSMC), -sin(obj.thetaSMC), 0];
        end
        function MSMCz = MSMCz(~)
            MSMCz = [0, 0, 1];
        end
        
        %Collagen Fibers
        function alphaj = alphaj(obj,j)
            alphaj = [obj.alpha1, -obj.alpha1, -obj.alpha4, obj.alpha4]; %deg to rad
            alphaj = alphaj(j);
        end
        function MECMj = MECMj(obj,j)
            MECMj = [0, sin(obj.alphaj(j)), cos(obj.alphaj(j))];
        end
        
        %I4 Functions
        function I4SMCe = I4SMCe(obj,r)
            if exist('r','var')
                lr = obj.lr(r); lt = obj.lt(r);
            else
                lr = obj.lr; lt = obj.lt;
            end
            
            I4SMCe = ((lr.^2).*sin(obj.thetaSMC).^2 + (lt.^2).*cos(obj.thetaSMC).^2)./obj.ufs.^2;
        end
        function I4SMC = I4SMC(obj,r)
            if exist('r','var')
                lr = obj.lr(r); lt = obj.lt(r);
            else
                lr = obj.lr; lt = obj.lt;
            end
            I4SMC = (lr.^2)*sin(obj.thetaSMC)^2 + (lt.^2)*cos(obj.thetaSMC)^2; %from #12
        end
        function I4SMCr = I4SMCr(obj,r)
            if exist('r','var')
                lr = obj.lr(r); lt = obj.lt(r);
            else
                lr = obj.lr; lt = obj.lt;
            end
            I4SMCr = (lr.^2)*cos(obj.thetaSMC)^2 + (lt.^2)*sin(obj.thetaSMC)^2; %from #12
        end
        function I4SMCz = I4SMCz(obj)
            I4SMCz = obj.lz.^2; %from #12
        end
        
        function [dMArMA0,dMAzMA0] = dMAiMA0(obj,r)
            I4SMCz = obj.I4SMCz;
            if exist('r','var')
                I4SMCr = obj.I4SMCr(r);
                
                dMArMA0 = obj.kMAi(1)*(sqrt(I4SMCr)-1) + 1;
                dMAzMA0 = obj.kMAi(2)*(sqrt(I4SMCz)-1) + 1;
            else
                I4SMCr = obj.I4SMCr;

                dMArMA0 = obj.kMAi(1)*(sqrt(I4SMCr)-1) + 1;
                dMAzMA0 = obj.kMAi(2)*(sqrt(I4SMCz)-1) + 1;
            end
        end
        function [LMr,LMz] = LMi(obj,r)
            if exist('r','var')
                [dMArMA0,dMAzMA0] = obj.dMAiMA0(r);
                LMr = obj.LMmax*(1-obj.kqp*dMArMA0); LMr(subs(LMr,obj.riNum)<0) = 0;
                LMz = obj.LMmax*(1-obj.kqp*dMAzMA0); LMz(subs(LMz,obj.riNum)<0) = 0;
            else
                [dMArMA0,dMAzMA0] = obj.dMAiMA0;
                LMr = obj.LMmax*(1-obj.kqp*dMArMA0); LMr(subs(LMr,obj.lrNum)<0) = 0;
                LMz = obj.LMmax*(1-obj.kqp*dMAzMA0); LMz(subs(LMz,obj.lrNum)<0) = 0;
            end
        end
        function [Lfor,Lfoz] = Lfoi(obj,r)
            if exist('r','var')
                [LMr,LMz] = obj.LMi(r);
                Ufs = obj.ufsfit(r);
            else
                [LMr,LMz] = obj.LMi;
                Ufs = obj.ufs;
            end
            Lfor = exp( -power(Ufs-obj.ufoOpt,2)./(2*power(obj.cfo.*LMr./obj.LMmax,2)) );
            Lfoz = exp( -power(Ufs-obj.ufoOpt,2)./(2*power(obj.cfo.*LMz./obj.LMmax,2)) );
        end
        function [eS2x,eS2y] = eS2(obj,i,r)
            if exist('r','var')
                [dMArMA0,dMAzMA0] = obj.dMAiMA0(r);
                dMArMA0 = subs(dMArMA0,obj.riNum);
                dMAzMA0 = subs(dMAzMA0,obj.riNum);
            else
                [dMArMA0,dMAzMA0] = obj.dMAiMA0;
                dMArMA0 = subs(dMArMA0,obj.lrNum);
                dMAzMA0 = subs(dMAzMA0,obj.lrNum);
            end
            
            if i==1
                dMAiMA0 = dMArMA0;
            else
                dMAiMA0 = dMAzMA0;
            end
            
            dMAi = dMAiMA0*obj.dMA0;
            
            eS2x = zeros(1,length(dMAi)); eS2y = eS2x;
            for j=1:length(dMAi)
                if dMAiMA0(j) < 1
                    LLAi = obj.LLA0.*dMAiMA0(j);
                    LS2i = obj.LS20 + obj.LLA0.*(1-dMAiMA0(j));
                else
                    LLAi = obj.LLA0;
                    LS2i = obj.LS20;
                end

                if (obj.lMD + obj.LLA0 + obj.LS20) > dMAi(j) && dMAi(j) >= obj.lMD
                    x2i = sqrt(LS2i.^2 - (dMAi(j)-obj.lMD-LLAi).^2);
                    y2i = dMAi(j) - obj.lMD - LLAi;

                    x2ci = x2i + LLAi.*sin(obj.alphaPS);
                    y2ci = dMAi(j) - obj.lMD - LLAi.*cos(obj.alphaPS);

                    eS2x(1,j) = (x2ci./sqrt(x2i.^2 + y2i.^2))*(1 - sqrt((x2i.^2 + y2i.^2)/(x2ci.^2 + y2ci.^2)));
                    eS2y(1,j) = (y2ci./sqrt(x2i.^2 + y2i.^2))*(1 - sqrt((x2i.^2 + y2i.^2)/(x2ci.^2 + y2ci.^2)));
                else
                    eS2x(1,j) = 0;
                    eS2y(1,j) = 0;
                end
            end
        end
        function obj = ufsUpdate(obj,LMr,Lfor,LMz,Lfoz,eS2xr,eS2xz,I4SMCe,lr,lt,Ufs)
            if ~exist('LMr','var')
                [LMr,LMz] = obj.LMi;
                [Lfor,Lfoz] = obj.Lfoi;
                [eS2xr,~] = obj.eS2(1);
                [eS2xz,~] = obj.eS2(2);
                I4SMCe = subs(obj.I4SMCe,obj.lrNum);
                lr = obj.lrNum;
                lt = obj.lt;
                Ufs = obj.ufs;
            end
            if length(LMr)>1
                NumVal = obj.riNum;
            else
                NumVal = obj.lrNum;
            end
            
            ufsdot = double(obj.beta*subs(obj.PCU(LMr,Lfor,LMz,Lfoz,I4SMCe,lr,lt,Ufs)-obj.PMM(LMr,Lfor,LMz,Lfoz,eS2xr,eS2xz),NumVal));
            if ufsdot<0
                Ufs = Ufs + ufsdot*obj.dt;
                obj.ufs = Ufs;
            end
        end
        
        %Stress Functions(:,1)-r (:,2)-theta (:,3)-z
        function sECM = sECM(obj,r)
            if exist('r','var')
                lr = obj.lr(r);
                lt = obj.lt(r);
            else
                lr = obj.lr;
                lt = obj.lt;
            end
            aj = obj.alphaj(1:4);
            C1 = (lt.^2).*sin(aj).^2 + (obj.lz.^2)*cos(aj).^2 - 1;
            
            sECM = [...
                obj.Cp.*lr.^2;
                lt.^2.*(obj.Cp + sum( obj.c(1,:).*exp( obj.c(2,:).*C1.^2).*(sin(aj).^2).*C1 ) );
                obj.lz.^2.*(obj.Cp + sum( obj.c(1,:).*exp( obj.c(2,:).*C1.^2).*(cos(aj).^2).*C1 ) );
                ];
        end
        function sSMC = sSMC(obj,I4SMCe,LMr,Lfor,LMz,Lfoz,lr,lt,Ufs)
            if ~exist('I4SMCe','var')
                I4SMCe = obj.I4SMCe;
                [LMr,LMz] = obj.LMi;
                [Lfor,Lfoz] = obj.Lfoi;
                lr = obj.lr;
                lt = obj.lt;
                Ufs = obj.ufs;
            end
            
            sSMC = ...
            	[(obj.nAMp*obj.EAMp + obj.nAM*obj.EAM)*obj.AS2*obj.NCF*(sqrt(I4SMCe)-1).*...
                    ( LMr.*Lfor.*obj.gammar + LMz.*Lfoz.*(1-obj.gammar) ).*...
                    (lr.^2).*power(sin(obj.thetaSMC),2)./...
                    ( 2.*sqrt(I4SMCe).*obj.deltam.*power(Ufs,2) ) ;
                (obj.nAMp.*obj.EAMp + obj.nAM.*obj.EAM).*obj.AS2.*obj.NCF.*(sqrt(I4SMCe)-1).*...
                    ( LMr.*Lfor.*obj.gammar + LMz.*Lfoz.*(1-obj.gammar) ).*...
                    (lt.^2).*power(cos(obj.thetaSMC),2)./...
                    ( 2.*sqrt(I4SMCe).*obj.deltam.*power(Ufs,2) ) ;
                zeros(1,length(LMr))];
        end
        function sMMy = sMMy(obj,LMr,Lfor,LMz,Lfoz,eS2yr,eS2yz,lr,lt)
            if ~exist('LMr','var')
                [LMr,LMz] = obj.LMi;
                [Lfor,Lfoz] = obj.Lfoi;
                [~,eS2yr] = obj.eS2(1);
                [~,eS2yz] = obj.eS2(2);
                lr = obj.lr;
                lt = obj.lt;
            end
        	sMMy = (obj.nAMp.*obj.EAMp + obj.nAM.*obj.EAM).*obj.AS2.*obj.NCU.*...
                [LMr.*Lfor.*obj.gammar.*eS2yr.*lr.*power(cos(obj.thetaSMC),2) ;
                 LMr.*Lfor.*obj.gammar.*eS2yr.*lt.*power(sin(obj.thetaSMC),2) ;
                 LMz.*Lfoz.*(1-obj.gammar).*eS2yz.*obj.lz ] ./obj.deltam;
        end
        
        function PMM = PMM(obj,LMr,Lfor,LMz,Lfoz,eS2xr,eS2xz)
            if ~exist('LMr','var')
                [LMr,LMz] = obj.LMi;
                [Lfor,Lfoz] = obj.Lfoi;
                [eS2xr,~] = obj.eS2(1);
                [eS2xz,~] = obj.eS2(2);
            end
            PMM = ((obj.nAMp.*obj.EAMp./obj.deltam).*obj.AS2.*obj.NCF)*( obj.gammar.*Lfor.*LMr.*eS2xr + (1-obj.gammar).*Lfoz.*LMz.*eS2xz );
        end
        function PCU = PCU(obj,LMr,Lfor,LMz,Lfoz,I4SMCe,lr,lt,Ufs)
            if ~exist('LMr','var')
                [LMr,LMz] = obj.LMi;
                [Lfor,Lfoz] = obj.Lfoi;
                I4SMCe = subs(obj.I4SMCe,obj.lrNum);
                lr = obj.lrNum;
                lt = obj.lt;
                Ufs = obj.ufs;
            end
            PCU = (obj.nAMp.*obj.EAMp+obj.nAM.*obj.EAM).*obj.AS2.*obj.NCF.*(sqrt(I4SMCe)-1).*...
                    ( LMr.*Lfor.*obj.gammar + LMz.*Lfoz.*(1-obj.gammar) ).*...
                    ( lr.*(sin(obj.thetaSMC).^2) + lt.*(cos(obj.thetaSMC).^2) )./...
                    ( 2.*sqrt(I4SMCe).*obj.deltam.*Ufs.^2 );
        end
    end
end
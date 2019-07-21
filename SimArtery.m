classdef SimArtery < handle
    properties
        
        ri
        ufs
        
        nAMp
        nAM
        
        dMArMA0
        dMAzMA0
        LMr
        LMz
        Lfor
        Lfoz
        eS2xr
        eS2yr
        eS2xz
        eS2yz
        
        sECM
        sSMC
        sMMy
        
        PMM
        PCU
        
        riVec
        ufsVec
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
            lr = obj.lr(r);
            lt = obj.lt(r);
            F = sym(zeros(3,3,r));
            for j=1:length(r)
                F = diag([lr(j),lt(j),lz]);
            end
        end
        function C = C(obj,r)
            F = obj.F(r);
            C = sym(zeros(3,3,length(r)));
            for j=1:length(r)
                C(:,:,j) = transpose(F(:,:,j)*F(:,:,j));
            end
        end
        function Ffs = Ffs(obj,r)
            ufsf = ufsFit(r);
            
            Ffs = sym(zeros(3,3,length(ufsf)));
            Mat = obj.MSMC'*obj.MSMC;
            for j=1:length(ufsf)
                Ffs(:,:,j) = ufsf(j)*Mat+(eye(3)-Mat)/sqrt(ufsf(j));
            end
        end
        function Fe = Fe(obj,r)
        	F = obj.F(r);
            Ffs = obj.Ffs(r);
            Fe = sym(zeros(3,3,length(r)));
            for j=1:length(r)
            	Fe(:,:,j) = F(:,:,j)/Ffs(:,:,j); %Fe=F*inv(Ffs)
            end
        end
        function Ce = Ce(obj,r)
            Fe = obj.Fe(r);
            Ce = sym(zeros(3,3,length(r)));
            for j=1:length(r)
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
            I4SMCe = ((obj.lr(r).^2).*sin(obj.thetaSMC).^2 + (obj.lt(r).^2).*cos(obj.thetaSMC).^2)./obj.ufsFit(r).^2;
        end
        function I4SMC = I4SMC(obj,r)
            I4SMC = (obj.lr(r).^2)*sin(obj.thetaSMC)^2 + (obj.lt(r).^2)*cos(obj.thetaSMC)^2; %from #12
        end
        function I4SMCr = I4SMCr(obj,r)
            I4SMCr = (obj.lr(r).^2)*cos(obj.thetaSMC)^2 + (obj.lt(r).^2)*sin(obj.thetaSMC)^2; %from #12
        end
        function I4SMCz = I4SMCz(obj)
            I4SMCz = obj.lz.^2; %from #12
        end
        
        function obj = dMAiMA0(obj,r)
            obj.dMArMA0 = obj.kMAi(1)*(sqrt(obj.I4SMCr(r))-1) + 1;
            obj.dMAzMA0 = obj.kMAi(2)*(sqrt(obj.I4SMCz)-1) + 1;
        end
        function obj = LMi(obj,r)
            obj.dMAiMA0(r);
        	obj.LMr = obj.LMmax*(1-obj.kqp*obj.dMArMA0); obj.LMr(obj.LMr<0) = 0;
        	obj.LMz = obj.LMmax*(1-obj.kqp*obj.dMAzMA0); obj.LMz(obj.LMz<0) = 0;
        end
        function obj = Lfoi(obj,r)
            obj.LMi(r);
            obj.Lfor = exp( -power(obj.ufsFit(r)-obj.ufoOpt,2)./(2*power(obj.cfo.*obj.LMr./obj.LMmax,2)) );
            obj.Lfoz = exp( -power(obj.ufsFit(r)-obj.ufoOpt,2)./(2*power(obj.cfo.*obj.LMz./obj.LMmax,2)) );
        end
        function obj = eS2(obj)
            for i=1:2
                %i=1 for r direction and i=2 for z
                if i==1
                    dMAiMA0 = obj.dMArMA0;
                else
                    dMAiMA0 = obj.dMAzMA0;
                end
                dMAi = dMAiMA0*obj.dMA0;
                
                x2cit = zeros(1,length(dMAi)); y2cit = x2cit;
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
                        
                        x2cit(1,j) = x2ci; y2cit(1,j) = y2ci;
                        
                        eS2x(1,j) = (x2ci./sqrt(x2i.^2 + y2i.^2))*(1 - sqrt((x2i.^2 + y2i.^2)/(x2ci.^2 + y2ci.^2)));
                        eS2y(1,j) = (y2ci./sqrt(x2i.^2 + y2i.^2))*(1 - sqrt((x2i.^2 + y2i.^2)/(x2ci.^2 + y2ci.^2)));
                    else
                        eS2x(1,j) = 0;
                        eS2y(1,j) = 0;
                    end
                end
                if i==1
                    obj.eS2xr = eS2x;
                    obj.eS2yr = eS2y;
                else
                    obj.eS2xz = eS2x;
                    obj.eS2yz = eS2y;
                end
            end
        end

        %Stress Functions(:,1)-r (:,2)-theta (:,3)-z
        function obj = sECMCalc(obj,r)
            aj = obj.alphaj(1:4);
            
            lr = obj.lr(r);
            lt = obj.lt(r);
            
            obj.sECM = sym(zeros(3,length(r)));
            for j=1:length(r)
                C1 = (lt(j)^2).*sin(aj).^2 + (obj.lz^2)*cos(aj).^2 - 1;
                
                obj.sECM(:,j) = [...
                    obj.Cp.*lr(j).^2;
                    lt(j)^2.*(obj.Cp + sum( obj.c(1,:).*exp( obj.c(2,:).*C1.^2).*(sin(aj).^2).*C1 ) );
                    obj.lz^2.*(obj.Cp + sum( obj.c(1,:).*exp( obj.c(2,:).*C1.^2).*(cos(aj).^2).*C1 ) );
                    ];
            end
        end
        function obj = sSMCCalc(obj,r)
            obj.sSMC = ...
            	[(obj.nAMp*obj.EAMp + obj.nAM*obj.EAM)*obj.AS2*obj.NCF*(sqrt(obj.I4SMCe(r))-1).*...
                    ( obj.LMr.*obj.Lfor.*obj.gammar + obj.LMz.*obj.Lfoz.*(1-obj.gammar) ).*...
                    (obj.lr(r).^2).*power(sin(obj.thetaSMC),2)./...
                    ( 2.*sqrt(obj.I4SMCe(r)).*obj.deltam.*power(obj.ufsFit(r),2) ) ;
                (obj.nAMp.*obj.EAMp + obj.nAM.*obj.EAM).*obj.AS2.*obj.NCF.*(sqrt(obj.I4SMCe(r))-1).*...
                    ( obj.LMr.*obj.Lfor.*obj.gammar + obj.LMz.*obj.Lfoz.*(1-obj.gammar) ).*...
                    (obj.lt(r).^2).*power(cos(obj.thetaSMC),2)./...
                    ( 2.*sqrt(obj.I4SMCe(r)).*obj.deltam.*power(obj.ufsFit(r),2) ) ;
                zeros(1,length(obj.LMr))];
        end
        function obj = sMMyCalc(obj,r)
        	obj.sMMy = ...
                (obj.nAMp.*obj.EAMp + obj.nAM.*obj.EAM).*obj.AS2.*obj.NCU.*...
                [obj.LMr.*obj.Lfor.*obj.gammar.*obj.eS2yr.*obj.lr(r).*power(cos(obj.thetaSMC),2) ;
                 obj.LMr.*obj.Lfor.*obj.gammar.*obj.eS2yr.*obj.lt(r).*power(sin(obj.thetaSMC),2) ;
                 obj.LMz.*obj.Lfoz.*(1-obj.gammar).*obj.eS2yz.*obj.lz ] ./obj.deltam;
        end
        
        function PMM = PMMCalc(obj)
            PMM = ((obj.nAMp.*obj.EAMp./obj.deltam).*obj.AS2.*obj.NCF)*( obj.gammar.*obj.Lfor.*obj.LMr.*obj.eS2xr + (1-obj.gammar).*obj.Lfoz.*obj.LMz.*obj.eS2xz );
            obj.PMM = PMM;
        end
        function PCU = PCUCalc(obj)
            PCU = (obj.nAMp.*obj.EAMp+obj.nAM.*obj.EAM).*obj.AS2.*obj.NCF.*(sqrt(obj.I4SMCe(r))-1).*...
                    ( obj.LMr.*obj.Lfor.*obj.gammar + obj.LMz.*obj.Lfoz.*(1-obj.gammar) ).*...
                    ( obj.lr(r).*(sin(obj.thetaSMC).^2) + obj.lt(r).*(cos(obj.thetaSMC).^2) )./...
                    ( 2.*sqrt(obj.I4SMCe(r)).*obj.deltam.*obj.ufsFit(r).^2 );
            obj.PCU = PCU;
        end
        
        %some shortcuts
        function ufs = ufsFit(obj,r)
            ufs = interp1(linspace(obj.ri,obj.ro,length(obj.ufs)),obj.ufs,r);
        end
    end
end
classdef Artery < handle
    properties
        cs = CurrentState;
        
        %Vectors
        V = SimulationVectors;
                
        %Stress-free Parameters
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
        %EAM = 0.072e3; %0.3*EAMp - According to the paper
        kMAi = [0.93 0.35]; %1-K_MAr, 2-K_MAz;
        cfo = 8.85;
        ufoOpt = 1.85;
        kqp = 0.91;
        gammar = 0.9;
        DKCL = 0.1e-4; %mm^2/s (Diffusion Rate) - needs to be same length units as h_adv
        lzOpt = 1.54; %in vivo axial stretch
    end
    methods
        function EAM = EAM(obj)
            EAM = 0.3*obj.EAMp;
        end
        
        %Defomation Tensor and etc.
        function F = F(obj)
            F = sym(zeros(3,3,length(obj.cs.lr)));
            for j=1:length(obj.cs.lr)
                F(:,:,j) = diag([obj.cs.lr(j),obj.cs.lt(j),obj.cs.lz]);
            end
        end
        function C = C(obj)
            F = obj.F;
            C = sym(zeros(3,3,length(obj.cs.lr)));
            for j=1:length(obj.cs.lr)
                C(:,:,j) = transpose(F(:,:,j)*F(:,:,j));
            end
        end
        function Ffs = Ffs(obj)
            Ffs = sym(zeros(3,3,length(obj.cs.ufs)));
            Mat = obj.MSMC'*obj.MSMC;
            for j=1:length(obj.cs.ufs)
                Ffs(:,:,j) = obj.cs.ufs(j)*Mat+(eye(3)-Mat)/sqrt(obj.cs.ufs(j));
            end
        end
        function Fe = Fe(obj)
        	F = obj.F;
            Ffs = obj.Ffs;
            Fe = sym(zeros(3,3,length(obj.cs.ufs)));
            for j=1:length(obj.cs.ufs)
            	Fe(:,:,j) = F(:,:,j)/Ffs(:,:,j); %Fe=F*inv(Ffs)
            end
        end
        function Ce = Ce(obj)
            Fe = obj.Fe;
            Ce = sym(zeros(3,3,length(obj.cs.ufs)));
            for j=1:length(obj.cs.ufs)
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
        function I4SMCe = I4SMCe(obj)
            I4SMCe = ((obj.cs.lr.^2).*sin(obj.thetaSMC).^2 + (obj.cs.lt.^2).*cos(obj.thetaSMC).^2)./obj.cs.ufs.^2;
        end
        function I4SMC = I4SMC(obj)
            I4SMC = (obj.cs.lr.^2)*sin(obj.thetaSMC)^2 + (obj.cs.lt.^2)*cos(obj.thetaSMC)^2; %from #12
        end
        function I4SMCr = I4SMCr(obj)
            I4SMCr = (obj.cs.lr.^2)*cos(obj.thetaSMC)^2 + (obj.cs.lt.^2)*sin(obj.thetaSMC)^2; %from #12
        end
        function I4SMCz = I4SMCz(obj)
            I4SMCz = obj.cs.lz.^2; %from #12
        end
        function I4f = I4f(obj,a)
            I4f = (obj.cs.lt.^2).*(sin(a)^2) + (obj.cs.lz.^2)*(cos(a)^2);
        end
        
        function obj = dMAiMA0(obj)
            obj.cs.dMArMA0 = obj.kMAi(1)*(sqrt(obj.I4SMCr)-1) + 1;
            obj.cs.dMAzMA0 = obj.kMAi(2)*(sqrt(obj.I4SMCz)-1) + 1;
        end
        function obj = LMi(obj)
            obj.dMAiMA0;
        	obj.cs.LMr = obj.LMmax*(1-obj.kqp*obj.cs.dMArMA0); 
        	obj.cs.LMz = obj.LMmax*(1-obj.kqp*obj.cs.dMAzMA0); obj.cs.LMz(obj.cs.LMz<0) = 0;
            for i=1:length(obj.cs.lrNum)
                if subs(obj.cs.LMr(i),obj.cs.lrNum(i))<0
                    obj.cs.LMr(i) = 0;
                end
            end
        end
        function obj = Lfoi(obj)
            obj.LMi;
            obj.cs.Lfor = exp( -power(obj.cs.ufs-obj.ufoOpt,2)./(2*power(obj.cfo.*obj.cs.LMr./obj.LMmax,2)) );
            obj.cs.Lfoz = exp( -power(obj.cs.ufs-obj.ufoOpt,2)./(2*power(obj.cfo.*obj.cs.LMz./obj.LMmax,2)) );
        end
        function obj = eS2(obj)
            for i=1:2
                %i=1 for r direction and i=2 for z
                if i==1
                    dMAiMA0 = obj.cs.dMArMA0;
                else
                    dMAiMA0 = obj.cs.dMAzMA0;
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
                        y2i = dMAi(j) - obj.lMD - LLAi;
                        x2i = sqrt(LS2i.^2 - y2i.^2);
                        
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
                    obj.cs.eS2xr = eS2x;
                    obj.cs.eS2yr = eS2y;
                    
                    obj.cs.x2cr = x2cit;
                    obj.cs.y2cr = y2cit;
                else
                    obj.cs.eS2xz = eS2x;
                    obj.cs.eS2yz = eS2y;
                    
                    obj.cs.x2cz = x2cit;
                    obj.cs.y2cz = y2cit;
                end
            end
        end
        function obj = ufs0(obj)
            %calcs numeric ufs values by I4SMCe=1
            obj.cs.ufs = sqrt( (obj.cs.lrNum.^2)*power(sin(obj.thetaSMC),2) + (obj.cs.ltNum.^2)*power(cos(obj.thetaSMC),2) );
        end
        function obj = ufsUpdate(obj)
            strain = vertcat(obj.cs.lr,obj.cs.lt,sym(ones(1,length(obj.cs.lt))*obj.lz));
            obj.cs.lr = obj.cs.lrNum;
            obj.cs.lt = obj.cs.ltNum;
            obj.cs.lz = obj.cs.lzNum;
            
            obj.Lfoi;
            obj.eS2;
            
            ufsdot = double(obj.beta*(obj.PCU-obj.PMM));
            if ufsdot<0
                obj.cs.ufs = obj.cs.ufs + ufsdot*obj.dt; %#ok<*MCNPN>
            end
            %Save numeric I4SMCe value for uniaxial Pisom calc
            obj.cs.I4SMCeNum = obj.I4SMCe;
            
            obj.cs.lr = strain(1,:);
            obj.cs.lt = strain(2,:);
            obj.cs.lz = strain(3,:);
            
            %obj.Lfoi; %To use symbolic values of these variables
        end
        
        %Stress Functions(:,1)-r (:,2)-theta (:,3)-z
        function obj = sECM(obj)
            aj = obj.alphaj(1:4);
            obj.cs.sECM = sym(zeros(3,length(obj.cs.lr)));
            z = 1;
            for j=1:length(obj.cs.lr)
                
                C1 = (obj.cs.lt(j)^2).*sin(aj).^2 + (obj.cs.lz(z)^2)*cos(aj).^2 - 1;
                
                %Calc numeric value of C1 to determine I4f=C1+1 
                C1Num = subs(C1,'ri',obj.cs.riNum);
                C1(C1Num<0) = 0;
%                 if ~isempty(C1(C1==0))
%                     vpa(C1Num,2)
%                 end
                
                obj.cs.sECM(:,j) = [...
                    obj.Cp.*obj.cs.lr(j).^2;
                    obj.cs.lt(j)^2.*(obj.Cp + sum( obj.c(1,:).*exp(obj.c(2,:).*C1.^2).*(sin(aj).^2).*C1 ) );
                    obj.cs.lz(z)^2.*(obj.Cp + sum( obj.c(1,:).*exp(obj.c(2,:).*C1.^2).*(cos(aj).^2).*C1 ) );
                    ];
                if length(obj.cs.lz) > 1
                    z = z+1;
                end
            end
        end
        function obj = sSMC(obj)
            obj.cs.sSMC = ...
            	[(obj.cs.nAMp*obj.EAMp + obj.cs.nAM*obj.EAM)*obj.AS2*obj.NCF*(sqrt(obj.I4SMCe)-1).*...
                    ( obj.cs.LMr.*obj.cs.Lfor.*obj.gammar + obj.cs.LMz.*obj.cs.Lfoz.*(1-obj.gammar) ).*...
                    (obj.cs.lr.^2).*power(sin(obj.thetaSMC),2)./...
                    ( 2.*sqrt(obj.I4SMCe).*obj.deltam.*power(obj.cs.ufs,2) ) ;
                (obj.cs.nAMp.*obj.EAMp + obj.cs.nAM.*obj.EAM).*obj.AS2.*obj.NCF.*(sqrt(obj.I4SMCe)-1).*...
                    ( obj.cs.LMr.*obj.cs.Lfor.*obj.gammar + obj.cs.LMz.*obj.cs.Lfoz.*(1-obj.gammar) ).*...
                    (obj.cs.lt.^2).*power(cos(obj.thetaSMC),2)./...
                    ( 2.*sqrt(obj.I4SMCe).*obj.deltam.*power(obj.cs.ufs,2) ) ;
                zeros(1,length(obj.cs.LMr))];
        end
        function obj = sMMy(obj)
        	obj.cs.sMMy = ...
                (obj.cs.nAMp.*obj.EAMp + obj.cs.nAM.*obj.EAM).*obj.AS2.*obj.NCU.*...
                [obj.cs.LMr.*obj.cs.Lfor.*obj.gammar.*obj.cs.eS2yr.*obj.cs.lr.*power(cos(obj.thetaSMC),2) ;
                 obj.cs.LMr.*obj.cs.Lfor.*obj.gammar.*obj.cs.eS2yr.*obj.cs.lt.*power(sin(obj.thetaSMC),2) ;
                 obj.cs.LMz.*obj.cs.Lfoz.*(1-obj.gammar).*obj.cs.eS2yz.*obj.cs.lz ] ./obj.deltam;
        end
        
        function PMM = PMM(obj)
            PMM = ((obj.cs.nAMp.*obj.EAMp./obj.deltam).*obj.AS2.*obj.NCF)*( obj.gammar.*obj.cs.Lfor.*obj.cs.LMr.*obj.cs.eS2xr + (1-obj.gammar).*obj.cs.Lfoz.*obj.cs.LMz.*obj.cs.eS2xz );
            obj.cs.PMM = PMM;
        end
        function PCU = PCU(obj)
            PCU = (obj.cs.nAMp.*obj.EAMp+obj.cs.nAM.*obj.EAM).*obj.AS2.*obj.NCF.*(sqrt(obj.I4SMCe)-1).*...
                    ( obj.cs.LMr.*obj.cs.Lfor.*obj.gammar + obj.cs.LMz.*obj.cs.Lfoz.*(1-obj.gammar) ).*...
                    ( obj.cs.lr.*(sin(obj.thetaSMC).^2) + obj.cs.lt.*(cos(obj.thetaSMC).^2) )./...
                    ( 2.*sqrt(obj.I4SMCe).*obj.deltam.*obj.cs.ufs.^2 );
            obj.cs.PCU = PCU;
        end
    end
end
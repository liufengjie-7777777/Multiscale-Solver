classdef SimArteryVessel < SimArtery
    properties
        dt = 1; %Seconds
        TotalTime = 15; %Minutes
        
        Pin = 90*133.322387415*1e-6; %Inner pressure (MPa)
        lambda = 1.5;
    end
    methods
        %Length Dimensions Functions 
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
        %Calculate Axial Force for a Given Stress
        function FT = FTCalc(obj)
            [x,w] = obj.xw(obj.ro);
            
            obj.Lfoi(x);
            obj.eS2;
            
            obj.sECMCalc(x);
            obj.sSMCCalc(x);
            obj.sMMyCalc(x);
            S = obj.sECM + obj.sSMC + obj.sMMy;
            Sr = S(1,:);
            St = S(2,:);
            Sz = S(3,:);
            
            f = subs(2*Sz-St-Sr,obj.ri);
            FT = double( pi*sum((f.*x).*w) );
        end
        
        %Calculate Lagrange Multiplier
        function p = pCalc(obj,r)
            obj.Lfoi(r);
            obj.eS2;
            
            obj.sECMCalc(r);
            obj.sSMCCalc(r);
            obj.sMMyCalc(r);
            s = obj.sECM + obj.sSMC + obj.sMMy;
            Srbar = s(1);
            
            if r > obj.ri
                [x,w] = obj.xw(r);
                obj.Lfoi(x);
                obj.eS2;

                obj.sECMCalc(x);
                obj.sSMCCalc(x);
                obj.sMMyCalc(x);

                s = obj.sECM + obj.sSMC + obj.sMMy;
                Sr = s(1,:);
                St = s(2,:);
                p = double(obj.Pin + Srbar - sum( (St-Sr)./x.*w));
            else
                p = double(obj.Pin + Srbar);
            end
        end
        
        function [p,S] = CauchyStress(obj,r)
            p = obj.pCalc(r);
            
            obj.Lfoi(r);
            obj.eS2;
            
            obj.sECMCalc(r);
            obj.sSMCCalc(r);
            obj.sMMyCalc(r);
            
            Sbar = double(obj.sECM + obj.sSMC + obj.sMMy);
            
            S = -p + Sbar;
        end
        
        %Calculate Legendre Points and Weights
        function [x,w] = xw(obj,r)
            [x,w] = lgwt(3,obj.ri,r);
        end
        
    end
end

classdef CurrentState < handle
    properties
        %Uniaxial
        lr = sym('lr','positive');
        lt = 1.69;
        
        %Biaxial
        ri = sym('ri','positive'); %Inner Radius
        Pin = 90*133.322387415*1e-6; %Inner Pressure
        lambda = 1.5;
        
        riNum = 0;
        FT
        
        %Both Simulations
        p = sym('p','real');
        
        ufs = 0;
        nAMp = 0;
        nAM = 0;
        
        lrNum
        ltNum
        lz
        
        LMr
        LMz
        Lfor
        Lfoz
        eS2xr
        eS2yr
        eS2xz
        eS2yz
        
        I4SMCe
        I4SMCeNum
        
        I4SMCNum
        I4SMCrNum
        
        sECM
        sSMC
        sMMy
        
    end
    
    methods
        
    end
end


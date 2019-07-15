classdef CurrentState < handle
    properties
        %Uniaxial
        lrG = sym('lr','positive');
        ltG = 1.69;
        
        %Biaxial
        Pin = 90*133.322387415*1e-6; %Inner Pressure
        lambda = 1.5;
        
        riG = sym('ri','positive'); %Inner Radius
        
        r
        R
        
        x
        w
        xNum
        wNum
        
        ri
        ro
        riNum = 0;
        roNum
        FT
        
        %Both Simulations
        p = sym('p','real');
        
        ufs
        nAMp
        nAM
        
        %Stretch Ratios
        lr
        lt
        lz
        
        lrNum
        ltNum
        lzNum
        
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
        
        I4SMCe
        I4SMCeNum
        
        sECM
        sSMC
        sMMy
        
        PMM
        PCU
        Pisom
    end
    
    methods
        
    end
end


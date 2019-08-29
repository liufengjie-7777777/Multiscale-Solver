function [sECM] = sECMCalc(aj,c1,c2,Cp,lr,lt,lz)
    sECM = sym(zeros(3,1));
    C1 = (lt^2).*cos(aj).^2 + (lz^2)*sin(aj).^2 - 1;
    
    sECM = [...
        Cp.*lr.^2;
        lt^2.*(Cp + sum( c1.*exp( c2.*C1.^2).*(cos(aj).^2).*2.*C1 ) );
        lz^2.*(Cp + sum( c1.*exp( c2.*C1.^2).*(sin(aj).^2).*2.*C1 ) );
        ];
end

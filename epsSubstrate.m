function [ epsOut ] = epsSubstrate( lamArr)
%EPSSUBSTRATE relative permittivity of substrate
% lamArr = operating wavelnegth, um


A = 6.255;
B= 2.316;
C = 0.6263^2;
D = 2.765;
E = 32.935^2;

epsOut = 1 + A + (B*lamArr.^2)./(lamArr.^2-C)+(D*lamArr.^2)./(lamArr.^2-E);

end


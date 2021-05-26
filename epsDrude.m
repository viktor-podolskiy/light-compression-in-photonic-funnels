function [ epsOut] = epsDrude( lamArr)
%EPSDRUDE returns relative permittivity of the doped semiconductor-based desingner metal
% lamArr = operating wavelnegth, um

c0=3e8; % speed of light in vacuum
omgArr=2e6*pi*c0./lamArr; % angular frequency

tau=1e13; %inelastic scattering frequency, 1/s
lamPl=8.5; % plasma wavelength, um
eps0=12.15; % permittivity of undoped semiconductor 

omgPl=2e6*pi*c0/lamPl; % plasma frequency 

epsOut=eps0*(1-omgPl^2./omgArr./(omgArr+1i*tau)); 


end


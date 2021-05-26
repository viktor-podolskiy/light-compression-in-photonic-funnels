function [rAvg] = rFrac(EEArrLam, r2,z2, rMax,dr,z0,frac)
%REXTENT Calculates the intensity confinement radius
%
% r2,z2 = meshgrid of radial and z coordinates
% EEArrLam = grid of |E| on the meshgrid above
% z0 defines the z coordinate where the confinement scale is calculated
% dr,rMax define the precision of calculation of confimenet radius
% frac represents the fractional cut-off intensity

% calculate the normalized radial intensity profile at a given z point
rArr=(0:dr:rMax); 
z0Arr=z0+0*rArr; 
I0Arr=smooth(interp2(r2,z2,abs(EEArrLam).^2,rArr,z0Arr,'linear',0)); 
I0Arr=I0Arr./max(I0Arr); 

for i0=length(I0Arr):-1:1
    if I0Arr(i0)>frac
        break
    end 
end 
rAvg=rArr(i0); 

end


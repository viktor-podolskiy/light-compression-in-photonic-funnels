clear

% --------
% (c) Viktor Podolskiy & Evan Simmons, U Mass Lowell
% --------
% part of the photonic funnels project - calculates propagation of light
% through conical arrays with composite layered cores
% note: model requires COMSOL/matlab livelink and uses ~55Gb of RAM with
% default settings
%
% see Adv. Opt. Mater. v.8, 2001321 (2020); developed with partial support
% from the National Science Foundation, grant # 2004298
% --------

lamArr=(4:0.25:12); %wavelength range for calculations

% material parameters
epsMArr=(conj(epsDrude(lamArr))); % permittivity of doped layers
epsD=10.23+0*epsMArr; % permittivity of undoped layers
epsZZArr=2./(1./epsD+1./epsMArr); % effective medium response of the composite
epsXYArr=(epsD+epsMArr)/2; 


% load comsol model 
% (model is built in COMSOL Multiphysics 5.5, uses axial symmetry mode of optics module)
model=mphload('funnelsAxial06.mph'); 
% ModelUtil.showProgress(true); 


%--- curved geometry config
xFunTop=0.25; % top radius of the funnel, um
xFunBot=5.; % bottom radius of the funnel, um 
hFun=4.0001; % height of the funnel, um 
hAu=3.2; % height of PEC (gold) sidewall, um

hSubFun=0.1; % vertical under-etch; does not affect composite model
funnelDR=0.05; % funnels radial over-etch 

% coordinates of the middle point of Bezier curve
rFit=1; 
zFit=1; 
wFit=1; 

% geometry of the simulation space
xLen=40; 
dr=0.01; 
z0=40; % 15
rArr=(0.01:dr:xLen); 

% mesh data to extract fields profiles from COMSOL
rFun=(0.01:0.025:15); 
zFun=(-5:0.025:10);
[rFun2,zFun2]=meshgrid(rFun,zFun); 

EELst3=zeros(size(rFun2,1),size(rFun2,2),length(lamArr)); % array to be used to store field profiles
tranArr=zeros(1,length(lamArr)); % array to be used to store transmission data

outFname=['./testAxial06.r=',num2str(xFunTop),'.hAu=',num2str(hAu),...
    '.rFit=', num2str(rFit),'.zFit=1.dr=',num2str(funnelDR),'.nMM.mat']; % file to be used to store transmission and fields data 

for il=1:length(lamArr) %iterate over wavelengths 
    
    % set model parameters
    lam0=lamArr(il); 

    model.param.set('geomLen',[num2str(xLen),'[um]']);
    model.param.set('subHt','15[um]');
    model.param.set('superHt','45[um]'); 
    model.param.set('undHt','0.5[um]');
    model.param.set('pmlSz','10[um]');
    
    model.param.set('funnelRtop',[num2str(xFunTop),' [um]']);
    model.param.set('funnelRbot',[num2str(xFunBot), ' [um]']);
    model.param.set('funnelHt',[num2str(hFun), ' [um]']);
    model.param.set('funnelHtB',[num2str(hSubFun), ' [um]']);
    model.param.set('funnelHtT',[num2str(hFun-hAu), ' [um]']);

    model.param.set('nTop','1');
    model.param.set('nSub',num2str(sqrt(epsSubstrate(lam0))));

    model.param.set('epsD',num2str(epsD(il)));
    model.param.set('epsM',num2str(epsMArr(il)));
    
    model.param.set('htD','80[nm]');
    model.param.set('htM','80[nm]');
    model.param.set('funnelDR',[num2str(funnelDR),'[um]']);
    

    model.param.set('epsUndZZ',num2str(epsSubstrate(lam0)));
    model.param.set('epsUndRP',num2str(epsSubstrate(lam0)));
    
    model.param.set('epsFunZZ',num2str(epsZZArr(il)));
    model.param.set('epsFunRP',num2str(epsXYArr(il)));
    model.param.set('epsFun1ZZ',num2str(epsZZArr(il)));
    model.param.set('epsFun1RP',num2str(epsXYArr(il)));
    
    model.param.set('lam0',[num2str(lam0),'[um]']);
    model.param.set('rFit',[num2str(rFit),'[um]']); 
    model.param.set('zFit',[num2str(zFit),'[um]']); 
    model.param.set('wFit',num2str(wFit)); 
    
    model.param.set('h0Max','lam0/45');
    model.param.set('h0Min','h0Max/5');
    model.param.set('h0Growth','1.1'); 
    
    figure(2) % display the geometry
    clf
    mphgeom(model); 
    xlim([0 7e-6])
    ylim([0 7e-6])
    daspect([1 1 1])

   
    % Run the simulation
    model.study('std1').run;
    
    % post-process the results
    rzLst=1e-6*[rArr(:),z0+0*rArr(:)]';
    PzLst=mphinterp(model,'ewfd.Poavz','coord',rzLst,'Dataset','dset2');
    tranArr(il)=sum(PzLst.*rArr)*dr*2e-12*pi; % calculate transmission (arb units)
    
    rzLst=1e-6*[rFun2(:),zFun2(:)]';
    EELst=mphinterp(model,'ewfd.normE','coord',rzLst,'Dataset','dset2');
    EELst=reshape(EELst,size(rFun2)); 
    EELst3(:,:,il)=EELst; % save intensity distribution 
    
    figure(1) % plot intensity distribution
    clf
    surf(rFun2,zFun2,EELst,'EdgeColor','none'); 
    colormap hot
    view(2)
    colorbar
    drawnow 
    

    % save current data 
    save(outFname,'lamArr','tranArr','epsXYArr','epsZZArr','epsMArr', ...
        'xFunTop', 'xFunBot', 'hFun','hAu', 'hSubFun',...
        'rFun2','zFun2','EELst3',...
        'rFit','zFit','wFit','funnelDR','epsD','epsMArr')
        
    disp([num2str(il),'/',num2str(length(lamArr))])
end 

% plot transmission
figure(10)
hold on 
semilogy(lamArr,tranArr)
set(gca,'YScale','log')

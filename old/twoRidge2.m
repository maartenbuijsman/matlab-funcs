function [an,bn,cn,dn,a1,a2,A0,B0,C0,AL,BL,CL]=twoRidge2(gam0,gamL,knL,dop);
% function [an,bn,cn,dn,a1,a2]=twoRidge(gam0,gamL,knL);
%
% an and bn are the amplitudes of the modes radiating from ridge 1 (west
% and east respectively).  cn and dn are the waves radiating from the
% eastern ridge (west and east respectively.
%
% gam0 is the frac of the WKB western ridge height, gam1 is
% the WKB frac of the eastern.  
%
% knL = k_n * x_1 for each mode n that we want in our solution.  This
% allows us to scpecify variable stratification since the only effect is
% that the wavelengths are not integer multiples of each other. 
% 
% a1 is the one-ridge solution
% a2 is the two-ridge.  
%
%


N = length(knL);
nn=N;


m=[1:N]';
m=m-1;
n = [1:N]';

[N,M]=meshgrid(n,m);

for iij=1:2
  if iij==1
    gam = 1-gamL;
  else
    gam = 1-gam0;
  end;
  A = -(-cos(m*pi*gam)*(n'.*sin(n'*pi*gam))+(m.*sin(m*pi*gam))*cos(n'*pi*gam))./(M.^2-N.^2)/pi;

  junk=((-1).^m*(-1).^n');
  B = (N-cos(m*pi*gam)*(n'.*cos(n'*pi*gam))-(m.*sin(m*pi*gam))*sin(n'*pi*gam))./(M.^2-N.^2)/pi;
  %
  del =gam;
    
  for i=1:nn-1;
    A(i+1,i)=(i*pi*(1-gam)-sin(i*pi*del).*cos(i*pi*del))./(2*i*pi);
    B(i+1,i)=-sin(i*pi*del)^2/2/i/pi;
  end;
  
  
  C = sin((m)*pi*gam)./(m)/pi;
  C(1) = -(1-gam);
  
  if iij==1
    AL = A;
    BL = B;
    CL = C;
  else
    A0 = A;
    B0 = B;
    C0 = C;
  end;
end

% Do the decoupled problem

DL = AL+BL;
D0 = A0+B0;
G = pinv(A0)*C0;
F = -pinv(A0)*D0;
% isolated ridge solutions
a1 = pinv(D0)*C0;
a2 = pinv(DL)*CL;

ep = diag(exp(knL*sqrt(-1)));
em = diag(exp(-knL*sqrt(-1)));

LHS = ((DL*ep)*F+AL*em);
RHS = CL-(DL*ep)*G;

LHS = (AL*em-DL*ep*pinv(A0)*D0);
RHS = (CL-DL*ep*G);

bn=pinv(LHS)*RHS;
cn = pinv(A0)*(C0-D0*bn);

LHS = (DL*ep-AL*em*inv(D0)*A0);
RHS = CL-AL*em*inv(D0)*C0;
cn2 = LHS\RHS;
bn2 = inv(D0)*(C0-A0*cn2);
an2 = bn2+cn2;
dn2 = ep*(ep*cn2+em*bn2);
an = an2;bn=bn2;cn=cn2;dn=dn2;



an = bn+cn;
dn = ep*(ep*cn+em*bn);

nargin;

if nargin>3
  
  nn=length(knL)
  [N,M]=meshgrid([1:nn],[0:nn-1]);
  
  UU=cos(M*pi.*N/nn);
  WW=sin(M*pi.*N/nn);
  
  Wb = -WW*(bn);
  Wc = WW*(cn);
  Wa = WW*(an);
  WbL = -WW*(bn).*exp(sqrt(-1)*knL);
  WcL = WW*(cn).*exp(-sqrt(-1)*knL);
  WaL = WW*(an).*exp(-sqrt(-1)*knL);
  
  %% the sides
  CC1=[];CC2=[];
  for j=linspace(0,pi*2,100)
  clear Ua Ud Uc Ub Wa Wb Wc Wd;
  CCstd = 0;
  Nx=400;
  for i=1:Nx;
    Ua(:,Nx-i+1) = UU*(an.*exp(sqrt(-1)*(knL*(-i)/Nx*4+j)));
    Ud(:,i)      = UU*(dn.*exp(-sqrt(-1)*(knL*(Nx/4+i)/Nx*4-j)));
    
%     Wa(:,Nx-i+1) = WW*(an.*exp(sqrt(-1)*knL*(-i)/Nx));
%     Wd(:,i) = -WW*(dn.*exp(-sqrt(-1)*knL*(100+i)/Nx));
%     
%     Wb(:,i) = -WW*(bn.*exp(-sqrt(-1)*knL*(i-1)/Nx));
%     Wc(:,i) = WW*(cn.*exp(sqrt(-1)*knL*(i-1)/Nx));
  end;
  
  %% central part
  Nx=100;  
  for i=1:Nx;
    Ub(:,i) = UU*(bn.*exp(-sqrt(-1)*(knL*(i-1)/Nx-j)));
    Uc(:,i) = UU*(cn.*exp(sqrt(-1)*(knL*(i-1)/Nx+j)));
    
%     Wa(:,Nx-i+1) = WW*(an.*exp(sqrt(-1)*knL*(-i)/Nx));
%     Wd(:,i) = -WW*(dn.*exp(-sqrt(-1)*knL*(100+i)/Nx));
%     
%     Wb(:,i) = -WW*(bn.*exp(-sqrt(-1)*knL*(i-1)/Nx));
%     Wc(:,i) = WW*(cn.*exp(sqrt(-1)*knL*(i-1)/Nx));
  end;
% %   
% %   
% %   
  %figure(dop);clf
  clf
  %CC = conv2(real([Ua Ub+Uc Ud]+1),ones(4,1)/4,'same'); %why plus 1?
  CC = conv2(real([Ua Ub+Uc Ud]),ones(4,1)/4,'same'); %% no BT
  %CC = conv2(abs([Ua Ub+Uc Ud]),ones(4,1)/4,'same'); %% amplitude
  %CC = angle(conv2([Ua Ub+Uc Ud]+1,ones(4,1)/4,'same'));
  %CC = conv2(real(Ud+1),ones(4,1)/4,'same');
  
  %CC = real([Ua Ub+Uc Ud]+1);
  %CC = angle([Ua Ub+Uc Ud]+1);
  %size(CC)
%  imagesc(CC);

%   imagesc([-400:500],[1 nn],CC);
%   hold on;
%   line([0 0],nn-[0 gam0]*nn,'linewi',2,'color','k')
%   line([0 0]+100,nn-[0 gamL]*nn,'linewi',2,'color','k')
%   caxis([-1 1]*5)
%   axis([-100 200 0 301])
%   %caxis([-pi pi])

ix=20;
iy=10;
xx = repmat(1:ix:size(CC,2),[length(1:iy:size(CC,1)) 1]); %size(xx)
yy = repmat([1:iy:size(CC,1)]',[1 length(1:ix:size(CC,2))]); %size(yy)
u = CC(1:iy:end,1:ix:end); %size(u)
v = zeros(size(CC(1:iy:end,1:ix:end))); %size(v)

imagesc(xx(1,:),yy(1,:),CC);
hold on
quiver_scale(xx,yy,u,v,4,'k')
axis([0 size(CC,2) 0 size(CC,1)])
caxis([-1 1]*5)
set(gca,'xlim',[300 600])


  hold off
  drawnow
  
  CCstd = CCstd + CC.^2;
  
  %% store data
  CC1 = [CC1 CC(97,375)];
  CC2 = [CC2 CC(97,410)];
  
  end
  
  figure
  subplot(2,1,1)
  plot(CC1)
  
  subplot(2,1,2)
  plot(CC2)
  
  
  figure
  plot([-399:500],CCstd(190,:))
  %plot([-799:1000],CCstd(200,:))
  
%   figure
%   imagesc([-400:500],[1 nn],sqrt(CCstd));
%   hold on;
%   line([0 0],nn-[0 gam0]*nn,'linewi',2,'color','k')
%   line([0 0]+100,nn-[0 gamL]*nn,'linewi',2,'color','k')
%   caxis([-1 1]*5)
%   axis([-100 200 0 301])
  
  

% % %% DEFAULT
% %   clear Ua Ud Uc Ub Wa Wb Wc Wd;
% %   Nx=300;
% %   for i=1:Nx;
% %     Ua(:,Nx-i+1) = UU*(an.*exp(sqrt(-1)*knL*(-i)/Nx));
% %     Ud(:,i) = UU*(dn.*exp(-sqrt(-1)*knL*(Nx+i)/Nx));
% %     
% %     Ub(:,i) = UU*(bn.*exp(-sqrt(-1)*knL*(i-1)/Nx));
% %     Uc(:,i) = UU*(cn.*exp(sqrt(-1)*knL*(i-1)/Nx));
% %     
% %     Wa(:,Nx-i+1) = WW*(an.*exp(sqrt(-1)*knL*(-i)/Nx));
% %     Wd(:,i) = -WW*(dn.*exp(-sqrt(-1)*knL*(100+i)/Nx));
% %     
% %     Wb(:,i) = -WW*(bn.*exp(-sqrt(-1)*knL*(i-1)/Nx));
% %     Wc(:,i) = WW*(cn.*exp(sqrt(-1)*knL*(i-1)/Nx));
% %   end;
% % 
% %   CC = conv2(real([Ua Ub+Uc Ud]+1),ones(4,1)/4,'same');
% %   size(CC)
% %   imagesc([-100:200],[1 nn],CC);
% %   hold on;
% %   line([0 0],nn-[0 gam0]*nn,'linewi',2)
% %   line([0 0]+100,nn-[0 gamL]*nn,'linewi',2)
% %   caxis([-1 1]*10)
% % %  colormap(colCog(32));


end;




return;
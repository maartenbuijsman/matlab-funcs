%%---------------------------------------------------------------------------------------------------
%%-- Program: tidefreq
%%-- Author: Herman Ridderinkhof
%%-- Modified by Maarten Buijsman 
%%-- Place: NIOZ
%%-- Date: 24-9-2002
%%-- Description: 
%%--determine tidal frequencies (based on 'analyse van getijden', Kalkwijk)
%%--coefficienten harmonische componenten (tabel 1, blz.22)
%%---------------------------------------------------------------------------------------------------
clear all
%basis getijden
csa = [0 0 1 0 0 -1];ct(1,:)=csa;
cssa =[0 0 2 0 0 0];ct(2,:)=cssa;
cmm = [0 -2 0 -1 0 0];ct(3,:)=cmm;
cmf = [0 2 0 0 0 0];ct(4,:)=cmf;
c2q1= [1 -3 0 2 0 0];ct(5,:)=c2q1;
csig1=[1 -3 2 0 0 0];ct(6,:)=csig1;
cq1=  [1 -2 0 1 0 0];ct(7,:)=cq1;
cro1= [1 -2 2 -1 0 0];ct(8,:)=cro1;
co1=  [1 -1 0 0 0 0];ct(9,:)=co1;
ctau1=[1 -1 2 0 0 0];ct(10,:)=ctau1;
cm1=  [1 0 0 -1 0 0];ct(11,:)=cm1;
cno1= [1 0 0 1 0 0];ct(12,:)=cno1;
cx1=  [1 0 2 -1 0 0];ct(13,:)=cx1;
cpi1= [1 1 -3 0 0 1];ct(14,:)=cpi1;
cp1=  [1 1 -2 0 0 0];ct(15,:)=cp1;
cs1=  [1 1 -1 0 0 1];ct(16,:)=cs1;
ck1=  [1 1 0 0 0 0];ct(17,:)=ck1;
cpsi1=[1 1 1 0 0 -1];ct(18,:)=cpsi1;
cphi1=[1 1 2 0 0 0];ct(19,:)=cphi1;
ctheta1=[1 2 -2 1 0 0];ct(20,:)=ctheta1;
cj1=  [1 2 0 -1 0 0];ct(21,:)=cj1;
coo1= [1 3 0 0 0 0];ct(22,:)=coo1;
ceps2=[2 -3 2 1 0 0];ct(23,:)=ceps2;
c2n2= [2 -2 0 2 0 0];ct(24,:)=c2n2;
cmu2= [2 -2 2 0 0 0];ct(25,:)=cmu2;
cn2=  [2 -1 0 1 0 0];ct(26,:)=cn2;
cnu2= [2 -1 2 -1 0 0];ct(27,:)=cnu2;
cm2=  [2 0 0 0 0 0];ct(28,:)=cm2;
clabda2=[2 1 -2 1 0 0];ct(29,:)=clabda2;
cl2=  [2 1 0 -1 0 0];ct(30,:)=cl2;
ct2=  [2 2 -3 0 0 1];ct(31,:)=ct2;
cs2=  [2 2 -2 0 0 0];ct(32,:)=cs2;
cr2=  [2 2 -1 0 0 -1];ct(33,:)=cr2;
ck2=  [2 2 0 0 0 0];ct(34,:)=ck2;
cksi2=[2 3 -2 1 0 0];ct(35,:)=cksi2;
ceta2=[2 3 0 -1 0 0];ct(36,:)=ceta2;
cm3=  [3 0 0 0 0 0];ct(37,:)=cm3;
%basis frequenties in graden/uur (blz 23)
%omegab=[14.492053 .549016 .041069 .004642 .002206 .000002]';
omegab=[14.49205211 0.54901653 0.04106864 0.00464184 0.00220641 0.00000196]';%%MCB02 http://www-ocean.tamu.edu/education/common/notes/chap18.htm#table1
%omegab=[14.4920535205	0.5490165205	0.0410686438	0.0046418219	0.0022063973	0.0000019726]';%%MCB02 Godin, 1972
nam1=char('Sa','Ssa','Mm','Mf','2Q1','sig1','q1','ro1','O1','tau1','M1','NO1');
nam2=char('x1','pi1','P1','S1','K1','psi1','phi1','the1','J1','OO1','eps2','2N2','mu2');
nam3=char('N2','nu2','M2','lab2','L2','T2','S2','R2','K2','ksi2','eta2','M3');
clear namesbt;
namesbt=nam1;
namesbt(13:25,1:4)=nam2;namesbt(26:37,1:4)=nam3;
%frequentie van de basis-componenten in radialen per dag
fac=(2*pi*24)/360.;
for i=1:length(ct)
   fqbt(i)=fac*ct(i,:)*omegab;
end
%nu ook de ondiep water getijden invoeren
%cno1=[1 0 0 1 0 0];cot(1,:)=cno1;%
cso1=[1 3 -2 0 0 0];cot(1,:)=cso1;
coq2=[2 -3 0 3 0 0];cot(2,:)=coq2;
%c2ms2=[2 -2 2 0 0 0];cot(4,:)=c2ms2;%
cop2=[2 0 -1 0 0 1];cot(3,:)=cop2;
cmks2=[2 0 2 0 0 0];cot(4,:)=cmks2;
%c2mn2=[2 1 0 -1 0 0];cot(7,:)=c2mn2;%
%cmsn2=[2 3 -2 -1 0 0];cot(8,:)=cmsn2;%
c2sm2=[2 4 -4 0 0 0];cot(5,:)=c2sm2;
cmo3=[3 -1 0 0 0 0];cot(6,:)=cmo3;
cso3=[3 1 -2 0 0 0];cot(7,:)=cso3;
cmk3=[3 1 0 0 0 0];cot(8,:)=cmk3;
csk3=[3 3 -2 0 0 0];cot(9,:)=csk3;
cmn4=[4 -1 0 1 0 0];cot(10,:)=cmn4;
cm4=[4 0 0 0 0 0];cot(11,:)=cm4;
csn4=[4 1 -2 1 0 0];cot(12,:)=csn4;
cms4=[4 2 -2 0 0 0];cot(13,:)=cms4;
cmk4=[4 2 0 0 0 0];cot(14,:)=cmk4;
cs4=[4 4 -4 0 0 0];cot(15,:)=cs4;
csk4=[4 4 -2 0 0 0];cot(16,:)=csk4;
c2mn6=[6 -1 0 1 0 0];cot(17,:)=c2mn6;
cm6=[6 0 0 0 0 0];cot(18,:)=cm6;
cmsn6=[6 1 -2 1 0 0];cot(19,:)=cmsn6;
c2ms6=[6 2 -2 0 0 0];cot(20,:)=c2ms6;
c2mk6=[6 2 0 0 0 0];cot(21,:)=c2mk6;
c2sm6=[6 4 -4 0 0 0];cot(22,:)=c2sm6;
csmk6=[6 4 -2 0 0 0];cot(23,:)=csmk6;
c3mn8=[8 -1 0 1 0 0];cot(24,:)=c3mn8;
cm8=[8 0 0 0 0 0];cot(25,:)=cm8;
c2msn8=[8 1 -2 1 0 0];cot(26,:)=c2msn8;
c3ms8=[8 2 -2 0 0 0];cot(27,:)=c3ms8;
c2ms8=[8 4 -4 0 0 0];cot(28,:)=c2ms8;
c2msk8=[8 -4 -2 0 0 0];cot(29,:)=c2msk8;
nam1=char('SO1  ','OQ2','OP2','MKS2','2SM2','MO3','SO3','MK3','SK3','MN4','M4','SN4');
nam2=char('MS4  ','MK4','S4','SK4','2MN6','M6','MSN6','2MS6','SMK6','2SM6','MSK6');
nam3=char('3MN8','M8','2MSN8','3MS8','2MS8','2MSK8');
clear namesot;
namesot=nam1;
namesot(13:23,1:5)=nam2;namesot(24:29,1:5)=nam3;

%frequentie van de basis-componenten in radialen per dag
fac=(2*pi*24)/360.;
for i=1:length(cot)
   fqot(i)=fac*cot(i,:)*omegab;
end
fq=[fqbt fqot];
save C:\data\TESO_data\dataout\tidalfrequencies\basisfreq namesbt fqbt namesot fqot fq

for i=1:length(fq)
   tdag(i)=(2*pi)/fq(i);tuur(i)=tdag(i)*24;
end

%lengte in dagen van benodigde tijdsreeks o.b.v. verschil tussen frequenties
for i=1:length(fqbt)-1
   dfbt(i)=fqbt(i+1)-fqbt(i);
   tbt(i)=2*pi/dfbt(i);
end
for i=1:length(fqot)-1
   dfot(i)=fqot(i+1)-fqot(i);
   tot(i)=2*pi/dfot(i);
end



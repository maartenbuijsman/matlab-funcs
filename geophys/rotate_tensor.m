function [Txxr,Txyr,Tyxr,Tyyr] = rotate_tensor(Txx,Txy,Tyx,Tyy,alp)
%% MCB, NRL, 2012-12-21
% T are tensor components and alp is the angle, CCW from positive x
% all input needs to be of the same size
% see comments 73
% Microsoft Word - Mathematics_of_CM_13_Coordinate_Transformation_Tensors.doc
% http://www.google.com/url?sa=t&rct=j&q=&esrc=s&source=web&cd=1&ved=0CDcQFjAA&url=http%3A%2F%2Fhomepages.engineering.auckland.ac.nz%2F~pkel015%2FSolidMechanicsBooks%2FPart_III%2FChapter_1_Vectors_Tensors%2FVectors_Tensors_13_Coordinate_Transformation_Tensors.pdf&ei=-nnUULPvBafV0gHAi4G4CA&usg=AFQjCNFZGANYkZMHHUvo__DbuD2Qr2gcrg&sig2=n1zBY79f3f2rzIK28aDBFw&bvm=bv.1355534169,d.dmQ

%% full form coord tran tensor
%  T11 T12   = Txx Txy   = T(1,1) T(1,2)
%  T21 T22     Tyx Tyy     T(2,1) T(2,2)

%% full form coord tran tensor
%  is the best
Txxr =  Txx.*cos(alp).^2 +Tyy.*sin(alp).^2 + 1/2*(Txy+Tyx).*sin(2*alp);
Txyr =  Txy.*cos(alp).^2 -Tyx.*sin(alp).^2 + 1/2*(Tyy-Txx).*sin(2*alp);
Tyxr = -Txy.*sin(alp).^2 +Tyx.*cos(alp).^2 + 1/2*(Tyy-Txx).*sin(2*alp);
Tyyr =  Txx.*sin(alp).^2 +Tyy.*cos(alp).^2 - 1/2*(Txy+Tyx).*sin(2*alp);

% %% alternatively
% %  Txy = Tyx
% Txxr =  Txx.*cos(alp).^2 +Tyy.*sin(alp).^2 + Txy.*sin(2*alp);
% Tyyr =  Txx.*sin(alp).^2 +Tyy.*cos(alp).^2 - Txy.*sin(2*alp);
% Txyr =  (Tyy-Txx).*sin(alp).*cos(alp)      + Txy.*cos(2*alp);
% Tyxr = Txyr;

return

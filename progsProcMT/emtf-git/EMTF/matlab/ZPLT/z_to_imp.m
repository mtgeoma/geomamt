%**********************************************************************
% Translates input from general Z file to one or more impedance matrices
% with signal and noise covariance matrices necessary for full error
% computation in any coordiante system.
%  ixy(2,# of impedance matrices) gives component numbers (first for Ex,
%    second for Ey) to be used for each impedance matrix to be extracted
%    NOTE: in a profilling setup with more dipoles along strike (x) than
%     across (y), there may be some y components used for more than one
%    impedance
%  orient(nch) gives orientations for all channels ... first two channels
%    are local reference (H) channels
%  ouput order is Zxx, Zxy, Zyx, Zyy
%
% Usage:   [Z2x2,SIG_S,SIG_E] = z_to_imp(Z,Sig_e,Sig_s,Nche,ixy,orient);
%  
%*********************************************************************

function [Z2x2,SIG_S,SIG_E] = z_to_imp(Z,Sig_e,Sig_s,Nche,ixy,orient)

ntemp = size(Z);
Nz = ntemp(2);
Nbt = ntemp(2)/Nche;
ntemp = size(ixy);
Nimp = ntemp(2);
Nt = Nbt*Nimp;
index_Z = [0:Nche:Nz-Nche];
orient = orient*pi/180;
%  make rotation matrix for H :
Th = [ cos(orient(1,1)) cos(orient(1,2));
       sin(orient(1,1)) sin(orient(1,2))];
Th = inv(Th);
Uh = kron(Th',Th');
Z2x2 = [];
SIG_E = [];
SIG_S = [];
for k=1:Nimp
%  make rotation matrix for E :
   Te = [ cos(orient(1,ixy(1,k)+2)) cos(orient(1,ixy(2,k)+2));
          sin(orient(1,ixy(1,k)+2)) sin(orient(1,ixy(2,k)+2))];
   Ue = kron(Te,Te);
   Uz = kron(Te,Th');
   Ztemp = [ Z(:,index_Z+ixy(1,k)); ...
                   Z(:,index_Z+ixy(2,k))] ;
   Z2x2 = [ Z2x2, Uz*Ztemp ];
   SIG_E = [ SIG_E , ...
     Ue*[ Sig_e(ixy(1,k),index_Z+ixy(1,k));
           Sig_e(ixy(1,k),index_Z+ixy(2,k));
           Sig_e(ixy(2,k),index_Z+ixy(1,k));
           Sig_e(ixy(2,k),index_Z+ixy(2,k))]];
   SIG_S = [SIG_S , Uh*[ Sig_s(:,1:2:2*Nbt-1);
                        Sig_s(:,2:2:2*Nbt) ] ];
  end
end

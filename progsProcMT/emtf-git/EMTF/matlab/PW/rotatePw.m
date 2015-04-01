
%  rotates pairs of data channels specified in array rot_ch(nrot,2),
%     modifying Pw.tf and residual covariance in Pw.cov
%  Channels listed in sing_ch are also included without rotation
%  Channels may be used for multiple rotations, so the total number
%   of channels output may exceed the input number of channels.
%   (but note that in this case the residual covariance matrix will
%   be singular)
%
%    Uses orientations in Pwhd.orient, and angle theta for new coordinate
%    axis.
%
%    Output is Pwrot, a data structure of the same form and number of 
%    bands as the input structure.
%
%  Usage :   [Pwrot] = rotatePw(Pw,Pwhd,rot_ch,sing_ch,theta);
%
function [PwRot] = rotatePw(Pw,Pwhd,rot_ch,sing_ch,theta)

%   set up full transformation matrix between output and input TF

nt = length(Pw.tf(1,:,1))
nrot = length(rot_ch(:,1));
T = Pw.T;
nbt = length(T);
nf = Pw.nf;
xxinv = Pw.xxinv;
nsing = length(sing_ch);
ornt = pi*(Pwhd.orient(1,:) - theta)/180;
co = cos(ornt);
so = sin(ornt);

row = []; col = []; s = []; s_inv = []; i1 = 1; i2 = 2;

for k=1:nrot
   row = [row i1 i1 i2 i2 ];
%  matrix to convert vector components from standard (rotated) orthogonal
%  coordinate system to arbitrary (not necessarily orthogonal)
%  measurement coordinate system
   ccmat = [co(rot_ch(k,1)) so(rot_ch(k,1)); ...
            co(rot_ch(k,2)) so(rot_ch(k,2))];
%  invert to get transformation from arbitrary measurement coordinates
%  into standard rotated orthogonal coordinate system
   ccmat = inv(ccmat);
   col = [col rot_ch(k,:) rot_ch(k,:) ];
   s = [s  reshape(ccmat',1,4) ];
   nout = i2;
   i1 = i1 + 2 ; i2 = i2 + 2;
end
for k = 1:nsing
   nout = nout + 1
   row = [ row nout ];
   col = [ col sing_ch(k) ];
   s = [s 1 ];
   s_inv = [s_inv 1 ];
end
U = sparse(row,col,s,nout,nt);
tf_rot = zeros(2,nout,nbt)+i*zeros(2,nout,nbt);
cov_rot = zeros(nout,nout,nbt)+i*zeros(nout,nout,nbt);
for k = 1:nbt
   tf_rot(:,:,k) = Pw.tf(:,:,k)*U';
   cov_rot(:,:,k) = U*Pw.cov(:,:,k)*U';
end

PwRot = struct('T',T,'nf',nf,'tf',tf_rot,...
   'xxinv',xxinv,'cov',cov_rot);

return

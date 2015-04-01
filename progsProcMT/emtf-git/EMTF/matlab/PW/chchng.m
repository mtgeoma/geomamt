%   script executed to change channel pairing choices

nrot = length(rot_ch(:,1));
nt = pwhd.nt;
lk = get(gco,'UserData');
l = lk(1);k = lk(2);
chnum1 = str2num(get(gco,'String'));
if(length(chnum1) == 0  & k <= nrot )
   % remove a channel pair
   ind1 = [1:nrot] ~= k;
   rc1 = rot_ch(k,:);
   rot_ch = rot_ch(ind1,:);
   if(sum( sum(rc1(1) == rot_ch))+sum(rc1(1) == sing_ch) == 0 )
      sing_ch = [ sing_ch rc1(1) ];
   end
   if(sum( sum(rc1(2) == rot_ch))+sum(rc1(2) == sing_ch) == 0 )
      sing_ch = [ sing_ch rc1(2) ];
   end
elseif(chnum1 <= nt & chnum1 > 0 & k <= nrot )
   % modify a channel pair
   rc1 = rot_ch(k,l);
   rot_ch(k,l) = chnum1;
   ind = (chnum1 ~= sing_ch);

   if(sum( sum(rc1 == rot_ch))+sum(rc1 == sing_ch) == 0 )
      sing_ch = [ sing_ch(ind) rc1 ];
   else
      sing_ch = sing_ch(ind);
   end

 elseif(chnum1 <= nt & chnum1 > 0 & k > nrot & l == 2 )
   % add channel pair
   rot_ch = [rot_ch ; sing_ch(k-nrot) chnum1 ];
   nrot = length(rot_ch(:,1));
   ind = sing_ch ~= rot_ch(nrot,1)  & sing_ch ~= rot_ch(nrot,2) ;
   sing_ch = sing_ch(ind);
end

ref_ch = rot_ch(1,:); 
nrot =  length(rot_ch(:,1));
nsing = length(sing_ch);

% paired channels
chnum(:,1:nrot) = rot_ch';
chname(1,1:nrot,:) = name_list(rot_ch(:,1),:);
chname(2,1:nrot,:) = name_list(rot_ch(:,2),:);
 
% single channels
i1 = nrot + 1;
i2 = nrot+nsing;
chnum(1,i1:i2) = sing_ch;
chnum(2,i1:i2) = 0;
for k=i1:i2
   chname(2,k,:) = name_list(nt+1);
end
chname(1,i1:i2,:) = name_list(sing_ch,:);
for k = i2+1:nlines
   chname(1,k,:) = name_list(nt+1);
   chname(2,k,:) = name_list(nt+1);
end
chnum(:,i2+1:nlines) = 0;
 
%  reset names and numbers in all checkboxes and edit fields
for k = 1:nlines
   for l = 1:2
      set(chhand(1,l,k),'string',setstr(chname(l,k,:)));
      set(chhand(2,l,k),'string',numstrbl(chnum(l,k)));
      set(chhand(1,l,k),'Value',0);
   end
end
chuse = zeros(2,nlines);

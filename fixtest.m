% vznum=length(dosmap2);
% nmax=0;
% for i=1:vznum
%     nn=length(dosmap2{i});
%     if nn>nmax
%         nmax=nn;
%     end
% end
% nmax=ceil(nmax/2)*2;
% dosmap=ones(nmax,vznum);
% for i=1:vznum
%     len=length(dosmap2{i});
%     if mod(len,2)==0
%         offset=nmax/2-len/2;
%         dosmap(1+offset:len+offset,i)=dosmap2{i};
%     else
%         offset=nmax/2-(len+1)/2;
%         dosmap(1+offset:(len+1)/2+offset,i)=dosmap2{i}(1:(len+1)/2);
%         dosmap(1+(len+1)/2+offset:1+len+offset,i)=dosmap2{i}((len+1)/2:len);
%     end
% end
% dosmap(dosmap==1)=nan;
%     

% for i=1:vznum
%     dosmap3{i}=dosmap2{i}.';
% end

fid = fopen('dosmap2.dat','w');
for i=1:vznum
    fprintf(fid,'%f ', dosmap3{i});
    fprintf(fid,'\n');
end
fclose(fid);
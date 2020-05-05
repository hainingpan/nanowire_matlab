function [re,val]=ispeak(x,y,tol)
%test if there is a peak near 0, within the tolerance of tol
[pks,loc]=findpeaks(y,x);
[pk,I]=mink(loc,1,'ComparisonMethod','abs');
if isempty(pk)
    re=false;
    val=nan;
else
    re=abs(pk)<tol;
    if re==true
        val=pks(I);
    else
        val=nan;
    end  
end
end

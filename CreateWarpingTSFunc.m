function [nL,LcYdelay,follSimVal] = CreateWarpingTSFunc(L,Y)
%CREATEWARPINGTSFUNC Summary of this function goes here
%   Detailed explanation goes here
% Y follows L
L=L(:)';
Y=Y(:)';
[~,wp]=DTW2(L,Y);
T=length(Y);
LcYdelay=nanmedian(wp);
follSimVal=nanmean(sign(wp));
if LcYdelay>0
    nL=ones(1,T).*median(L);
    wp=wp(T:end);
    for t=1:length(wp)
        currWarpBack=wp(t);
        if t-currWarpBack>0 && t-currWarpBack<=length(wp)
        inx=t-wp(t-currWarpBack);
            if inx >0 && inx<T
                nL(t)=L(inx);
            end
        end
    end
elseif LcYdelay<0
    tmp=L;
    L=Y;
    Y=tmp;
    nL=ones(1,T).*median(L);
    wp=wp(T:end);
    for t=1:length(wp)
        if t+wp(t)<=length(wp)
        inx=t+wp(t);
            if inx >0 && inx<T
                nL(t+inx)=L(t);
            end
        end
    end
else 
    LcYdelay=0;
    nL=L;
end

end


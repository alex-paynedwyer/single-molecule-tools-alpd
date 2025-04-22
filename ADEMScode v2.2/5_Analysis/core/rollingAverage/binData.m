function AvData=binData(X,Y,numLim,num,Colour)
if nargin<3
    numLim=0;
end
if nargin<5
    Colour=rand(1,3);
end
    if nargin<4
        [counts, edges,index]=histcounts(X);
    else
        [counts, edges,index]=histcounts(X,num);
    end
for n=1:max(index)
    if sum(index==n)>0
meanY(n)=mean(Y(index==n));
%meanX(n)=mean(X(index==n));
meanX(n)=min(edges)+n*mean(edges(2:end)-edges(1:end-1));
stdY(n)=std(Y(index==n));
numY(n)=length(Y(index==n));
errorY(n)=stdY(n)/numY(n)^0.5;
    else
        meanY(n)=0;
        meanX(n)=0;
stdY(n)=0;
numY(n)=0;
errorY(n)=0;
    end
end

errorbar(meanX(numY>numLim), meanY(numY>numLim),errorY(numY>numLim),'LineWidth',2,'color',Colour);
AvData=[meanX(numY>numLim)', meanY(numY>numLim)',errorY(numY>numLim)',numY(numY>numLim)'];

end
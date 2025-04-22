function cdf=CDFeval(DiffusionConstNuc)

cdf(1,:)=unique(DiffusionConstNuc);
for i=1:length(unique(DiffusionConstNuc))
    cdf(2,i)=sum((DiffusionConstNuc<cdf(1,i)))/length(DiffusionConstNuc);
end
end
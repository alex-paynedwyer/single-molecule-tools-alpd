function spotsOut=SpotCellCropper(spots,CellObject)
spotsOut=[];
if sum(CellObject,'All')>0
for c=1:size(CellObject,3)
    cellProps=regionprops(CellObject(:,:,c),'boundingbox');
    spotC=spots;
    spotC(spotC(:,1)<cellProps.BoundingBox(1),:)=[];
    spotC(spotC(:,1)>(cellProps.BoundingBox(1)+cellProps.BoundingBox(3)),:)=[];
    spotC(spotC(:,2)<cellProps.BoundingBox(2),:)=[];
    spotC(spotC(:,2)>(cellProps.BoundingBox(2)+cellProps.BoundingBox(4)),:)=[];
    spotsOut=cat(1,spotsOut,spotC);
end
end

end
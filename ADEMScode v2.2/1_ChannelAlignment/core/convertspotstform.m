function spotstform=convertspotstform(tformin)
spotstform=tformin;
xsteps=tformin.T(3,1); %swap x and y steps to compensate for MATLAB y,x convention
ysteps=tformin.T(3,2);
spotstform.T(2,1)=tformin.T(2,1)*-1;
spotstform.T(1,2)=tformin.T(1,2)*-1;
spotstform.T(3,1)=-ysteps;
spotstform.T(3,2)=-xsteps;
end
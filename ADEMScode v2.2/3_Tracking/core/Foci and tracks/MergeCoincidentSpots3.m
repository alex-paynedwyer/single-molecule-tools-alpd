function [spots2, new_spots]=MergeCoincidentSpots3(spots, framenum1, d_min)
 %Function calculated pairwise distances and merges any closer
  %than dmin. Returns a spot array with merged spot centre and any
  %redundant spot numbers set to very high co-ordinates. Also the spot
  %numbers for the new spots so their centres can be determined in main
  %code.

if isempty(spots(spots(:,9)==framenum1))==0  % at least 1 accepted spot in frame
    spots_int=spots;
    spots_to_refind=[];
    
    N1 = size(spots(spots(:,9)==framenum1),1);  % no. of accepted spots in current frame
 
    

    %counters=bad programming
    m1=1;
    % loop though accepted spots in previous frame.
   m2=0;
   % Keep going until no spots are linked
   while m1>m2
       m2=m1;
       %loop over all spots in image
       for q0 = min(find(spots_int(:,9)==framenum1)):max(find(spots_int(:,9)==framenum1))
           if  spots_int(q0,1)~=100000
               %loop over all spots in image
               for q1 = min(find(spots_int(:,9)==framenum1)):max(find(spots_int(:,9)==framenum1))
                   if spots_int(q1,1)~=100000
                       % don't try to link spots to themselves
                       if q0==q1
                       else
                           % d_01: distance between spot Centres
                           d_01 = sqrt((spots_int(q0,1)-spots_int(q1,1))^2+(spots_int(q0,2)-spots_int(q1,2))^2);
                           
                           
                           %If spots are too close, calculate average and assign
                           %to q0, assign q1 to giant number to ignore it
                           if d_01 < d_min
                               spots_int(q0,1)=(spots_int(q0,1)+spots_int(q1,1))/2;
                               spots_int(q0,2)=(spots_int(q0,2)+spots_int(q1,2))/2;
                               %Set the now redundant spots centre to 100000
                               spots_int(q1,1)=100000;
                               %store the spot number
                               spots_to_refind(m1,:)=[q0,q1];
                               %   spots_to_remove(m1)=q1;
                               m1=m1+1;
                           end
                       end
                   end
               end
           end
       end
       % Remove all the coincident spots
       %   spots_int(spots_int(:,1)==100000,:)=[];
   end
    %assign spot matrix without coincident spots
    spots2=spots_int;
    if spots_to_refind>0
    %output only spot numbers which will need re-fitting
    new_spots=unique(spots_to_refind(:,1));
    else
        new_spots=0;
    end
else
    %if no spots output inputs
    spots2=spots;
new_spots=0;
end

end
function LBW = Watershed(I)
    %% 
    Distance=bwdist(~I);
    Distance=-Distance;
    Distance = imhmin(Distance,1);
    L=watershed(Distance);
    L(~I) = 0;
    LBW=logical(L);
end

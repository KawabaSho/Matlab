function det = det33(RR)
    det =  (RR(1,1)*RR(2,2)*RR(3, 3) + RR(1, 2)*RR(2, 3)*RR(3, 1) + RR(1, 3)*RR(2, 1)*RR(3, 2))...
           -(RR(1, 3)*RR(2, 2)*RR(3, 1) + RR(1, 2)*RR(2, 1)*RR(3, 3) + RR(1, 1)*RR(2, 3)*RR(3, 2));
end
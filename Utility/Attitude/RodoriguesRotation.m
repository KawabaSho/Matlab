function Rotation = RodoriguesRotation(x, rad)
    cs = cos(rad); sn = sin(rad);
    Rotation = ...
    [cs+x(1)*x(1)*(1-cs), x(1)*x(2)*(1-cs)-x(3)*sn, x(1)*x(3)*(1-cs)+x(2)*sn;
     x(1)*x(2)*(1-cs)+x(3)*sn, cs+x(2)*x(2)*(1-cs), x(2)*x(3)*(1-cs)-x(1)*sn;
     x(1)*x(3)*(1-cs)-x(2)*sn, x(2)*x(3)*(1-cs)+x(1)*sn, cs+x(3)*x(3)*(1-cs)];
    Rotation = Rotation/det33(Rotation);
end
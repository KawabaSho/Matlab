function Vertices = generateCuboid_Vertices(L,W,H)
    % patch
    Vertices = ...
        [ L, W, H;
         -L, W, H;
         -L,-W, H;
          L,-W, H;
          L,-W,-H;
         -L,-W,-H;
         -L, W,-H;
          L, W,-H; ]*0.5;
end
function model = generateCube(L)

    model.Vertices = ...
        [ 1, 1, 1;
         -1, 1, 1;
         -1,-1, 1;
          1,-1, 1;
          1,-1,-1;
         -1,-1,-1;
         -1, 1,-1;
          1, 1,-1; ].*L;
    model.Faces = ...
        [   1	2	3
            1	3	4
            5	4	3
            5	3	6
            6	3	2
            6	2	7
            7	8	5
            7	5	6
            8	1	4
            8	4	5
            7	2	1
            7	1	8 ];
    % patch(model,'FaceColor', [0.8 0.8 1.0])
    % view([20 20])
end
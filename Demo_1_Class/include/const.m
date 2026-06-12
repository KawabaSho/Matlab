classdef const
   properties (Constant)
   %% Geo params
   g     = 9.81           % m/s^2
   % C = 2.99792458e8     % m/s
   % G = 6.67259          % m/kgs
   % Me = 5.976e24        % Earth mass (kg)
   % Re = 6.378e6         % Earth radius (m)

   %% Drone params
   m     = 1.608            % mass kg
   II    = [0.0366, 0, 0.0065;
            0, 0.0327, 0.0046;
            0.0065, 0.0046, 0.0187];
   invII = inv(const.II)
                            % Inertia kgm^2
   L     = 0.216            % Length of arm
   u_min = 0.1046           % thrust N
   u_max = 6                
   k     = 1.6e-2           % Coefficient of rotors m
   T     = [-const.L/sqrt(2), -const.L/sqrt(2),  const.L/sqrt(2),  const.L/sqrt(2);
            -const.L/sqrt(2),  const.L/sqrt(2),  const.L/sqrt(2), -const.L/sqrt(2);
                     const.k,         -const.k,          const.k,         -const.k]
   T_vec = [1,-1,1,-1];

   M     = [1,1,1,1; const.T]  % output matrix
   invM  = inv([1,1,1,1; const.T])  % mixing
   
   n_arm = [ 1/sqrt(2), -1/sqrt(2), -1/sqrt(2),  1/sqrt(2);
            -1/sqrt(2), -1/sqrt(2),  1/sqrt(2),  1/sqrt(2);
            0,0,0,0]

   % Render Drone
   axis_L = 0.5
   arm = const.n_arm.*const.L % vector
   arm_num = 4
   propeller_radius = 0.05
   body_width  = 0.07
   body_depth  = 0.20
   body_height = 0.03
   coord = [...
            0                         0    0;
            const.body_depth  0    0;
            const.body_depth  const.body_width  0;
            0                         const.body_width  0;
            0           0                         const.body_height;
            const.body_depth  0                         const.body_height;
            const.body_depth  const.body_width  const.body_height;
            0           const.body_width  const.body_height]...
      - [const.body_depth const.body_width const.body_height].*0.5;
   coord_x = const.coord(:,1)
   coord_y = const.coord(:,2)
   coord_z = const.coord(:,3)
   idx = [4 8 5 1 4; 1 5 6 2 1; 2 6 7 3 2; 3 7 8 4 3; 5 8 7 6 5; 1 4 3 2 1]'
   idx_num = numel(const.idx)
   ratio = 1 - const.propeller_radius/const.L
   teta = -pi:0.02:pi
   motor_plot = [const.propeller_radius*cos(const.teta);...
                 const.propeller_radius*sin(const.teta);...
                 zeros(1,numel(const.teta))];
   body_color = [0, 0.2, 0.4, 0.6, 0.8, 1]

   % camera params
   camera_cmos_ratio = 24.0 / 36.0      % full size
   FOV               = 40/180*pi        % [rad]
   camera_pos  = [const.body_depth*0.5; 0; 0]
   camera_vec  = [1;0;0]
   unit_vec    = RodoriguesRotation(cross([0;-1;-const.camera_cmos_ratio],const.camera_vec), const.FOV*0.5)*const.camera_vec
   camera_coord = [0,0,0;
                  const.unit_vec';
                  const.unit_vec(1),  const.unit_vec(2), -const.unit_vec(3);
                  const.unit_vec(1), -const.unit_vec(2), -const.unit_vec(3);
                  const.unit_vec(1), -const.unit_vec(2),  const.unit_vec(3);].*0.25'
   camera_idx   = [1,2,3,1;1,3,4,1;1,4,5,1;1,2,5,1]'


   %% Aero params
   gamma = 1e-2             % the coefficient of air resistance N m s/rad
   D     = [0,0,0;0,0,0;0,0,const.gamma]

   end
end













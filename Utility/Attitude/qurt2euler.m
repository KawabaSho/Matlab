function Euler = qurt2euler(q)
    % Euler = [phi(b_x); theta(b_y); psi(b_z)]
    q = q/(q(1)*q(1) + q(2)*q(2) + q(3)*q(3) + q(4)*q(4));
    Euler = [asin(2*(q(1)*q(2) + q(3)*q(4)));
             atan(2*(q(1)*q(3) - q(2)*q(4))/(q(1)^2 - q(2)^2 - q(3)^2 + q(4)^2));
             atan(2*(q(1)*q(4) - q(2)*q(3))/(q(1)^2 - q(2)^2 + q(3)^2 - q(4)^2))];
    % Euler = real(Euler);
end
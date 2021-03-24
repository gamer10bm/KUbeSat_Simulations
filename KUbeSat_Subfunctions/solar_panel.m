function SP = solar_panel(SC, R, Q, RS)
%  This function estimates solar panel collection based on the relative
%  geometry and vehicle model.
%  No efficiency of power conversion is accounted for in this model, and
%  basic estimate is derived from maximum rated power and cosine losses
%  from orientation.
%  Also no self-shadowing for s/c is assumed.

    % Determine transformation from body axes to inertial frame
    C = [1-2*(Q(2)^2+Q(3)^2), 2*(Q(1)*Q(2)+Q(3)*Q(4)), 2*(Q(1)*Q(3)-Q(2)*Q(4));
    2*(Q(1)*Q(2)-Q(3)*Q(4)), 1-2*(Q(1)^2+Q(3)^2), 2*(Q(2)*Q(3)+Q(1)*Q(4));
    2*(Q(1)*Q(3)+Q(2)*Q(4)), 2*(Q(2)*Q(3)-Q(1)*Q(4)), 1-2*(Q(1)^2+Q(2)^2)];
    C = C';

    % extract number of s/c sides
    n = length(SC.side);
    
    % initialize sum of solar power generation
    SP = 0;
    % test for eclipse by finding minimum of sun LOS to earth center
    tt = -dot(RS,R); 
    if (tt > 0) % directions to Earth (-R) and Sun (rS) are the same plane
        % of sky. Risk of Eclipse. Need to check if a ray starting at s/c
        % passes through Earth (eclipse) on the way toward Sun.
        % calculate minimum distance of such a ray.
        d = sqrt((R(1)+RS(1)*tt)^2 + (R(2)+RS(2)*tt)^2 + ...
            (R(3)+RS(3)*tt)^2);  
    else % directions to Sun and Earth are opposite direction
        % no risk of eclipse
        d = norm(R);
    end
    if (d > SC.Re) % in case of no eclipse
        for i = 1:n
            s = C*SC.side{i}; % convert side normal to inertial coords
            del = dot(s,RS); % cosine of side normal and sun vector
            if (del > 0) % side is in sunlight
                SP = SP + SC.SP(i)*del; % add power to total with cosine loss
            end
        end % panel loop k
    end % no eclipse check
end
X = Rocker_range(1)


function [X] = Rocker_Range(VTh7, rocker_link_angle, r_rocker_shockside)
    D2R = pi/180;
    num_changes = 5;
    
    figure()
    clf
    untouched_rocker_link_angle = rocker_link_angle;
    untouched_r_rocker_shockside = r_rocker_shockside;

    rock_angle = linspace(abs(rocker_link_angle) - 10, abs(rocker_link_angle) + 7, num_changes);     % The angle between the pullrod hole, the rocker axis, and the rocker shock mounting hole
    r_shock = linspace(r_rocker_shockside - 10,r_rocker_shockside + 10, num_changes);        % The radius from the rocker axis to the shock hole
    
    
%   This section iterates over an array change_vals which contains values for each optimization variable. Each 'slice' 
%   in the third dimension represents a different set of changes. 
%   For example, using the array below will plot the motion ratio vs. a varying rocker angle in the 1st iter, then a
%   varying r_shock in the 2nd iter.
%       change_vals = [mean(r_shock)*ones(size(rock_angle)); rock_angle];     % constant shock radius, varying rocker angle
%       change_vals(:,:,2) = [mean(rock_angle)*ones(size(r_shock)); r_shock];

%       optimization variables: [r_rocker_shockside, rocker angle]
    change_vals =        [untouched_r_rocker_shockside*ones(size(rock_angle)); rock_angle];     % constant shock radius, varying rocker angle
    change_vals(:,:,2) = [r_shock;                              untouched_rocker_link_angle*ones(size(r_shock))];
    
    clf % clear current figure
    changevars = 2; % number of variables to plot over
    for iter = 1:changevars
        subplot(2,1,iter)
        hold on
        for i = 1:length(change_vals)
            % Iterate variables
            r_rocker_shockside = change_vals(1,i,iter);
            rocker_link_angle = change_vals(2,i,iter);    % the angle between the two rocker links
            
            % Recalculate results
            shock_mount_pos = O_R -.75.*r_rocker_shockside + 1i*(25.4*7.7);   %Shock mount pos wrt to rocker rotation axis; shock max len = 260mm
            R9 = exp(1i*(VTh7 - D2R.*rocker_link_angle)) .* r_rocker_shockside;  % Rocker Shock-Side link
            G = O_R + R9;
            shock_length = ( real(shock_mount_pos - G).^2 + imag(shock_mount_pos - G).^2 ).^ 0.5;
            spring_dx = diff(shock_length);
            bump_velocity = diff(bump);
            motion_ratio = abs(spring_dx) ./ abs(bump_velocity);
            
            % plot
            plot(bump(2:end),motion_ratio)
            grid on
        end
        xticks(0-60:10:60)%, yticks(linspace( (round(min(spring_velocity),1,'significant')), max(spring_velocity) ,10 ) )
        
        if iter == 1
            legend(string(round(change_vals(2,:,iter),3,'significant') ) + " deg")
            title("MR v. Wheel Disp. - Varying rocker angle")
        elseif iter == 2
%             change_vals = [r_shock; 110*ones(size(r_shock))];
            legend(string(round(change_vals(1,:,iter),3,'significant') ) + " mm")
            title("MR v. Wheel Disp. - Varying shock rocker radius")
%         elseif iter == 3
        end
    end
end
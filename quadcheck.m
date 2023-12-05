function gamma = quadcheck(vec)
% vec must be at least a 1x4 array

    % if 1.e-10 > vec(1,1)
    %     gamma = 0;
    % elseif ((pi-1.e-10) < vec(1,1)) && ((pi+1.e-10) > vec(1,1))
    %     gamma = pi;
    % else
    % Primary loop
    gamma = [];
            % Finds the absolute differnce between each column in the array
    diff1 = abs(vec(1)-vec(2));
    diff2 = abs(vec(1)-vec(3));
    diff3 = abs(vec(1)-vec(4));
    diff4 = abs(vec(2)-vec(3));
    diff5 = abs(vec(2)-vec(4));
    diff6 = abs(vec(3)-vec(4));
    
    % If the absolute differnce meets the criteria, one of the columns
    % associated with it is saved to gamma
    if diff1 < 1.e-8
        gamma = vec(1);
    elseif diff2 < 1.e-8
        gamma = vec(1);
    elseif diff3 < 1.e-8
        gamma = vec(1);
    elseif diff4 < 1.e-8
        gamma = vec(2);
    elseif diff5 < 1.e-8
        gamma = vec(2);
    elseif diff6 < 1.e-8
        gamma = vec(3);
    else
        low = min([diff1 diff2 diff3 diff4 diff5 diff6]);
        if low == diff1
            gamma = vec(1);
        elseif low == diff2
            gamma = vec(1);
        elseif low == diff3
            gamma = vec(1);
        elseif low == diff4
            gamma = vec(2);
        elseif low == diff5
            gamma = vec(2);
        elseif low == diff6
            gamma = vec(3);
        end
    end     % end
end
function [ERA, P_v_check, count] = ERA_adjustment_to_P_v(ERA_total, UI_k, SB, P_v, n_PB, n_runs)

% Check probability of violation P_v
P_v_check = zeros(n_PB,1);
for i=1:n_PB
    n=0;
    for j=1:n_runs
        if ERA_total*UI_k(1,i,j)>SB(i,1,j)
            n=n+1;
        end
    end
    P_v_check(i) = n/n_runs;
end
count = 0;
if max(P_v_check)< P_v
    % increase ERA until P_v is reached, SoP stays constant
    while max(P_v_check)< P_v
        if ERA_total == 0
            break
        end
        ERA_total = ERA_total * (1+P_v-max(P_v_check));
        for i=1:n_PB
            n=0;
            for j=1:n_runs
                if ERA_total*UI_k(1,i,j)>SB(i,1,j)
                    n=n+1;
                end
            end
            P_v_check(i) = n/n_runs; 
        end
        count = count + 1;
    end
elseif max(P_v_check)> P_v
    % decrease ERA until P_v is reached, SoP stays constant
    while max(P_v_check)> P_v
        if ERA_total == Inf
            ERA_total = 0;
            break
        end
        ERA_total = ERA_total * (1-P_v * 5 * (1 + max(P_v_check)));
        for i=1:n_PB
            n=0;
            for j=1:n_runs
                if ERA_total*UI_k(1,i,j)>SB(i,1,j)
                    n=n+1;
                end
            end
            P_v_check(i) = n/n_runs; 
        end
        count = count - 1;
    end
    while max(P_v_check)< P_v
        if ERA_total == 0
            break
        end
        ERA_total = ERA_total * (1+P_v-max(P_v_check));
        for i=1:n_PB
            n=0;
            for j=1:n_runs
                if ERA_total*UI_k(1,i,j)>SB(i,1,j)
                    n=n+1;
                end
            end
            P_v_check(i) = n/n_runs; 
        end
        count = count - 1;
    end
else
    count = 'not scaled';
end

ERA = ERA_total;
end
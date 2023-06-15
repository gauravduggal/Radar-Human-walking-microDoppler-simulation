close all
clear all
clc
%pulses per cycle of human (1 left step and 1 right step)
nt = 2048;
radar_loc = [10,10,1];
scale = 3;
numcyc = 50;

moves = scale*[0.707,0.707;
            0,1;
            1,0;
            0,-1;
            -1,0];
boundary = [0,10;
            0,10];

start_pos = [0,0];
cur_pos = [0,0];
key_points = zeros(numcyc,2);
for nc = 1:numcyc

    key_points(nc,:) = cur_pos;
    new_pos = move(boundary,cur_pos, moves);
    
    cur_pos = new_pos;

%     plot(cur_pos(1),cur_pos(2),'*')
%     hold on

end

t = 0:1/(nt*scale):numcyc-1/(nt*scale);
t = t(1:numcyc*nt-1);
trajectory_x = spline( 0:numcyc-1,key_points(:,1),t);
% 
trajectory_y = spline( 0:numcyc-1,key_points(:,2),t);
% trajectory_x = trajectory_x(1:numcyc*nt-1);
% trajectory_y = trajectory_y(1:numcyc*nt-1);
plot(trajectory_x,trajectory_y)
xlim(boundary(1,:))
ylim(boundary(2,:))
trajectory = [trajectory_x; trajectory_y];
d = 0;
for nc=1:range(length(t))
    d = d + sqrt((trajectory_x(nc)-trajectory_x(nc+1))^2+(trajectory_y(nc)-trajectory_y(nc+1))^2);
end
d

save('trajectory.mat','trajectory','nt','numcyc','radar_loc')
function [new_pos] = move(boundary,cur_pos, moves)
        while(true)
        m = moves(randperm(size(moves,1),1),:);
        temp = cur_pos + m
        if (temp(1)>boundary(1,1) && temp(1)<boundary(1,2) && temp(2)>boundary(2,1) && temp(2)<boundary(2,2))
            new_pos = temp;
            break
        end
        end


end
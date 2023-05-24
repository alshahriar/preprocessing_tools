function [newX,newY] = findSingleLineXY(x,y,smallestDis)
nx = length(x);
count = 1;
for i = 1:nx-1
    distance(count) = sqrt((x(i) - x(i+1))^2 + (y(i) - y(i+1))^2);
    count = count + 1;
end
distance(count) = sqrt((x(nx) - x(1))^2 + (y(nx) - y(1))^2); % last point is connected to the 1st point
maxDis = max(distance);
if (maxDis>=2*smallestDis) %distance between the points should be large enough
    count=1;
    for i=1:nx-1
        % angle for inclined line
        angle = atand((1)*(y(i) - y(i+1))/(x(i) - x(i+1)));
        nNewPoint(i) = floor(distance(i)/smallestDis);% how many grid for line i
        nNewPoint = double(nNewPoint);
        tempDis = double(nNewPoint(i))*smallestDis; % length = gridSize*nOfGrid
        differ = distance(i)-tempDis;
        if (differ<(smallestDis/2) && differ>0) % if makeup portion is less than, this portion will be added to the last grid
            nNewPoint(i)= nNewPoint(i)-1;
        elseif (differ==0)
            nNewPoint(i)= nNewPoint(i)-1;
        end
        newX(count)=x(i); % starting point % one end
        newY(count)=y(i);
        countOld=count;
        if (x(i)>x(i+1))%current end is > the next end
            for j=count:count+nNewPoint(i)-1
                count=count+1;
                %so we will minus each distance
                newX(j+1)=newX(j) - smallestDis*cosd(angle);
                % we have to move this smallestDis*cosd(angle) amount
                %in x-direction
            end
        else
            for j=count:count+nNewPoint(i)-1
                newX(j+1)=newX(j) + smallestDis*cosd(angle);
                count=count+1;
            end
        end
        if (x(i)<x(i+1))
            for j=countOld:countOld+nNewPoint(i)-1
                newY(j+1)=newY(j) + smallestDis*sind(angle);
            end
        else
            for j=countOld:countOld+nNewPoint(i)-1
                newY(j+1)=newY(j) - smallestDis*sind(angle);
            end
        end
        count=count+1;
    end
end
end
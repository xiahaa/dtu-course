function mouseClick2(object, eventdata)
    global count m keyframez;
    if(count < m)
        count = count + 1;
        C = get(gca, 'CurrentPoint');
        keyframez(1,count) = C(1,1);
        keyframez(2,count) = C(1,2);               %z
    end
end
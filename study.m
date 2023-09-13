clear all
close all
clc
k = 0;
% shape = ["Sine", "I", "QSW", "QSW"; "HSW1", "S", "Sine", "U"];
shape = ["QSW", "I", "U", "QSW"; "HSW1", "S", "QSW", "L"];
for j = 1:size(shape,2)
    initial = shape(1,j);
    desired = shape(2,j);
    for i = 6:2:20
        k = k +1
        sochi_data(k) = dlo_update(i, initial, desired);
    end
end


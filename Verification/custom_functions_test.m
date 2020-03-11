clear; clc;
%   Test all converter fucntions


% Test seconds converter function
time = [1 05 32];
seconds = tosecs(time);
assert(seconds == 3932)

time = [0 00 -02];
seconds = tosecs(time);
assert(seconds == -2)

time = [0 00 00];
seconds = tosecs(time);
assert(seconds == 0)


% Test kg converter function
lbs = 1;
kg = tokg(lbs);
assert(kg == 0.45359237)

lbs = 0;
assert(tokg(0) == 0);


% Test knots to m/s converter function
assert(round(toms(1),4) == 0.5144);
assert(round(toms(-1),4) == -0.5144);


% Test feet to m converter function

assert(tom(1) == 0.3048);
assert(tom(-1) == -0.3048);


function [secs] = tosecs(time)
%TOSECS Calculate the seconds from a list of hours,mins,seconds
%   take each index of index and multiply by correct unit conversion, sum
%   to get total seconds
hh = time(1);
mm = time(2);
ss = time(3);
secs = hh * 60 * 60 + mm * 60 + ss;
end
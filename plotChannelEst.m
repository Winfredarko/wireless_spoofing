function plotChannelEst(chanest, thresh, varargin)
% first varargin should be how many channel estimates to plot
%
% channel estimate vector should be a vector of 512 elements, where the
% first 256 are real while latter 256 are the imaginary components

if nargin == 2
    max = size(chanest, 1);   % amount of rows in channelest
    mag = zeros(max, size(chanest, 2) / 2);
    
    for i = 1:max
        cs = chanest(i, :);
        len = size(chanest, 2) / 2;
        a = cs(1:len);
        b = cs(len+1:2*len);
        complexCS = complex(a, b);
        
        mag(i, :) = normalize(abs(complexCS));     % use angle(...) to get the phase angle
       
        % thresh is threshold above which channel vectors are excluded
        if isempty(find(mag(i, :) > thresh, 1))
            hold on
            plot(mag(i, :))
        end
    end
else
    cs = chanest(varargin{1}, :);
    len = size(chanest, 2) / 2;
    a = cs(1:len);
    b = cs(len+1:2*len);
    complexCS = complex(a, b);
    
    mag = normalize(abs(complexCS));
%     mag = abs(complexCS);
    plot(mag(128:192));
end
    
end
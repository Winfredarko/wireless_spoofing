%% Load in data
load 'channelstates.mat';

len = size(channelstates, 2) / 2;
n = size(channelstates, 1);

% processing channel states so they can be represented in complex form
complexCS = zeros(size(channelstates, 1), len);
magCS = zeros(size(channelstates, 1), len / 4);
for i = 1:n
    cs = channelstates(i, :);
    complexCS(i, :) = complex(cs(1:len), cs(len+1:2*len));
    magCS(i, :) = normalize(abs(complexCS(i, 129:192)));
end

% only focus on the correct channel
random = rand(4, len / 4) * 7 - 1;
total = [magCS; random];                  % insert spoofed channel states
shuffledTotal = total(randperm(size(total, 1)), :);

for i = 1:size(shuffledTotal, 1)
    hold on;
    plot(shuffledTotal(i, :));
end

%% Initializing relevant variables
G1 = 6;
C1 = 6;
G0 = 9;
C0 = 4;
y = 0.25;
mu = 0.8;
delta = 0.7;
eps = 0.1;
X = 29;

%% Frobenius norm
% go through all the channel vectors
rhat = magCS(1, :);     % assuming first vector is legitimate

for i = 1:size(shuffledTotal, 1)
    r = shuffledTotal(i, :);
    
    % calculates the Frobenius norm, denoted with L
    if isequal(rhat, r)
        continue
    end
    L = (norm(r - rhat))^2 / (norm(rhat))^2
end
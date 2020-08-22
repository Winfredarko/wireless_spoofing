% 1: Initialize: ?, ?, ?, Q(s,x)=0, V (s)=0,?x? {l/K}0?l?K
parseInput();

% receiverProb = ????, Probability that receiver chooses the suboptimal actions
receiverProb = 0.1;

% learnRate = ?, Learning rate of Q-learning
learnRate = 0.1;

% discountFactor = ?, Discount factor of Q-learning
discountFactor = 0.7;

% Q-function of the receiver choosing threshold x(quantized into 10 values)
% in s,(where P_f and P_m are also quantized into 10 values)
s_quantNum = 500;
x_quantNum = 500;

Q = zeros(s_quantNum, s_quantNum, x_quantNum);

% Best threshold to maximize Q value of the receiver in s
V = zeros(s_quantNum, s_quantNum);

% Expected Sum Utility
U_n = 0;

% Channel record for transmitter i
[vecRows, vecCols] = size(channelest);

rTrueSize = 10;
rTrue = zeros(rTrueSize, vecCols / 2);

channels = zeros(vecRows, vecCols/2);
for n = 1:vecRows
    currRow = complex(channelest(n, 1: vecCols/2), channelest(n, (vecCols/2 + 1): (vecCols)));
    channels(n, :) = currRow(:);
end

% Variales for Utility u_r

% Go, Payoff from accepting legitimate packet
G0 = 6;

% G1, Payoff from rejecting a spoofing one
G1 = 6;

% Co, Cost of rejecting legitimate packet
C0 = 4;

% C1, Cost of accepting spoofing packet
C1 = 2;

% Number of nodes M
M = 64;

% y, Spoofing frequency (There's only one spoofer, so arbritrary at this point)
y = 0.1;

% Number of Spoofers
numSpoofers = 1;

% Number of time slots
timeSlotNum = floor(vecRows / 10);

s_next = [0, 0];
s_n_range = 5;

vectorNum = 1;

rSpooferVec = zeros(rTrueSize, vecCols / 2);

spooferNum = zeros(timeSlotNum, 1);
thresholds = zeros(timeSlotNum, 1);
utilities = zeros(timeSlotNum, 1);

spoofers = 0;

Pf_max = 0;

 
% 2: % for n = 1,2,3,...do
for n = 1:timeSlotNum
    % 3: Observe the current state sn
    
    if n == 1
        X_n = 0.3;
        [~, X_index] = quantize(X_n, x_quantNum, s_n_range);
        [s_n, s_n_Index] = observeState(X_n, channels(vectorNum, :), rSpooferVec);%, numSpoofers);
        
    else
        s_n = s_next;
        [~, s_n_Index(1)] = quantize(s_n(1), s_quantNum - 1, s_n_range);
        [~, s_n_Index(2)] = quantize(s_n(2), s_quantNum - 1, s_n_range);
        % disp(s_n);
        
        disp("S_n index is: " + s_n_Index);
    end
    
    disp("-----");
    % 4: Select a threshold xn via(20)
    if n ~= 1
        X_n = X_next;
        disp(X_n);
        [X_n, X_index] = quantize(X_n, x_quantNum, s_n_range);
    end

    
    
    % 5: for k = 1 to T do
    
    % Number of packets received in time slot n
    T = 10; % getNumPackets()
    
    for k = 1:T
        % Channel vector of the k-th packet claming to be sent by the
        % transmitter i
        rClaim = channels(vectorNum, : );
        vectorNum = vectorNum + 1;
        
        % Assumes first time slot has legit packets
        if (n == 1)
            % disp(T);
            rGiven = rClaim;
            % if (k == 1) mac = readMac(k); end %#ok<SEPEX>
            if (k <= rTrueSize) rTrue(k, :) = rGiven; end %#ok<SEPEX>
        end
        
        % 6: Read the MAC address ?i of the kth packet
        % macGiven = readMac(k);
        
        % 7: Extract ri and ˆ ri
        
        % Gives an incorrect vector to test spoofer packet
%         if ((k == 10) && (rem(n, 10) == 0))
%             spoofers = spoofers + 1;
%             rGiven = randn(1, vecCols/2);
%         end
        
        % 8: Calculate L via(2)
        % L is the statistic of the hypothesis test in the spoofing detector
        [~, rClosestIndex] = min(rTrue(:) - rGiven);
        rClosest = rTrue(rClosestIndex);
        L = norm((rGiven  - rClosest), 'fro')^2 / norm(rClosest, 'fro')^2;
        % disp(L);
        
        % 9: if (L <= X_n and the kth packet passes the higher layer
        % authentication) then,
        if (L <= X_n) % && (passAuthenitication(mac, macGiven) == true)
            % 10: ˆ ri <-- ri
            if rem(vectorNum, 5) == 0
            [~, index] = max(abs(rTrue(:) - rClosest));
            rTrue(index) = rClaim;
            end
            
            % 11: Accept the kth packet
            %acceptPacket(k);
            
        % 12: else
        else
            % 13: Send spoofing alarm
            spooferNum(n) = spooferNum(n) + 1;
            % rSpooferVec = rGiven - rClosest;
        % 14: end if
        end
    % 15: end for
    end
    
    
    
    % 16: Observe sn+1 and Un
    
    [X_next] = max(Q(s_n_Index(1), s_n_Index(2), : ));
    thresholds(n) = X_next;
    [~, s_next] = observeState(X_next, channels(vectorNum, :), rSpooferVec);
    if s_n(1) > Pf_max
        Pf_max = s_n(1);
    end
    disp(s_n(1));

    % Adds onto expected sum utility
    U_n = U_n + getUtility(G0, G1, C0, C1, y, s_n(1), s_n(2));
    utilities(n) = U_n;
    
    % 17: Update Q(s_n,x_n) via(17) 
    Q(s_n_Index(1), s_n_Index(2), X_index) = (1 - learnRate) * Q(s_n_Index(1), s_n_Index(2), X_index) + learnRate*(U_n + discountFactor*V(s_next(1), s_next(2)));

    % 18: Update V (s_n) via(18)
    V(s_n_Index(1), s_n_Index(2)) = max(Q(s_n_Index(1), s_n_Index(2), : ));

    %19: end for
end

plot(spooferNum(:));
% disp("Number of spoofers: " + spoofers);
% pause(10);
% plot(thresholds);
% disp(utilities);
% [num, ~] = quantize(0.6, 100);
% disp("Num is: " + num);

function [s_n, s_n_Index] = observeState(X_n, rChanVec, rSpooferVec)%, numSpoofers)
% b, relative change in the channel gain due to 
    % environmental changes, assumed to be 0
    b = 0;
    
    % Power of received channel
    % channelPowerLegit = ;
    
    % Power of noise
    % noisePowerLegit = ;
    
    % SINR = ?, signal-to-interference-plus-noise-ratio (SINR) of packets 
    % sent by legitimate transmitter
    SINR =  0; % channelPowerLegit / noisePowerLegit;
    
    % o2, average power gain from legitimate transmitter at the receiver
    o2_f = norm(rChanVec)^2;
    o2_m = o2_f; % + (norm(rSpooferVec)^2);
    disp("o2_f: " + o2_f);
    disp("o2_m: " + o2_m);
    
    % ?, ratio of channel gain of spoofer to that of the legitimate
    % transmitter
    k =  0; % channelPowerReceived / channelPowerLegit;
    
    M = 64;
    
    s_quantNum = 499;
    s_n = [P_f(X_n, o2_f, SINR, b, M, s_quantNum), P_m(X_n, o2_m, SINR, k, M, s_quantNum)] * 100;
    s_n_Index = [0,0];
    s_n_range = 5;
    
    
    disp("P_f: " + s_n(1));
    disp("P_m: " + s_n(2));
    
    % Quantizing s_n
     [~, s_n_Index(1)] = quantize(s_n(1), s_quantNum, s_n_range);
     [~, s_n_Index(2)] = quantize(s_n(2), s_quantNum, s_n_range);
     % disp(s_n);
      disp(s_n_Index);
    % s_n_Index(:) = s_n(:);
end

function rM = readMac(packetK)

end

function gNP = getNumPackets()

end

function sSA = sendSpoofingAlarm()

end

function pA = passAuthenitication(mac, givenMac)
if mac == givenMac
    pA = true;
else
    pA = false;
end
end

function aP = acceptPacket(k)

end

function gU = getUtility(G0, G1, C0, C1, y1, Pf, Pm)
gU = ((G0 - G1) * y1 -  (G0 + C0) * Pm * y1) - ((G1 + C1) * Pf * (1 - y1)) + G1;
end

function Pf = P_f(x, o2, p, b, M, s_quantNum)
Pf = chi2cdf((x * 10e15)/ o2, 2*M); % 1 - chi2cdf(2*x*p / (2*o2 + b*p*o2), 2*M);
end

function Pm = P_m(x, o2, p, k, M, s_quantNum)
Pm = chi2cdf((x * 10e15) / o2, 2*M); % 1 - chi2cdf(2*x*p / (2*o2 + (1+k)*p*o2), 2*M);
end

function [quantizedValue, quantizedIndex] = quantize(value, quantNum, range)
edge = range/quantNum;
quantizedValue = 0;
quantizedIndex = 1;
 % value = value + 1;
for i = 0:(quantNum - 1)
    quant = edge * (quantNum - i);
    if (value >= quant)
        quantizedValue = quant;
        quantizedIndex = quantNum - i + 1;
        break
    end
end
end
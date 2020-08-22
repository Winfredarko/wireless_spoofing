clear;
close all;

fileID1 = fopen('TimeStamps10');
times = fread(fileID1, 'uint32');
fclose(fileID1);

fileID2 = fopen('MACs10');
macs = vec2mat(fread(fileID2, 'uint32'), 6);
fclose(fileID2);

fileID3 = fopen('ChannelEstimation10');
channelest = vec2mat(fread(fileID3, 'int32'), 512);
fclose(fileID3);

% vector of 512 elements, first 256 are real while latter 256 are the
% imaginary components

% try simple udp to figure out what is being sent
% 32 bit integers, presumably unsigned

% plotChannelEst(channelest, 3);

for i = 1:size(channelest,1)
    plotChannelEst(channelest, 100, i);
    title(['Channel Estimation ' num2str(i)]);
    pause(1);
    close;
end


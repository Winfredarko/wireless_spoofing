fileID1 = fopen('TimeStamps10');
times = fread(fileID1, 'uint32');
fclose(fileID1);

fileID2 = fopen('MACs10');
macs = vec2mat(fread(fileID2, 'uint32'), 6);
fclose(fileID2);

fileID3 = fopen('ChannelEstimation10');
channelest = fread(fileID3, 'uint8');
fclose(fileID3);
histogram(channelest(1:128));
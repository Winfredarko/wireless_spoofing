clear
clc
close all
hudpr = dsp.UDPReceiver('RemoteIPAddress','127.0.0.1', 'LocalIPPort', 25000, 'MessageDataType', 'uint8')
string=[];
count=1;
eof=0;
while(~eof)
    dataReceived = step(hudpr);
    if isempty( dataReceived )
        string=[];
        continue
    end
    if strcmp(char(dataReceived),'2')
        eof=1; %! signified the end
    else
        string=char(dataReceived); %fill the string
        count=count+1;
    end
    string
end
release(hudpr)
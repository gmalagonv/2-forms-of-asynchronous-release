%% 10.09.2020 function to genereate differnt IDs for each frame
function frameID = frameIDgenerator(locEv)
frameID = zeros(length(locEv),1);
for ev = 1:length(locEv)
    frameID(ev,1) =  str2double([num2str([locEv(ev).synID]), num2str([locEv(ev).frame])]);
    
end
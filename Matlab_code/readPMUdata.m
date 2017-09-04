function [bus_v,bus_freq ,name,ilf]=readPMUdata(path,allnames,num_bus) 
    name=fullfile(path,allnames{1,num_bus});
    data=load(name);
    bus_v= data.bus_v;
    bus_freq=data.bus_freq;
    bus_theta=data.theta;  
end

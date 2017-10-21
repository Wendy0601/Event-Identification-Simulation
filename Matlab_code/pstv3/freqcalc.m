% function [] = freqcalc(k,t,vmod )
% 
% global  basmva sys_freq bus_freq freqavg bus_freqf bf_hpf dbf_hpf o2p kpx %theta_un_k
% global  bus_v theta n_bus mac_con mac_spd
% 
% Nbf_hpf = 70;
% K1bf_hpf = 1;
%  
% if vmod == 0
%     kpx = ones(size(t));
%     bus_freq = zeros(n_bus+1,length(t));
%     bus_freq(:,1) = 1;% - to be CONFIRMED 
%     
%     freqavg = zeros(1,length(t));
%     %=====================================================================%
%     bf_hpf = zeros(size(theta));
%     dbf_hpf = zeros(size(theta));
%     bus_freqf = zeros(size(theta));
%     bus_freqf(:,1) = 1;
%     %bus_freqf = zeros(n_bus+1,length(t));
%     o2p = zeros(size(theta,1),1);
%     bf_hpf(:,1) = K1bf_hpf.*theta(:,1); % not sure it is here the initialization but just try
%     
%     %theta_un_k = zeros(size(theta));
% elseif vmod == 1
%     
%     if k~=1
% %         dVangle = angle(bus_v(1:n_bus,k)) - angle(bus_v(1:n_bus,k-1));
%         
%         dVangle = theta(:,k)-theta(:,k-1);
%         dVangle(dVangle >= pi) = dVangle(dVangle >= pi) - 2*pi;
%         dVangle(dVangle <= -pi) = dVangle(dVangle <= -pi) + 2*pi;
% 
%         dbus_freq = dVangle/(t(k)-t(k-1))*(1/(2*pi*sys_freq(1)));
% 
% %         dbus_freq(abs(dbus_freq) > 0.01) = 0;
%             
%         bus_freq(:,k) = bus_freq(:,1) + dbus_freq; 
%         
% %         bfr_d = bus_freq(:,k) - bus_freq(:,k-1);
% %         ix_dfr = abs(bfr_d) > 0.010;
% %         bus_freq(ix_dfr,k) = bus_freq(ix_dfr,1);
% 
% 
% 
%         %=====================================================================%
% %         if t(k)>3.95
% %             disp('st')
% %         end
%         
%         ixp2p = theta(:,k)-theta(:,k-1) >= pi;
%         ixm2p = theta(:,k)-theta(:,k-1) <= -pi;
%         o2p(ixp2p) = o2p(ixp2p) + 1*kpx(k);
%         o2p(ixm2p) = o2p(ixm2p) - 1*kpx(k);
%         theta_un_k = theta(:,k) - 2*pi.*o2p;
%                 
%         % theta_un = unwrap(theta,[],2);
%         % dbf_hpf(:,k) = Nbf_hpf.*K1bf_hpf.*theta_un_k(:,k) - Nbf_hpf.*bf_hpf(:,k);
%         dbf_hpf(:,k) = Nbf_hpf.*K1bf_hpf.*theta_un_k - Nbf_hpf.*bf_hpf(:,k);
%         bus_freqf(:,k) = bus_freqf(:,1) + dbf_hpf(:,k)*(1/(2*pi*sys_freq(1)));
%         kpx(k) = 0;
%     else
%         dbf_hpf(:,k) = Nbf_hpf.*K1bf_hpf.*theta(:,k) - Nbf_hpf.*bf_hpf(:,k);
%         bus_freqf(:,k) = bus_freqf(:,1) + dbf_hpf(:,k)*(1/(2*pi*sys_freq(1)));
%         % when k==1, the vector is already initialized
%     end
%     
%     TotH = sum(mac_con(:,16));
%     freqavg(:,k) = sum(mac_spd(:,k).*mac_con(:,16))/TotH;
%     
% end




function [] = freqcalc(k,t,vmod)

global  basmva sys_freq bus_freq freqavg bus_freqf bf_hpf dbf_hpf o2p kpx %theta_un_k
global  bus_v theta n_bus mac_con mac_spd
% definir el numero de buses como n_bus
% probar con theta
Nbf_hpf = 70;
K1bf_hpf = 1;

if vmod == 0
    kpx = ones(size(t));
    bus_freq = zeros(n_bus+1,length(t));
    bus_freq(:,1) = 1;% - to be CONFIRMED 
    
    freqavg = zeros(1,length(t));
    %=====================================================================%
    bf_hpf = zeros(size(theta));
    dbf_hpf = zeros(size(theta));
    bus_freqf = zeros(size(theta));
    bus_freqf(:,1) = 1;
    %bus_freqf = zeros(n_bus+1,length(t));
    o2p = zeros(size(theta,1),1);
    bf_hpf(:,1) = K1bf_hpf.*theta(:,1); % not sure it is here the initialization but just try
    
    %theta_un_k = zeros(size(theta));
elseif vmod == 1
    
    if k~=1        
        dVangle = theta(:,k)-theta(:,k-1);
        dVangle(dVangle >= pi) = dVangle(dVangle >= pi) - 2*pi;
        dVangle(dVangle <= -pi) = dVangle(dVangle <= -pi) + 2*pi;

        dbus_freq = dVangle/(t(k)-t(k-1))*(1/(2*pi*sys_freq(1)));
        if k==51
            dbus_freq = dVangle/(t(k)-t(k-1))*(1/(2*pi*60));%%%%%%%%%% it is some modification when event occurs when k=50
        end
            
        bus_freq(:,k) = bus_freq(:,1) + dbus_freq; 
        

        
        ixp2p = theta(:,k)-theta(:,k-1) >= pi;
        ixm2p = theta(:,k)-theta(:,k-1) <= -pi;
        o2p(ixp2p) = o2p(ixp2p) + 1*kpx(k);
        o2p(ixm2p) = o2p(ixm2p) - 1*kpx(k);
        theta_un_k = theta(:,k) - 2*pi.*o2p;
                
        % theta_un = unwrap(theta,[],2);
        % dbf_hpf(:,k) = Nbf_hpf.*K1bf_hpf.*theta_un_k(:,k) - Nbf_hpf.*bf_hpf(:,k);
        dbf_hpf(:,k) = Nbf_hpf.*K1bf_hpf.*theta_un_k - Nbf_hpf.*bf_hpf(:,k);
        bus_freqf(:,k) = bus_freqf(:,1) + dbf_hpf(:,k)*(1/(2*pi*sys_freq(1)));
        kpx(k) = 0;
    else
        dbf_hpf(:,k) = Nbf_hpf.*K1bf_hpf.*theta(:,k) - Nbf_hpf.*bf_hpf(:,k);
        bus_freqf(:,k) = bus_freqf(:,1) + dbf_hpf(:,k)*(1/(2*pi*sys_freq(1)));
        % when k==1, the vector is already initialized
    end
    
    TotH = sum(mac_con(:,16));
    freqavg(:,k) = sum(mac_spd(:,k).*mac_con(:,16))/TotH;
    
end

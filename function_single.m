%___________________________________________________________________%
%  PSO source codes version                                         %
%                                                                   %
%  Developed in MATLAB R2022b                                       %
%                                                                   %
%  Author and programmer: Dessalegn Bitew                           %
%                                                                   %
%         e-Mail: dessalegnbitew29@gmail.com                        %
%                 dessalegn_bitew@dmu.edu.et                        %
%                                                                   %
%       Homepage: http://www.alimirjalili.com                       %
%                                                                   %
%   Main paper: Aeggegn, Dessalegn Bitew, George Nyauma Nyakoe, and %
% Cyrus Wekesa. "Optimal sizing of grid connected multi-microgrid   %
% system using grey wolf optimization." Results in Engineering 23   %
% (2024): 102421.,                                                  %
%     DOI: https://doi.org/10.1016/j.rineng.2024.102421             %
%                                                                   %
%___________________________________________________________________%

function y = function_single(x)

 Npv = x(1);
 Nwt = x(2);
 Pelectro = x(3);
 Pfcu = x(4);
 Mtank = x(5);
 Ebatn = x(6);
 
persistent Cost_pan    Cost_el Cost_tank Cost_FC  Cost_bat Cost_elOM HHVh2
persistent  Rv  ir Cost_OMpv  Cost_OMVFC Kinv Kfc Cost_RFC Cost_Rinv Cost_OMtank
persistent p   DOD  Cost_batOM Vc Vn Vof Pr Cost_inv Cost_OMinv sell buy nef
persistent  Gmes  Ntemps   NOCT beta  Tj  Ta   Kb DOD1 Ve Cost_OMwind S Data Rpv
persistent  nfc nel   load_M   ndis nc ve Xe alpha ninv Cost_wind CRF
if isempty(Rv)
    %Financial aspect
  Cost_wind = 900;
    Cost_pan = 380  ;
    Cost_el = 1500/1e3 ;
    Cost_tank = 415.46 ;
    Cost_FC = 600/1e3;
    Cost_RFC = 600/1e3 ;
    Cost_inv = 400/1e3 ;
    Cost_Rinv = 400/1e3 ;
    Cost_OMinv = 10/1e3;
    Cost_OMwind = 2/1e3;
    Cost_elOM = 80/1e3;
    Cost_OMVFC = 0.08/1e3;
    Cost_OMpv = 4/1e3;
    Cost_bat = 203000/(1.67e6);
    Cost_batOM = 10/1e6;
    Cost_OMtank = 10;
    DOD = 0.9;
    DOD1 = 0.8;
    
    ndis = 0.95;
    nc = 0.95;
    ninv = 0.95;
    ir = 0.06;
   
    Rv = 20;
    
Kb = (1/(1+ir))^11; %+ (1/(1+ir))^20;
Kinv = (1/(1+ir))^16;
Kfc = (1/(1+ir))^11;
p = (( 1 + ir )^Rv -1)/(ir*(1+ir)^Rv);
CRF = (ir*(1+ir)^Rv)/(( 1 + ir )^Rv -1);

   %Sources data
    Data = xlsread('HOURLY DATA');
    Gmes = Data(:,2)';
    Ta =  Data(:,7)';% temperature
    
    sell = Data(:,6)';
    buy = Data(:,5)';
    ve = Data(:,3);
    alpha = 1/7;
    Xe = ((50/10).^alpha).*ve;
    Ve = Xe';
    HHVh2 = 39.4*1000;
    
    NOCT = 45;
    Ntemps=length(Gmes);
    
     beta =  0.41/100;
     Tj =(Ta+(Gmes*(NOCT-20)/(800)));
   
    load_M = 1e6*Data(:,4)';
    
     S = 1.976*0.991 ;
     Rpv = 0.1941 ;
     nef = Rpv.*(1 - beta.*(Tj - 25));
   
     %Storage data 
     nel = 0.85;
     nfc = 0.65;
     Pr = 1e3;
     %     Wind parameters
     Vc = 2;
     Vof = 40;
    
     Vn = 9; 
         
end
     
Pw = zeros(1,Ntemps);
ESurplus = zeros(1,Ntemps);
Hour = zeros(1,Ntemps);
Pfch2 = zeros(1,Ntemps);
Pel = zeros(1,Ntemps);

Ebuy = zeros(1,Ntemps);
Esell = zeros(1,Ntemps);
Etankmax = repmat(Mtank*HHVh2,1,Ntemps);
Etank = (1-DOD)*Etankmax; 
Etankmin = (1-DOD)*Etankmax; 
% Etank(:,1) = Etankmin;
  
Pdeficit = zeros(1,Ntemps);
Pfcg = zeros(1,Ntemps);

Pdis = zeros(1,Ntemps);
Pch = zeros(1,Ntemps);
R = zeros(1,Ntemps);

Ebatmax = repmat(Ebatn,1,Ntemps);
Ebatmin = (1-DOD1)*Ebatmax;
Ebat = Ebatmin;
Ebat(:,1) = Ebatn;
% Solar model
Pmodule = S.*Gmes.*nef;
Ppv = Npv.*Pmodule;

for i=1:1
    for j=1:Ntemps

        if Ve(i,j) < Vc || Ve(i,j) > Vof
            Pw(i,j) = 0;
        end
        if Ve(i,j)> Vc && Ve(i,j)< Vn
            
         Pw(i,j) = Nwt*Pr*((Ve(i,j)^3 -Vc^3)/(Vn^3-Vc^3));
         
        elseif Ve(i,j)>= Vn && Ve(i,j)< Vof
             
         Pw(i,j) = Nwt*Pr;
         
        end
        
        
    end
end

EB = Ppv + Pw - load_M/ninv;
for i=1:1
    for j=1:Ntemps

        if Ve(i,j) < Vc || Ve(i,j) > Vof
            Pw(i,j) = 0;
        end
        if Ve(i,j)> Vc && Ve(i,j)< Vn
            
         Pw(i,j) = Nwt*Pr*((Ve(i,j)^3 -Vc^3)/(Vn^3-Vc^3));
         
        elseif Ve(i,j)>= Vn && Ve(i,j)< Vof
             
         Pw(i,j) = Nwt*Pr;
         
        end
        
        
    end
end

EB = Ppv + Pw - load_M/ninv;
for i=1:1
       
    for j=2:Ntemps
       if EB(i,j) < 0
          
               
         if  (Etank(i,j-1)-Etankmin(i,j))*nfc == 0 || abs(EB(i,j)) < 5e6
             Pdis(i,j) = min((Ebat(i,j-1)- Ebatmin(i,j))*ndis,(abs(EB(i,j))));
             Ebat(i,j) = Ebat(i,j-1) - (1/ndis)*Pdis(i,j);
             Pfcg(i,j) = 0;
         elseif (Ebat(i,j-1)- Ebatmin(i,j))*ndis ==0
             Pfcg(i,j) = min((Etank(i,j-1)-Etankmin(i,j))*nfc,min(abs(EB(i,j)),Pfcu));
             Pfch2(i,j) = (1/nfc)*Pfcg(i,j);   
             Etank(i,j) = Etank(i,j-1)- Pfch2(i,j);
             Pdis(i,j) = 0;
         end
          if (Etank(i,j-1)-Etankmin(i,j))*nfc > 0 
              if(Ebat(i,j-1)- 0.8*Ebatmax(i,j))*ndis > 0
%          
             Pdis(i,j) = min((Ebat(i,j-1)- Ebatmin(i,j))*ndis,(abs(EB(i,j))));
             Ebat(i,j) = Ebat(i,j-1) - (1/ndis)*Pdis(i,j);
             Pfcg(i,j) = min((Etank(i,j-1)-Etankmin(i,j))*nfc,min(abs(EB(i,j))-Pdis(i,j),Pfcu));
             Pfch2(i,j) = (1/nfc)*Pfcg(i,j);   
             Etank(i,j) = Etank(i,j-1)- Pfch2(i,j);
%          elseif (Etank(i,j-1)-Etankmin(i,j))*nfc > 0 && (Ebat(i,j-1)- Ebatmin(i,j))*ndis > 0
%         
%               Pdis(i,j) = min((Ebat(i,j-1)- Ebatmin(i,j))*ndis,(abs(EB(i,j))- Pfcg(i,j)));
%               Ebat(i,j) = Ebat(i,j-1) - (1/ndis)*Pdis(i,j);
%              Pfcg(i,j) = min((Etank(i,j-1)-Etankmin(i,j))*nfc,min(abs(EB(i,j))-Pdis(i,j),Pfcu));
%              Pfch2(i,j) = (1/nfc)*Pfcg(i,j);   
%              Etank(i,j) = Etank(i,j-1)- Pfch2(i,j);
%          end
              else
        Pfcg(i,j) = min((Etank(i,j-1)-Etankmin(i,j))*nfc,min(abs(EB(i,j)),Pfcu));
        Pfch2(i,j) = (1/nfc)*Pfcg(i,j);   
        Etank(i,j) = Etank(i,j-1)- Pfch2(i,j);
        
        Pdis(i,j) = min((Ebat(i,j-1)- Ebatmin(i,j))*ndis,(abs(EB(i,j))- Pfcg(i,j)));
        Ebat(i,j) = Ebat(i,j-1) - (1/ndis)*Pdis(i,j);
               
        Pel(i,j) = 0; 
              end
          end

        Pdeficit(i,j) = abs(EB(i,j) + ninv*Pdis(i,j)+ ninv*Pfcg(i,j)) ;
        Ebuy(i,j) = buy(i,j)*Pdeficit(i,j)./1e6;
        
        
           if Pfcg(i,j)>0
               Hour(i,j)= 1;
           end     
                      
            
       else
       
        Pel(i,j) = min((Etankmax(i,j)-Etank(i,j-1))/nel,min(Pelectro,EB(i,j)));
        Etank(i,j) = Etank(i,j-1)+ nel*Pel(i,j);
        R(i,j) = EB(i,j)- Pel(i,j); 
        
        Pch(i,j) = min((Ebatmax(i,j)-Ebat(i,j-1))*(1/nc),(R(i,j)));
        Ebat(i,j) = Ebat(i,j-1) + nc*Pch(i,j);
        
        ESurplus(i,j) = R(i,j)- Pch(i,j);
        
        Hour(i,j)= 0;  
          
        Esell(i,j) = sell(i,j)*(ESurplus(i,j)/1e6);
        
        
       end
       if Ppv(i,j) == 0 && load_M(i,j)==0
           if Pw(i,j) == 0
        Etank(i,j) = Etank(i,j-1);
        Ebat(i,j) = Ebat(i,j-1);
           end
       end

    end
    
end
Pinvn = max(Ppv + Pw + Pfcg + Pdis);

lps = load_M - ninv*(Ppv + Pw + Pfcg +Pdis);
LPS= lps(lps>0);

 LPSP = (sum(LPS)/sum(load_M))*100;

Etotal = sum(load_M)/1000;

TNPCpv = Npv*Cost_pan+ Npv*Cost_OMpv*p;
TNPCfc = Pfcu*Cost_FC + Pfcu*Cost_RFC*Kfc + sum(Hour)*Pfcu*Cost_OMVFC*p  ;
TNPCB = Ebatn*Cost_bat + Ebatn*Cost_bat*Kb + Ebatn*Cost_batOM*p ;
TNPCt = Mtank*Cost_tank + Mtank*Cost_OMtank*p;
TNPCel = Pelectro*Cost_el + Pelectro*Cost_elOM*p;
TNPCinv = Pinvn*(Cost_inv + Cost_Rinv*Kinv + Cost_OMinv*p);
TNPCw = Nwt*(Cost_wind + Cost_OMwind*p);
Totselling = sum(Esell);
Totalbuying = sum(Ebuy);

 
Ysell = zeros(1,20);
Ybuy = zeros(1,20);

for i=1:20
Ysell(i) = Totselling*(1/((1+ir)^i));
Ybuy(i) = Totalbuying*(1/((1+ir)^i));
end

TNPC = TNPCpv + TNPCfc + TNPCt + TNPCel + TNPCB + TNPCw + TNPCinv  + sum(Ybuy) - sum(Ysell);

  penalty = zeros(1,2);
  c0 = [];
  c0(:,1) = LPSP-1;%-0.05

  PP = 10^20;
  
  for i=size(c0,1)
      for j=size(c0,2)
          if(c0(i,j)>0)
              penalty(i,j)=  PP*c0(i,j);
          end
      end
  end
y = TNPC + sum(penalty,2);

          
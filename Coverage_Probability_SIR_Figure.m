clear all;
close all;

e1 = 0.6;
e2 = 0.8;
sir_db = -10:20;
sir = 10.^(sir_db./10);


for i = 1 : length(sir_db);
    Ptotal=1;
    p11 = e1.*Ptotal;
    p12 = e2.*Ptotal;

    p21 = (1 - e1).*Ptotal;
    p22 = (1 - e2).*Ptotal; %1 - p1;  %(1 - e(i)).*Ptotal.*0.01;
    T1 = sir(i);
    T2 = sir(i);
        
        c = T1./(p11 - p21.*T1);
        fc = 1 - (1 ./ (1 + (sqrt(c.*Ptotal).*(pi./2 - atan(1 ./ sqrt(c.*Ptotal))))));
        fci = 2.*fc.*(1-fc)+fc.^2;
        fun(i) = 1 - fci;
        if(fun(i) <= 0)
            final(i) = 0;
        else
            final(i) = fun(i);
        end
        %function for coverage probability for user 1 at 0.6 e

   
        c1 = T1./(p12 - p22.*T2);
        fc1 = 1 - (1 ./ (1 + (sqrt(c1.*Ptotal).*(pi./2 - atan(1 ./ sqrt(c1.*Ptotal))))));
        fci1 = 2.*fc1.*(1-fc1)+fc1.^2;
        fun1(i) = 1 - fci1;
        if(fun1(i) <= 0)
            final1(i) = 0;
        else
            final1(i) = fun1(i);
        end
        %function for coverage probability for user 1 at 0.8 e
    
        fun2(i) = 1 - (1 - (1 ./ (1 + sqrt((T1./p21).*Ptotal).*(pi./2 - atan(1./sqrt((T1./p21).*Ptotal)))))).^2;
    %function for coverage probability for user 2 at 0.6 e
    
        fun21(i) = 1 - (1 - (1 ./ (1 + sqrt((T2./ p22).*Ptotal).*(pi./2 - atan(1./sqrt((T2./ p22).*Ptotal)))))).^2;
    %function for coverage probability for user 2 at 0.8 e
        
        co1 = T1./Ptotal;
        fco1 = 1 - (1 ./ (1 + (sqrt(co1.*Ptotal).*(pi./2 - atan(1 ./ sqrt(co1.*Ptotal))))));
        fcio1 = 2.*fco1.*(1-fco1)+fco1.^2;
        funo1(i) = 1 - fcio1;
        if(funo1(i) <= 0)
            finalo1(i) = 0;
        else
            finalo1(i) = funo1(i);
        end
        
        co2 = T1./Ptotal;
        fco2 = 1 - (1 ./ (1 + (sqrt(co2.*Ptotal).*(pi./2 - atan(1 ./ sqrt(co2.*Ptotal))))));
        fcio2 = fco2.^2;
        funo2(i) = 1 - fcio2;
        if(funo2(i) <= 0)
            finalo2(i) = 0;
        else
            finalo2(i) = funo2(i);
        end
end

for i = 1 : length(sir_db);
    Ptotal_1=1;
    p11_1 = e1.*Ptotal_1;
    p12_1 = e2.*Ptotal_1;

    p21_1 = (1 - e1).*Ptotal_1;
    p22_1 = (1 - e2).*Ptotal_1;
     for mk=1:10000 % Number of Monte Carlo Simulations
        
         T1_1 = sir(i);
         T2_1 = sir(i);
        
        c_1 = T1_1./(p11_1 - p21_1.*T1_1); 
        fc_1 = 1 - (1 ./ (1 + (sqrt(c_1.*Ptotal_1).*(pi./2 - atan(1 ./ sqrt(c_1.*Ptotal_1))))));
        fci_1 = 2.*fc_1.*(1-fc_1)+fc_1.^2;
        fun_1(i) = 1 - fci_1;
        if(fun_1(i) <= 0)
            final_1(i) = 0;
        else
            final_1(i) = fun_1(i);
        end
        
        c1_1 = T1_1./(p12_1 - p22_1.*T2_1);
        fc1_1 = 1 - (1 ./ (1 + (sqrt(c1_1.*Ptotal_1).*(pi./2 - atan(1 ./ sqrt(c1_1.*Ptotal_1))))));
        fci1_1 = 2.*fc1_1.*(1-fc1_1)+fc1_1.^2;
        fun1_1(i) = 1 - fci1_1;
        if(fun1_1(i) <= 0)
            final1_1(i) = 0;
        else
            final1_1(i) = fun1_1(i);
        end
        %function for coverage probability for user 1 at 0.8 e
    
        fun2_1(i) = 1 - (1 - (1 ./ (1 + sqrt((T1_1./p21_1).*Ptotal_1).*(pi./2 - atan(1./sqrt((T1_1./p21_1).*Ptotal_1)))))).^2;
        %function for coverage probability for user 2 at 0.6 e
    
        fun21_1(i) = 1 - (1 - (1 ./ (1 + sqrt((T2_1./ p22_1).*Ptotal_1).*(pi./2 - atan(1./sqrt((T2_1./ p22_1).*Ptotal_1)))))).^2;
        %function for coverage probability for user 2 at 0.8 e
     end
        plot1(i) = final_1(i) ;
        plot2(i) = final1_1(i) ;
        plot3(i) = fun2_1(i) ;
        plot4(i) = fun21_1(i);
   
end



plot(sir_db,final);
hold on;
plot(sir_db,final1);
plot(sir_db,fun2);
plot(sir_db,fun21);
 plot(sir_db,plot1,'^');
 hold on;
 plot(sir_db,plot2,'>');
 plot(sir_db,plot3,'o');
 plot(sir_db,plot4,'<');
plot(sir_db,finalo1,'Marker','x');
plot(sir_db,finalo2,'Marker','*');
legend('UE1 value e 0.6','UE1 value e 0.8','UE2 value e 0.6','UE2 value e 0.8','UE1 e 0.6 sim','UE1 e 0.8 sim','UE2 e 0.6 sim','UE2 e 0.8 sim','OFDMA UE1','OFDMA UE2','north','east');
xlabel('SIR (db)');
ylabel('coverage probability');
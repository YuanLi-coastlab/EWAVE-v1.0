%  calculate the breaking-induced wave energy dissipation rate

function [db] = wave_breaking_dissipation_functions (Hsw,Tsw,depth,mwl,gamma,kh,k,wave_type,Db_formula)
grav=9.81; pi=3.14159; rho=1000;

 if strcmpi(wave_type,'random') == 1
     if strcmpi(Db_formula,'BTG2020') == 1
            b=1;
            fp=1/Tsw;
            a=0.25*rho*grav*fp*b;  
            pfhb=0.88*tanh(gamma.*kh./0.88)./k;
            rhb=pfhb./Hsw;
            db=a.*Hsw.^2*(1+rhb.^2)*exp(-rhb.^2);  % based on the formula of Baldock 1998   
            if   rhb < 1  % use TG83 formula in case of very shallow water stability e.g., on the reef edge
                a=(3*sqrt(pi)*rho*grav*fp*b^3)./16;
                db=a*Hsw.^5./gamma.^2./(depth+mwl).^3*(1-1/(1+(Hsw./gamma./(depth+mwl)).^2).^2.5);
            end
            db=max(db,0);        
 
     elseif strcmpi(Db_formula,'B1998') == 1
            b=1;
            fp=1/Tsw;
            a=0.25*rho*grav*fp*b;  
            pfhb=0.88*tanh(gamma.*kh./0.88)./k;
            rhb=pfhb./Hsw;
            db=a.*Hsw.^2*(1+rhb.^2)*exp(-rhb.^2);
            db=max(db,0);
        elseif strcmpi(Db_formula,'JB2007') == 1
            b=1;
            fp=1/Tsw;
            a=0.25*rho*grav*fp*b;
            pfhb=0.88*tanh(gamma*kh./0.88)./k; 
            rhb=pfhb./Hsw;
            %%rhb=erf(rhb);
            db1=rhb.^3+1.5*rhb;
            db2=exp(-rhb.^2);
            db3=0.75*sqrt(pi)*(1-erf(rhb));
            db=a*Hsw.^3./(depth+mwl)*(db1*db2+db3);
            db=max(db,0);
        elseif strcmpi(Db_formula,'TG1983') == 1
            b=1;
            fp = 1/Tsw;
            a=(3*sqrt(pi)*rho*grav*fp*b^3)./16;
            db=a*Hsw.^5./gamma.^2./(depth+mwl).^3*(1-1/(1+(Hsw./gamma./(depth+mwl)).^2).^2.5);
            db=max(db,0);
        elseif strcmpi(Db_formula,'BJ1978') == 1
            b=1;
            fp=1/Tsw;
            a=0.25*rho*grav*fp*b;
            pfhb=0.88*tanh(gamma.*kh/0.88)./k;

            rhb=Hsw./pfhb;
            q=0;
            if (rhb==0);
            rhb=NaN;
            end
            for iter=1:999;         
            xxx=q;
            q=real(exp((q-1)./rhb.^2));
            q=0.5*xxx+0.5*q;
            if ((q-xxx)/xxx<0.01); 
            continue; end;
            end  

            db=a*pfhb.^2*q;
            db=max(db,0);
       
        elseif strcmpi(Db_formula,'0') == 1
            b=1;
            db=0;
     end
      

    elseif strcmpi(wave_type,'regular') == 1

        b=1;
        fp=1/Tsw;
        a=0.25*rho*grav*fp*b; 
        pfhb=0.88*tanh(gamma.*kh/0.88)./k;
        rhb=Hsw./pfhb;
        db=a*Hsw.^3/(depth+mwl)*(rhb>1);
     
end
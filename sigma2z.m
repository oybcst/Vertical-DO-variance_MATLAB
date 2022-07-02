function z = sigma2z(hc,sc_r,Cs_r,h,zeta,N)

for k=1:N
        z0=(hc.*sc_r(k)+Cs_r(k).*h)./(hc+h);
        z(k)=zeta+(zeta+h).*z0;
end
      
end
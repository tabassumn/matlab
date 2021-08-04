%all critical parameters
Tc = [535.5;591.8];
Pc = [41.5;41.1];
Vc = [267;316];
Zc = [0.249;0.264];
omega = [0.323;0.262];
n = input('no. of components: ');
P = 25; %in kPa
T = 323.15; %in K
%defining blank matrices
a = n*(n-1)/2;
Tcij = zeros(a,1);
Pcij = zeros(a,1);
Zcij = zeros(a,1);
Oij = zeros(a,1);
Vcij = zeros(a,1);
% ij parameters of critical values
k = 1;
for i = 1:n
    for j = 1:n
        if i~=j
         Tcij(k) = sqrt(Tc(i)*Tc(j));
         Oij(k) = (omega(i)+omega(j))*0.5;
         Zcij(k) = (Zc(i)+Zc(j))*0.5;
         Vcij(k) = (0.5*((Vc(i))^(1/3)+(Vc(j))^(1/3)))^3;
         Pcij(k) = Zcij(k) *83.14*Tcij(k)/Vcij(k); %R is 83.14
         k = k+1;
        else 
            continue
        end
    end
end
Tcij = unique(Tcij);
Pcij = unique(Pcij);
Zcij = unique(Zcij);
Oij = unique(Oij);
Vcij = unique(Vcij);
Tc = vertcat(Tc,Tcij);
Pc = vertcat(Pc,Pcij);
Vc = vertcat(Vc,Vcij);
Zc = vertcat(Zc,Zcij);
omega = vertcat(omega,Oij);
Trij = 323.15./Tc;
B_o = 0.083-0.422./(Trij).^1.6;
B_1 = 0.139-0.172./(Trij).^4.2;
B = 83.14*Tc.*(B_o+omega.*B_1)./(Pc);
Bij = zeros(n,n);
for i = 1:n
    Bij(i,i) = B(i);
    for j = 1:n
        if i==j
            continue
        else
            if (i+n)<=length(B)
            Bij(i,j) = B(i+n);
            Bij(j,i) = B(i+n);
            end
        end
    end
end
Bij
del_ik = zeros(a)
del_ij = zeros(a);
for k = 1:n
    for i = 1:n
        if i~=k
            del_ik(i,k) = (2*Bij(i,k)-Bij(i,i)-Bij(k,k));
            del_ik(k,i) = (2*Bij(i,k)-Bij(i,i)-Bij(k,k));
        else
            del_ik(i,k) = 0;
        end
        for j = 1:n
            if i~=j
             del_ij(i,j) = (2*Bij(i,j)-Bij(i,i)-Bij(j,j));
            else
                del_ij(i,j) = 0;
            end
        end
    end
end
del_ij
del_ik
sumnum = zeros(n,1);
y = zeros(n,1);
for m = 1:n-1
    X = ['specify ','mole fraction ', num2str(k)];
    disp(X)
    y(m) = input('mole fraction: ');
end
n
y(n) = 1-sum(y(1:n,:));
y
ln_phi = zeros(n,1);
for k = 1:n
    for i = 1:n
        for j = 1:n
            sumnum(k) = sumnum(k) + y(i)*y(j)*(2*del_ik(i,k)-del_ij(i,j));
            
        end
    end
    ln_phi(k) = (P/(8314*T))*(Bij(k,k)+0.5*sumnum(k));
end
phi = exp(ln_phi)
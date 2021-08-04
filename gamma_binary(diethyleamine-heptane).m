k = [1;2;33];
R = [0.9011;0.6744;1.2070];
Q = [0.848;0.540;0.936];
nu = [2 1 1;2 5 0];
sg = length(k);
n = input('components: ');
r = zeros(n,1);
for j = 1:n
    for l = 1:sg
        r(j) = r(j)+nu(j,l)*R(l);
    end
end
r;
q = zeros(n,1);
for j = 1:n
    for l = 1:sg
        q(j) = q(j)+nu(j,l)*Q(l);
    end
end
q;
e = zeros(sg,n)
for j = 1:n
    for l = 1:sg
        e(l,j) = Q(l)*nu(j,l)/q(j)
    end
end
e
a = zeros(length(k),length(k));
a(1:2,3) = 255.7;
a(3,1:2) = 65.33;
a
length(a)
tao = zeros(length(k),length(k))
for row_ind = 1: length(k)
    for col_ind = 1:length(k)
        tao(row_ind, col_ind) = exp(-a(row_ind,col_ind)/T);
    end
end
tao

beta = zeros(n,sg)
for i = 1:n
    for k = 1:sg
        for m = 1:sg
            beta(i,k) = beta(i,k)+e(m,i)*tao(m,k)
        end
    end
end
beta

%component mole fractions
z = zeros(1,n);
for i = 1:n-1
    str = ['component number ' num2str(i)];
    disp(str);
    z(i) = input('mole fraction i: ');
end
z(n) = 1-sum(z(1:n-1))
% define x as z
x = z;
%evaluating the numerator quantities in theta
theta_num = zeros(1,sg);
for k_ind = 1:sg
    for i = 1:n
        theta_num(k_ind) = theta_num(k_ind) + x(i)*q(i)*e(k_ind,i);
    end
end
theta_den = zeros(1,1);
for i = 1:n
    theta_den = theta_den + x(i)*q(i);
end
theta = theta_num./theta_den
s = zeros(1,sg);
for k_ind = 1:sg
    for m = 1:sg
        s(k_ind) = s(k_ind) + theta(m)*tao(m,k_ind);
    end
end
s
J = zeros(1,n);
for i = 1:n
    for j = 1:n
        J(i) = J(i) + r(j)*x(j)
    end
    J(i) = r(i)/J(i);
end
J;
L = zeros(1,n);
for i = 1:n
    for j = 1:n
        L(i) = L(i) + q(j)*x(j);
    end
    L(i) = q(i)/L(i);
end
L
ln_C = zeros(1,n);
ln_R = zeros(1,n);
for i = 1:n
    ln_C(i) = 1-J(i)+log(J(i))-5*q(i)*(1-(J(i)/L(i))+log((J(i)/L(i))));
end
for i = 1:n
    for k = 1:sg
        ln_R(i) = ln_R(i) + ((theta(k)*beta(i,k)/s(k))-(e(k,i)*log(beta(i,k)/s(k))));
    end
    ln_R(i) = q(i)*(1-ln_R(i));
end
ln_R
%evaluation of gamma
ln_gamma = zeros(1,n);
gamma = zeros(1,n);
for i = 1:n
    ln_gamma(i) = ln_R(i)+ln_C(i);
    gamma(i) = exp(ln_gamma(i));
end
gamma


        
    
        

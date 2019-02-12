function Quasilinearization_1D_PBE
clf; clear; format long; 
psi = 1;
x_0 = 1; x_N = 5;
% phi_0 = psi; phi_N = psi*exp(-4*x_N);
phi_0 = psi; phi_N = 2*atanh(tanh(psi/2)*exp(-x_N+x_0));
% phi_0 = psi;  phi_N = 1;
Err_tol = 10^(-12);
%
% exact solution
%
r = x_0:0.001:x_N;
% phi_exact = psi*exp(-4*r);
phi_exact = 2*atanh(tanh(psi/2)*exp(-r+x_0));
figure(1); clf(1); set(gca,'fontsize',24); hold on
plot(r,phi_exact); 
xlabel('x');
ylabel('\phi');
% string = sprintf('Err-tol = %d',Err_tol); title(string);
% 
% N=1 interval
% 
N = 1;
x_mid = [];
phi_mid = [];
plot_linear(1,x_mid,phi_mid,x_0,x_N,phi_0,phi_N);
% legend('exact solution','one-interval solution')
x = [x_0 x_N];
phi = [phi_0 phi_N];
% 
% N=2^m intervals
% 
for m = 1:3
    N = 2^m;
    %
    % project uniform phi points onto x-axis, initial guess
    % 
    if N == 2
        phi_unif = (phi_0 + phi_N)/2;
    elseif N >= 4
        d_phi = (phi_0 - phi_N)/N;
        phi_unif = phi_0-d_phi : -d_phi : phi_N+d_phi;
    end
    x_mid = zeros(1,N-1);
    for i = 1:N-1
        if mod(i,2) == 0
            x_mid(i) = x((i/2)+1);
        end
        if mod(i,2) == 1  
            for j = 1:(N/2)
                if(phi_unif(i)<=phi(j) && phi_unif(i)>=phi(j+1))
                    x_mid(i) = solve_x_new(x(j),x(j+1),phi(j),phi(j+1),phi_unif(i));
                end
            end
        end
    end        
    x = [x_0 x_mid x_N];
    
    %
    % outer iteration - solve for new x point
    %
    phi_mid = phi_unif;
    diff_x = x_N;
    while diff_x > Err_tol
    %
    % inner iteration - solve for new phi point
    %   
        if N == 2
            phi_mid = update_1_point_phi(x_mid,phi_unif,x_0,x_N,phi_0,phi_N,Err_tol);
        elseif N >=4
            phi_mid = update_phi(N,x_mid,phi_mid,x_0,x_N,phi_0,phi_N,Err_tol);
        end
    %
    % update the x points
    %
        phi = [phi_0 phi_mid phi_N];
        for i = 1:N-1
            for j = 1:N
                if(phi_unif(i)<=phi(j) && phi_unif(i)>=phi(j+1))
                    x_mid(i) = solve_x_new(x(j),x(j+1),phi(j),phi(j+1),phi_unif(i));
                end
            end
        end
        diff_x = norm(x_mid-x(2:N),inf);
        x = [x_0 x_mid x_N]; 
    end
    plot_linear(N,x_mid,phi_mid,x_0,x_N,phi_0,phi_N);
    phi = [phi_0 phi_mid phi_N];
    phi_exact = 2*atanh(tanh(psi/2)*exp(-x+x_0));
    table(m,1) = N;
    table(m,2) = norm(phi-phi_exact,inf);
    table(m,3) = table(m,2)*(N^2);
end
table
end


function y = f(phi)
% y = 16*phi; 
y = 0.5*sinh(phi*2);
% y = sinh(phi)/(1+2*0.5*sinh(phi/2)*sinh(phi/2));
end

function plot_linear(N,x_mid,phi_mid,x_0,x_N,phi_0,phi_N)
if N == 1
    dr = (x_N-x_0)/1000;
    r = x_0:dr:x_N;
    alpha = sqrt((f(phi_N) - f(phi_0))/(phi_N - phi_0));
    ui = sinh(alpha*(r-x_0))  /sinh(alpha*(x_N-x_0));
    vi = sinh(alpha*(x_N-r))/sinh(alpha*(x_N-x_0));
    beta = f(phi_0)/(alpha^2);
    phi_linear_homogeneous = phi_N*ui + phi_0*vi;
    phi_linear_particular = (phi_0 - beta)*(1-ui-vi);    
    phi_linear = phi_linear_particular + phi_linear_homogeneous;
    plot(r,phi_linear)
end
if N >= 2
    dr = (x_mid(1)-x_0)/1000;
    r = x_0:dr:x_mid(1);
    alpha(1) = sqrt((f(phi_mid(1)) - f(phi_0))/(phi_mid(1) - phi_0));
    ui = sinh(alpha(1)*(r-x_0))     /sinh(alpha(1)*(x_mid(1)-x_0));
    vi = sinh(alpha(1)*(x_mid(1)-r))/sinh(alpha(1)*(x_mid(1)-x_0));
    betai = f(phi_0)/(alpha(1)^2);
    phi_linear_homogeneous = phi_mid(1)*ui + phi_0*vi;
    phi_linear_particular = (phi_0 - betai)*(1-ui-vi);    
    phi_linear = phi_linear_particular + phi_linear_homogeneous;
    plot(r,phi_linear)
if N > 2; for i = 2:N-1
    dr = (x_mid(i)-x_mid(i-1))/1000;
    r = x_mid(i-1):dr:x_mid(i);
    alpha(i) = sqrt((f(phi_mid(i)) - f(phi_mid(i-1)))/(phi_mid(i) - phi_mid(i-1)));
    ui = sinh(alpha(i)*(r-x_mid(i-1)))/sinh(alpha(i)*(x_mid(i)-x_mid(i-1)));
    vi = sinh(alpha(i)*(x_mid(i)-r))  /sinh(alpha(i)*(x_mid(i)-x_mid(i-1)));
    betai = f(phi_mid(i-1))/(alpha(i)^2);
    phi_linear_homogeneous = phi_mid(i)*ui + phi_mid(i-1)*vi;
    phi_linear_particular = (phi_mid(i-1) - betai)*(1-ui-vi);    
    phi_linear = phi_linear_particular + phi_linear_homogeneous;
    plot(r,phi_linear)
end; end
    dr = (x_N-x_mid(N-1))/1000;
    r = x_mid(N-1):dr:x_N;
    alpha(N) = sqrt((f(phi_N) - f(phi_mid(N-1)))/(phi_N - phi_mid(N-1)));
    ui = sinh(alpha(N)*(r-x_mid(N-1)))/sinh(alpha(N)*(x_N-x_mid(N-1)));
    vi = sinh(alpha(N)*(x_N-r))       /sinh(alpha(N)*(x_N-x_mid(N-1)));
    betai = f(phi_mid(N-1))/(alpha(N)^2);
    phi_linear_homogeneous = phi_N*ui + phi_mid(N-1)*vi;
    phi_linear_particular = (phi_mid(N-1) - betai)*(1-ui-vi);    
    phi_linear = phi_linear_particular + phi_linear_homogeneous;
    plot(r,phi_linear)
    plot(x_mid,phi_mid,'o') 
end
end

function phi_mid = check_one_point(x_l,x_r,phi_l,phi_r,x_mid)
alpha = sqrt((f(phi_r) - f(phi_l))/(phi_r - phi_l));
ui = sinh(alpha*(x_mid-x_l))  /sinh(alpha*(x_r-x_l));
vi = sinh(alpha*(x_r-x_mid))/sinh(alpha*(x_r-x_l));
beta = f(phi_l)/(alpha^2);
phi_linear_homogeneous = phi_r*ui + phi_l*vi;
phi_linear_particular = (phi_l - beta)*(1-ui-vi);    
phi_mid = phi_linear_particular + phi_linear_homogeneous;
end

function x_mid = solve_x_new(x_l,x_r,phi_l,phi_r,phi_mid)
% 
% project the uniform phi points on to the current solution
% which depends on the updated phi,
% this function recieves the new x
%
alpha = sqrt((f(phi_l) - f(phi_r))/(phi_l- phi_r));
gamma = sinh(alpha*(x_r-x_l));
omega = cosh(alpha*(x_r-x_l));
betai = f(phi_l)/(alpha^2);

a = (betai*gamma - (betai*omega - betai - phi_r + phi_l))* ...
    (betai*gamma + (betai*omega - betai - phi_r + phi_l));
b = -2*(betai*omega - betai - phi_r + phi_l)* ...
             gamma*(phi_mid - phi_l + betai);
c = (betai*gamma - (gamma*(phi_mid - phi_l + betai)))* ...
    (betai*gamma + (gamma*(phi_mid - phi_l + betai)));

if b < 0
    X1 = (-b + sqrt(b^2 - 4*a*c))/(2*a);
    X2 = 2*c/(-b+sqrt(b^2 - 4*a*c));
end
if b > 0
    X1 = 2*c/(-b - sqrt(b^2 - 4*a*c));
    X2 = (-b - sqrt(b^2 - 4*a*c))/(2*a);
end

x1 = x_l + asinh(X1)/alpha;
x2 = x_l + asinh(X2)/alpha; 

if ((x_l - x1)<1e-15) && ((x1 - x_r)<1e-15)
    x_mid = x1;
elseif ((x_l - x2)<1e-15) && ((x2 - x_r)<1e-15)
    x_mid = x2;
end
phi_check = check_one_point(x_l,x_r,phi_l,phi_r,x_mid);
% [phi_mid phi_check]'
% plot(x_mid,phi_mid,'o'); 
end

function phi = update_phi(N,x,phi,x_0,x_N,phi_0,phi_N,Err_tol)
%
% This function updates the phi values for given x points.
%
    index = 0;
    phi = phi';
    diff_phi = 1;
    alpha = zeros(N,1);
    while diff_phi > Err_tol 
        index = index +1;
        
        alpha(1) = sqrt( (f(phi(1)) - f(phi_0))   /(phi(1) - phi_0)    );
        alpha(N) = sqrt( (f(phi_N) -  f(phi(N-1)))/(phi_N -  phi(N-1)) );
        for i = 2:N-1
            alpha(i) = sqrt( (f(phi(i)) - f(phi(i-1)))/(phi(i) - phi(i-1)) );
        end
        
        u_ip(1) = alpha(1)/tanh(alpha(1)*(x(1)-x_0));
        for i = 2:N-1
            u_ip(i) = alpha(i)/tanh(alpha(i)*(x(i)-x(i-1)));
        end
        
        v_ip(1) = -alpha(1)/sinh(alpha(1)*(x(1)-x_0));
        for i = 2:N-1
            v_ip(i) = -alpha(i)/sinh(alpha(i)*(x(i)-x(i-1)));
        end
        
        for i = 1:N-2
            u_iplus1p(i) = alpha(i+1)/sinh(alpha(i+1)*(x(i+1)-x(i)));
        end
        u_iplus1p(N-1) = alpha(N)/sinh(alpha(N)*(x_N-x(N-1))); 
        
        for i = 1:N-2
            v_iplus1p(i) = -alpha(i+1)/tanh(alpha(i+1)*(x(i+1)-x(i)));
        end
        v_iplus1p(N-1) = -alpha(N)/tanh(alpha(N)*(x_N-x(N-1)));
        
        a = v_ip;
        b = u_ip - v_iplus1p;
        c = -u_iplus1p;
        
        A = zeros(N-1);
        for i = 2: N-2
            A(i,i-1) = a(i);
            A(i,i) = b(i);
            A(i,i+1) = c(i);
        end
        A(1,1) = b(1);
        A(1,2) = c(1);
        A(N-1,N-2) = a(N-1);
        A(N-1,N-1) = b(N-1);

        beta(1) = f(phi_0)/(alpha(1)^2);   
        for i = 2:N
           beta(i) = f(phi(i-1))/(alpha(i)^2); 
        end
        
        rhs = zeros(N-1,1);
        for i = 2:N-1
            rhs(i) = (phi(i-1)-beta(i))*(u_ip(i) + v_ip(i)) - (phi(i) - beta(i+1))*(u_iplus1p(i) + v_iplus1p(i));
        end
        rhs(1) = (phi_0-beta(1))*(u_ip(1) + v_ip(1)) - (phi(1) - beta(2))*(u_iplus1p(1) + v_iplus1p(1)) - a(1)*phi_0;
        rhs(N-1) = rhs(N-1) - c(N-1)*phi_N;
        phi_new = A\rhs;
        diff_phi = norm(phi-phi_new,inf);
        phi = phi_new;
    end
    phi = phi';
end

function phi_mid = update_1_point_phi(x_mid,phi_mid,x_0,x_N,phi_0,phi_N,Err_tol)
diff_phi = 1;
while diff_phi > Err_tol
    alpha(1) = sqrt( (f(phi_mid)-f(phi_0))  /(phi_mid-phi_0) );
    alpha(2) = sqrt( (f(phi_N)  -f(phi_mid))/(phi_N-phi_mid) );
    a = -alpha(1)/tanh(alpha(1)*(x_mid-x_0));
    b =  alpha(1)/tanh(alpha(1)*(x_mid-x_0)) + ...
         alpha(2)/sinh(alpha(2)*(x_N-x_mid));
    c = -alpha(2)/sinh(alpha(2)*(x_N-x_mid));
    beta(1) = f(phi_0)/(alpha(1)^2); 
    beta(2) = f(phi_mid)/(alpha(2)^2);
    d = -alpha(1)/sinh(alpha(1)*(x_mid-x_0));
    e = -alpha(2)/tanh(alpha(2)*(x_N-x_mid));
    rhs = -beta(1)*(-a + d) + beta(2)*(-c + e);
    rhs = rhs - a*phi_0 - c*phi_N;
    phi_new = rhs/b;
    diff_phi = abs(phi_new-phi_mid);
    phi_mid = phi_new;
end
end

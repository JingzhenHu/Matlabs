function quasilinear_inf
clf; clear; format long; 
% coord = 0 means the 1D cartesian coordinate
% coord = 1 means the 2D polar coordinate
coord = 2;
psi = 2;
x_0 = 1; xplot_N = 10;
ratio = xplot_N/x_0;
phi_0 = psi; phi_N = 0;
Err_tol = 10^(-10);
% up to the 2^par interval of solution
par = 9; 
% use to record previous solutions - x values
% and store the derived projection x values 
sol_record = zeros(par,2^par-1); 
similar_sol_record_2D_x = zeros(par,2^par-1);
similar_sol_record_2D_phi = zeros(par,2^par-1);
list_color = ['b','r','k','c','m','g','r'];
% mid-phi points to project and test for error convergence
roof = 2^par;
if roof == 2
    phi_pro = (phi_0 + phi_N)/2;
elseif roof >= 4
    d_phi = (phi_0 - phi_N)/roof;
    phi_pro = phi_0-d_phi : -d_phi : phi_N+d_phi;
end
if coord == 0
    %
    % exact solution
    %
    r = x_0:0.001:xplot_N;
    phi_exact = 4*atanh(tanh(psi/4)*exp(-r+x_0));
    figure(1); clf(1); set(gca,'fontsize',24); hold on
    plot(r,phi_exact); 
    axis([x_0 xplot_N min(phi_N,phi_0) max(phi_0,phi_N)]);
    xlabel('x');
    ylabel('\phi');
end

if coord == 1 || coord == 2
    figure(1); clf(1); set(gca,'fontsize',24); hold on
    xlabel('x');
    ylabel('\phi');
    axis([x_0 xplot_N min(phi_N,phi_0) max(phi_0,phi_N)]);
end

% string = sprintf('Err-tol = %d',Err_tol); title(string);
% 
% N=1 interval
% 
N = 1;
x_mid = [];
phi_mid = [];
plot_linear(1,x_mid,phi_mid,x_0,xplot_N,phi_0,phi_N,coord,list_color(1));
if coord == 0
    legend('exact solution','one-interval solution')
elseif coord == 1
    legend('one-interval solution')
end
x = [x_0];
phi = [phi_0];

table2(1,1) = 1;
num = length(phi_pro);
one_intr_x = zeros(1,num);
for i = 1:num
    one_intr_x(i) = slove_x_last(coord,x_0,phi_0,phi_N,phi_pro(i),Err_tol,ratio);
end
sol_record(1,:) = one_intr_x;

if coord == 0
   exact_err_analysis = zeros(1,num);
   for i = 1:num
%        disp(phi_pro(i))
       exact_err_analysis(i) = err_proj(xplot_N,x_0,phi_0,phi_pro(i),Err_tol);
       % err_proj_newton(phi_0,x_0,xplot_N,phi_pro(i),Err_tol);
   end
end
% 
% N=2^m intervals
% 
for m = 1:par
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
            if N > 2
                for j = 1:(N/2)-1
                    if(phi_unif(i)<=phi(j) && phi_unif(i)>=phi(j+1))
                        % x_mid(i) = outter_newton(coord,x(j),x(j+1),phi(j),phi(j+1),phi_unif(i),Err_tol,0);
                        x_mid(i) = projection(coord,x(j),x(j+1),phi(j),phi(j+1),phi_unif(i),Err_tol,0);
                    end
                end
            end
            if (phi_unif(i)<=phi(N/2) && phi_unif(i)>=phi_N)
                x_mid(i) = slove_x_last(coord,x(N/2),phi(N/2),phi_N,phi_unif(i),Err_tol,ratio);
            end
        end
    end        
    x = [x_0 x_mid];
    
    %
    % outer iteration - solve for new x point
    %
    phi_mid = phi_unif;
    diff_x = xplot_N;
    while diff_x > Err_tol
        %
        % inner iteration - solve for new phi point
        %   
        if N == 2
            phi_mid = update_1_point_phi(coord,x_mid,phi_unif,x_0,phi_0,phi_N,Err_tol);
        elseif N >=4
            phi_mid = update_phi(coord,N,x_mid,phi_mid,x_0,phi_0,phi_N,Err_tol);
        end
        %
        % update the x points
        %
        phi = [phi_0 phi_mid];
        for i = 1:N-1
            for j = 1:N-1
                if(phi_unif(i)<=phi(j) && phi_unif(i)>=phi(j+1))
                    % x_mid(i) = outter_newton(coord,x(j),x(j+1),phi(j),phi(j+1),phi_unif(i),Err_tol,0);
                    x_mid(i) = projection(coord,x(j),x(j+1),phi(j),phi(j+1),phi_unif(i),Err_tol,0);
                end
            end
            if (phi_unif(i)<=phi(N) && phi_unif(i)>=phi_N)
                x_mid(i) = slove_x_last(coord,x(N),phi(N),phi_N,phi_unif(i),Err_tol,ratio);
            end
        end
        diff_x = norm(x_mid-x(2:N),inf);
        % disp(diff_x);
        x = [x_0 x_mid]; 
    end
    % print(list_color(m))
    plot_linear(N,x_mid,phi_mid,x_0,xplot_N,phi_0,phi_N,coord,list_color(mod(m,7)+1));
    phi = [phi_0 phi_mid];
    for indexxx = 1:length(x_mid)
        similar_sol_record_2D_x(m,indexxx) = x_mid(indexxx);
        similar_sol_record_2D_phi(m,indexxx) = phi_mid(indexxx);
    end
    % error analysis
    if m < par
        table2(m+1,1) = N;
    end
    if m < par
        up_pts = (2^par-1)-(N-1);
        ea = up_pts/N;
        index = 1;
        for j = 1:N-1
            for i = 1:ea
                test_x(index) = projection(coord,x(j),x(j+1),phi(j),phi(j+1),phi_pro(index),Err_tol,0);
                index = index + 1;
            end
            test_x(index) = x(j+1);
            index = index + 1;
        end
        for i = 1:ea
            test_x(index) = slove_x_last(coord,x(N),phi(N),phi_N,phi_pro(index),Err_tol,ratio);
            index = index + 1;
        end
        % disp(test_x)
        sol_record(m+1,:) = test_x;
    else
        sol_record(m+1,:) = x(2:2^par);
    end
    if coord == 1 && m < par
        table3(m,1) = N;
    end
    if coord == 2 && m < par
        table3(m,1) = N;
    end
    if coord == 0
        phi_exact = 4*atanh(tanh(psi/4)*exp(-x+x_0));
        table(m,1) = N;
        table(m,2) = norm(phi-phi_exact,inf); 
        % norm(exact_err_analysis-sol_record(m,:),inf);  
        table(m,3) = table(m,2)*(N^2);
    end
end
% disp(sol_record)
for i = 1:par
    table2(i,2) = max(abs(sol_record(i,:)-sol_record(par+1,:)));
    table2(i,3) = max(abs(sol_record(i,:)-sol_record(i+1,:)))*(table2(i,1)^2);
    table2(i,4) = table2(i,2)*(table2(i,1)^2);
end
if coord == 0
    table
end
table2
if coord == 1 || coord == 2
    for i = 1:par-1
        len = 2^i-1;
         phi_eva = similar_sol_record_2D_phi(par,:);
        for j = 1:len
            for k = 1:(2^par-2)
                left = similar_sol_record_2D_x(par,k);
                right = similar_sol_record_2D_x(par,k+1);
                phi_l = similar_sol_record_2D_phi(par,k);
                phi_r = similar_sol_record_2D_phi(par,k+1);
                if(similar_sol_record_2D_x(i,j) >= left && similar_sol_record_2D_x(i,j) <= right)
                    phi_eva(j*2^(par-i)) = sol(coord,similar_sol_record_2D_x(i,j),left,right,phi_l,phi_r);
                end
            end
        end
        table3(i,2) = max(abs(phi_eva-similar_sol_record_2D_phi(par,:)));
        table3(i,3) = table3(i,2)*(table3(i,1)^2);
    end
    if par > 1
        table3
    end
end
end

function evaluation = sol(coord,x,x_l,x_r,phi_l,phi_r)
alpha = sqrt((f(phi_r) - f(phi_l))/(phi_r - phi_l));
beta = f(phi_l)/(alpha^2);
ui = u(coord,alpha,x,x_l,x_r);
vi = v(coord,alpha,x,x_l,x_r);
evaluation = phi_r*ui + phi_l*vi + (phi_l - beta)*(1-ui-vi);
end

function x = slove_x_last(coord,x_l,phi_l,phi_N,phi_last,Err_tol,ratio)
% here assume x_0 and x_plotN are both positive since for coordinate, 
% the origin is a singularity
x_r = ratio*x_l/2;
x = false;
while x == false
    x = projection(coord,x_l,x_r,phi_l,phi_N,phi_last,Err_tol,1);
    x_r = 2*x_r;
end
end

function y = f(phi)
% y = 16*phi; 
y = sinh(phi);
end

function res = u(coord,alpha,x,x_l,x_r)
if coord == 0
    res = sinh(alpha*(x-x_l))  /sinh(alpha*(x_r-x_l));
end
if coord == 1
    Q = besseli(0,alpha*x_l)*besselk(0,alpha*x_r) - ...
        besseli(0,alpha*x_r)*besselk(0,alpha*x_l);
    res =(-besseli(0,alpha*x)*besselk(0,alpha*x_l) ...
          +besseli(0,alpha*x_l)*besselk(0,alpha*x))/Q; 
end
if coord == 2
    res = x_r*sinh(alpha*(x-x_l))./ (x*sinh(alpha*(x_r-x_l)));
end
end

function res = v(coord,alpha,x,x_l,x_r)
if coord == 0
    res = sinh(alpha*(x_r-x))  /sinh(alpha*(x_r-x_l));
end
if coord == 1
    Q = besseli(0,alpha*x_l)*besselk(0,alpha*x_r) - ...
        besseli(0,alpha*x_r)*besselk(0,alpha*x_l);
    res =(besseli(0,alpha*x)*besselk(0,alpha*x_r) ...
         -besseli(0,alpha*x_r)*besselk(0,alpha*x))/Q;
end
if coord == 2
    res = x_l*sinh(alpha*(x_r-x))./ (x*sinh(alpha*(x_r-x_l)));
end
end

function res = up(coord,alpha,x,x_l,x_r)
if coord == 0
    res = alpha*cosh(alpha*(x - x_l))/ ...
                sinh(alpha*(x_r-x_l));
end
if coord == 1
    Q = besseli(0,alpha*x_l)*besselk(0,alpha*x_r) - ...
        besseli(0,alpha*x_r)*besselk(0,alpha*x_l);
    res =-alpha*(besseli(1,alpha*x)*besselk(0,alpha*x_l) ...
               +besseli(0,alpha*x_l)*besselk(1,alpha*x))/Q; 
end
if coord == 2
    res = x_r*alpha*cosh(alpha*(x-x_l))./ (x*sinh(alpha*(x_r-x_l))) ...
         -x_r*sinh(alpha*(x-x_l))./ (x^2*sinh(alpha*(x_r-x_l)));
end
end

function res = vp(coord,alpha,x,x_l,x_r)
if coord == 0
    res = -alpha*cosh(alpha*(x_r - x))/ ...
                 sinh(alpha*(x_r-x_l));
end
if coord == 1
    Q = besseli(0,alpha*x_l)*besselk(0,alpha*x_r) - ...
        besseli(0,alpha*x_r)*besselk(0,alpha*x_l);
    res =alpha*(besseli(1,alpha*x)*besselk(0,alpha*x_r) ...
               +besseli(0,alpha*x_r)*besselk(1,alpha*x))/Q;
end
if coord == 2
    res = -x_l*alpha*cosh(alpha*(x_r-x))./ (x*sinh(alpha*(x_r-x_l))) ...
          -x_l*sinh(alpha*(x_r-x))./(x^2*sinh(alpha*(x_r-x_l)));
end
end

function plot_linear(N,x_mid,phi_mid,x_0,xplot_N,phi_0,phi_N,coord,color)
if N == 1
    dr = (xplot_N-x_0)/1000;
    r = x_0:dr:xplot_N;
    alpha = sqrt((f(phi_N) - f(phi_0))/(phi_N - phi_0)); 
    if coord == 0
        phi_linear = phi_0*exp(alpha*(x_0-r));
    end
    if coord == 1
        phi_linear = phi_0*besselk(0,alpha*r)/besselk(0,alpha*x_0);
    end
    if coord == 2
        phi_linear = phi_0*x_0*exp(alpha*(x_0-r))./r;
    end
    plot(r,phi_linear,color)
end
if N >= 2
    dr = (x_mid(1)-x_0)/1000;
    r = x_0:dr:x_mid(1);
    alpha(1) = sqrt((f(phi_mid(1)) - f(phi_0))/(phi_mid(1) - phi_0));
    ui = u(coord,alpha(1),r,x_0,x_mid(1));
    vi = v(coord,alpha(1),r,x_0,x_mid(1));
    betai = f(phi_0)/(alpha(1)^2);
    phi_linear_homogeneous = phi_mid(1)*ui + phi_0*vi;
    phi_linear_particular = (phi_0 - betai)*(1-ui-vi);    
    phi_linear = phi_linear_particular + phi_linear_homogeneous;
    plot(r,phi_linear,color)
if N > 2; for i = 2:N-1
    dr = (x_mid(i)-x_mid(i-1))/1000;
    r = x_mid(i-1):dr:x_mid(i);
    alpha(i) = sqrt((f(phi_mid(i)) - f(phi_mid(i-1)))/(phi_mid(i) - phi_mid(i-1)));
    ui = u(coord,alpha(i),r,x_mid(i-1),x_mid(i));
    vi = v(coord,alpha(i),r,x_mid(i-1),x_mid(i));
    betai = f(phi_mid(i-1))/(alpha(i)^2);
    phi_linear_homogeneous = phi_mid(i)*ui + phi_mid(i-1)*vi;
    phi_linear_particular = (phi_mid(i-1) - betai)*(1-ui-vi);    
    phi_linear = phi_linear_particular + phi_linear_homogeneous;
    plot(r,phi_linear,color)
end; end
    dr = (xplot_N-x_mid(N-1))/1000;
    r = x_mid(N-1):dr:xplot_N;
    alpha(N) = sqrt(f(phi_mid(N-1))/phi_mid(N-1));
    if coord == 0
        phi_linear = phi_mid(N-1)*exp(alpha(N)*(x_mid(N-1)-r));
    end
    if coord == 1
%         disp('coord = 1 case for 2-interval')
%         disp(phi_mid(N-1)/besselk(0,alpha(N)*x_mid(N-1)))
        coeff = phi_mid(N-1)/besselk(0,alpha(N)*x_mid(N-1));
        phi_linear = coeff*besselk(0,alpha(N)*r);
    end
    if coord == 2
         phi_linear = phi_mid(N-1)*x_mid(N-1)*exp(alpha(N)*(x_mid(N-1)-r))./r;
    end
    plot(r,phi_linear,color)
    plot(x_mid,phi_mid,'o') 
end
end


function target = s(x,coord,x_l,x_r,phi_l,phi_r,phi_mid,last_int)
alpha = sqrt((f(phi_l) - f(phi_r))/(phi_l- phi_r));
if last_int == 0
    beta = f(phi_l)/(alpha^2);
    ui = u(coord,alpha,x,x_l,x_r);
    vi = v(coord,alpha,x,x_l,x_r);
    target = phi_r*ui + phi_l*vi + (phi_l-beta)*(1-ui-vi) - phi_mid;
end
if last_int == 1
    if coord == 0
        target = phi_l*exp(-alpha*(x-x_l)) - phi_mid;
    end
    if coord == 1
        target = phi_l*besselk(0,alpha*x)/besselk(0,alpha*x_l) - phi_mid;
    end
    if coord == 2
        target = phi_l*x_l*exp(alpha*(x_l-x))./x - phi_mid;
    end
end
end

function x_mid = projection(coord,x_l,x_r,phi_l,phi_r,phi_mid,Err_tol,last_int)
if s(x_l,coord,x_l,x_r,phi_l,phi_r,phi_mid,last_int)...
  *s(x_r,coord,x_l,x_r,phi_l,phi_r,phi_mid,last_int) > 0
    disp(s(x_l,coord,x_l,x_r,phi_l,phi_r,phi_mid,last_int))
    disp(s(x_r,coord,x_l,x_r,phi_l,phi_r,phi_mid,last_int))
    disp('bisection cannot work for this problem!')
    x_mid = false;
else
    a = x_l;
    b = x_r;
    x_mid = (a + b)/2;
    err = abs(s(x_mid,coord,x_l,x_r,phi_l,phi_r,phi_mid,last_int));
    iter_num = 0;
    while err > (Err_tol*1e-2)
        left = s(a,coord,x_l,x_r,phi_l,phi_r,phi_mid,last_int);
%         x_mid
        mid = s(x_mid,coord,x_l,x_r,phi_l,phi_r,phi_mid,last_int);
%         b
%         right = s(b,coord,x_l,x_r,phi_l,phi_r,phi_mid,last_int)
       if (left <= 0 && mid >= 0) ||(left >= 0 && mid <= 0)
           b = x_mid;
       else%if mid*right < 0
           a = x_mid;
       end
       x_mid = (a+b)/2;
       % disp(x_mid);
       iter_num = iter_num+1;
       err = abs(s(x_mid,coord,x_l,x_r,phi_l,phi_r,phi_mid,last_int));
    end
%     disp(iter_num)
end
end

function target = sp(x,coord,x_l,x_r,phi_l,phi_r,last_int)
alpha = sqrt((f(phi_l) - f(phi_r))/(phi_l- phi_r));
if last_int == 0
    beta = f(phi_l)/(alpha^2);
    uip = up(coord,alpha,x,x_l,x_r);
    vip = vp(coord,alpha,x,x_l,x_r);
    target = phi_r*uip + phi_l*vip - (phi_l-beta)*(uip+vip);
end
if last_int == 1
    if coord == 0
        target = -alpha*phi_l*exp(-alpha*(x-x_l));
    end
    if coord == 1
        target = -alpha*phi_l*besselk(1,alpha*x)/besselk(0,alpha*x_l);
    end
end
end

function err_proj_1D = e(x,x_0,phi_0,phi_mid)
err_proj_1D = 4*atanh(tanh(phi_0/4)*exp(-x+x_0))-phi_mid;
end

function x = err_proj(x_N,x_0,phi_0,phi_mid,Err_tol)
while e(x_0,x_0,phi_0,phi_mid)...
  *e(x_N,x_0,phi_0,phi_mid) > 0
    disp(e(x_0,x_0,phi_0,phi_mid))
    disp(e(x_N,x_0,phi_0,phi_mid))
    disp('extend bisection right bound to make it works!');
    x_N = 2*x_N;
end
a = x_0;
b = x_N;
x = (a + b)/2;
err = abs(e(x,x_0,phi_0,phi_mid));
iter_num = 0;
while err > (Err_tol)
   left = e(x_0,x_0,phi_0,phi_mid);
   mid = e(x,x_0,phi_0,phi_mid);
   if (left <= 0 && mid >= 0) ||(left >= 0 && mid <= 0)
       b = x;
   else%if mid*right < 0
       a = x;
   end
   x = (a+b)/2;
   % disp(x)
   iter_num = iter_num+1;
   err = abs(e(x,x_0,phi_0,phi_mid));
end
end

function err_proj_1D_1_der = ep(phi)
err_proj_1D_1_der = -2*sinh(phi/2);
end

function x = err_proj_newton(phi_0,x_0,x_plotN,phi,Err_tol)
itmax = 100;
del=2*Err_tol;
k=0;
% initial guess - x_plotN
x= x_plotN;
y=e(phi_0,x_0,x,phi);

while (abs(del) > Err_tol) && (k < itmax)
    del=-y/ep(phi);
    fold=abs(y);
    xold=x;
    alph=1;
    fnew=2*fold;
    while (fnew >= fold) && (k < itmax)
        x=xold+alph*del;
        y=e(phi_0,x_0,x,phi);
        fnew=abs(y);
        alph=alph/2;
        k=k+1;
    end 
end
if itmax == 100
    disp('Newton method might not be convergent!')
end
end

function phi = update_phi(coord,N,x,phi,x_0,phi_0,phi_N,Err_tol)
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
        u_ip(1) = up(coord,alpha(1),x(1),x_0,x(1));
        for i = 2:N-2
            u_ip(i) = up(coord,alpha(i),x(i),x(i-1),x(i));
        end
        
        v_ip(1) = vp(coord,alpha(1),x(1),x_0,x(1));
        for i = 2:N-2
            v_ip(i) = vp(coord,alpha(i),x(i),x(i-1),x(i));
        end
        
        for i = 1:N-2
            u_iplus1p(i) = up(coord,alpha(i+1),x(i),x(i),x(i+1));
        end
        
        for i = 1:N-2
            v_iplus1p(i) = vp(coord,alpha(i+1),x(i),x(i),x(i+1));
        end
        
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
        u_Nmins1p = up(coord,alpha(N-1),x(N-1),x(N-2),x(N-1));
        A(N-1,N-2) = -u_Nmins1p;
        if coord == 0
            A(N-1,N-1) = u_Nmins1p + alpha(N);
        end
        if coord == 1
            A(N-1,N-1) = u_Nmins1p + ...
                alpha(N)*besselk(1,alpha(N)*x(N-1))/besselk(0,alpha(N)*x(N-1));
        end
        if coord == 2
            A(N-1,N-1) = u_Nmins1p + alpha(N) + 1/x(N-1);
        end

        beta(1) = f(phi_0)/(alpha(1)^2);   
        for i = 2:N
           beta(i) = f(phi(i-1))/(alpha(i)^2); 
        end
        
        rhs = zeros(N-1,1);
        for i = 2:N-2
            rhs(i) = (phi(i-1)- beta(i)  )*(u_ip(i)      + v_ip(i)) - ...
                     (phi(i)  - beta(i+1))*(u_iplus1p(i) + v_iplus1p(i));
        end
        rhs(1) = (phi_0  - beta(1))*(u_ip(1)      + v_ip(1)) - ...
                 (phi(1) - beta(2))*(u_iplus1p(1) + v_iplus1p(1)) - a(1)*phi_0;
        rhs(N-1) = -beta(N-1)*(u_Nmins1p + vp(coord,alpha(N-1),x(N-1),x(N-2),x(N-1)));
        phi_new = A\rhs;
        diff_phi = norm(phi-phi_new,inf);
        phi = phi_new;
    end
    phi = phi';
    % check for accuracy 
    if coord == 0
        for i = 1:N-1
            if i == 1
                u_1p = up(coord,alpha(i),x(i),x_0,x(i));
                v_1p = vp(coord,alpha(i),x(i),x_0,x(i));
                u_2p = up(coord,alpha(i+1),x(i),x(i),x(i+1));
                v_2p = vp(coord,alpha(i+1),x(i),x(i),x(i+1));
                dev_pre = phi(i)*u_1p + phi_0*v_1p - (phi_0 - beta(i))*(u_1p + v_1p);
                dev_pos = phi(i)*u_2p + phi(i+1)*v_2p - (phi(i) - beta(i+1))*(u_2p + v_2p);
                if(dev_pre - dev_pos > Err_tol)
                    disp(dev_pre)
                    disp(dev_pos)
                    disp('fail')
                end
            end
            if (i > 1) && (i< (N-1))
                u_1p = up(coord,alpha(i),x(i),x(i-1),x(i));
                v_1p = vp(coord,alpha(i),x(i),x(i-1),x(i));
                u_2p = up(coord,alpha(i+1),x(i),x(i),x(i+1));
                v_2p = vp(coord,alpha(i+1),x(i),x(i),x(i+1));
                dev_pre = phi(i)*u_1p + phi(i-1)*v_1p - (phi(i-1) - beta(i))*(u_1p + v_1p);
                dev_pos = phi(i)*u_2p + phi(i+1)*v_2p - (phi(i) - beta(i+1))*(u_2p + v_2p);
                if(dev_pre - dev_pos > Err_tol)
                    disp(dev_pre)
                    disp(dev_pos)
                    disp('fail')
                end
            end
            if i == N-1
                u_1p = up(coord,alpha(i),x(i),x(i-1),x(i));
                v_1p = vp(coord,alpha(i),x(i),x(i-1),x(i));
                dev_pre = phi(i)*u_1p + phi(i-1)*v_1p - (phi(i-1) - beta(i))*(u_1p + v_1p);
                dev_pos = -alpha(i+1)*phi(i);
                if(dev_pre - dev_pos > Err_tol)
                    disp(dev_pre)
                    disp(dev_pos)
                    disp('fail')
                end
            end
        end
    end
end

function phi_mid = update_1_point_phi(coord,x_mid,phi_mid,x_0,phi_0,phi_N,Err_tol)
diff_phi = 1;
if coord == 0
    while diff_phi > Err_tol
        % slope for the linearized right hand side function in each interval
        alpha(1) = sqrt( (f(phi_mid)-f(phi_0))  /(phi_mid-phi_0) );
        alpha(2) = sqrt( (f(phi_N)  -f(phi_mid))/(phi_N-phi_mid) );
        u_1p = up(coord,alpha(1),x_mid,x_0,x_mid);
        v_1p = vp(coord,alpha(1),x_mid,x_0,x_mid);
        b =  u_1p + alpha(2);
        beta = f(phi_0)/(alpha(1)^2);
        rhs = -beta*(u_1p + v_1p);
        rhs = rhs + u_1p*phi_0;
        phi_new = rhs/b;
        diff_phi = abs(phi_new-phi_mid);
        phi_mid = phi_new;
    end
    % check the derivatives:
    der_1 = phi_mid*u_1p + phi_0*v_1p - (phi_0 - beta)*(u_1p + v_1p);
    der_2 = -alpha(2)*phi_mid;
    if der_1 - der_2 > Err_tol
        disp('derivative of the first: ')
        disp(der_1)
        disp('derivative of the second: ')
        disp(der_2)
    end
end
if coord == 1
    while diff_phi > Err_tol
        % slope for the linearized right hand side function in each interval
        alpha(1) = sqrt( (f(phi_mid)-f(phi_0))  /(phi_mid-phi_0) );
        alpha(2) = sqrt( (f(phi_N)  -f(phi_mid))/(phi_N-phi_mid) );
        u_1p = up(coord,alpha(1),x_mid,x_0,x_mid);
        v_1p = vp(coord,alpha(1),x_mid,x_0,x_mid);
        b =  u_1p + alpha(2)*besselk(1,alpha(2)*x_mid)/besselk(0,alpha(2)*x_mid);
        beta = f(phi_0)/(alpha(1)^2);
        rhs = -beta*(u_1p + v_1p);
        rhs = rhs + u_1p*phi_0;
        phi_new = rhs/b;
        diff_phi = abs(phi_new-phi_mid);
        phi_mid = phi_new;
    end
    % check the derivatives:
    der_1 = phi_mid*u_1p + phi_0*v_1p - (phi_0 - beta)*(u_1p + v_1p);
    der_2 = -alpha(2)*phi_mid*besselk(1,alpha(2)*x_mid)/besselk(0,alpha(2)*x_mid);
    if der_1 - der_2 > Err_tol
        disp('derivative of the first: ')
        disp(der_1)
        disp('derivative of the second: ')
        disp(der_2)
    end
end
if coord == 2
    while diff_phi > Err_tol
        % slope for the linearized right hand side function in each interval
        alpha(1) = sqrt( (f(phi_mid)-f(phi_0))  /(phi_mid-phi_0) );
        alpha(2) = sqrt( (f(phi_N)  -f(phi_mid))/(phi_N-phi_mid) );
        u_1p = up(coord,alpha(1),x_mid,x_0,x_mid);
        v_1p = vp(coord,alpha(1),x_mid,x_0,x_mid);
        b =  u_1p + alpha(2) + 1/x_mid;
        beta = f(phi_0)/(alpha(1)^2);
        rhs = -beta*(u_1p + v_1p);
        rhs = rhs + u_1p*phi_0;
        phi_new = rhs/b;
        diff_phi = abs(phi_new-phi_mid);
        phi_mid = phi_new;
    end
end
end


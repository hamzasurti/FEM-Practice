function myfem 

clc
clear
close all

% set up double precision

format long e

R=1000;
A= 10^-4;
E = 70000000000;

f = inline('-x*1000/(70000000000*(10^-4))');
a = 0;
b = 1;
g = 0;
h = R/(E*A);



numelm = 10;

[x, u] = fem1d_case1(f,g,h,a,b,numelm);

plot(x,u,'-ro','LineWidth',4)
title('Finite Element Result')
xlabel('x')
ylabel('u')
hold on
%fun = inline(dsolve('D2y = -x', 'Dy(0) = 0', 'y(1) = 1','x'));
fun = inline(dsolve(['D2y = ' char(f)], ['Dy(1) = ' num2str(h)], ['y(0) =' num2str(g)] ,'x'));
u_exact = fun(x);
plot(x,u_exact,'-b','LineWidth',2)
hold off
legend('FEM Solution', ['Exact Solution ' char(fun)],0);

function [x, u] = fem1d_case1(f,g,h,a,b,numelm)
%
%   solve the following ODE using FEM
%
%   u,xx = f,   on [a,b]  
%   u(b) = g                    
%   u,x(a) = h                
%
%   input:
%   f:  real function
%   g:  real number
%   h:  real number
%   a:  real number
%   b:  real number
%   n:  number of element
%
%   output:
%   u:  solution
%   x:  node coordinates
%
%   check input data

    if(a>=b)
        error(' parameter ''a'' needs to be smaller than parameter b')
    end
    
    if(numelm <= 0)
        error('number of element needs to be greater than zero')
    end

% construct mesh data

    % set number of nodes

    numnod = numelm + 1;
    
    % set number of equations
    
    numeqs = numelm;
    
    % set node coordinates
    
    x = linspace(a,b,numnod)';
    
    % set mesh connectivity array
    
    en = zeros(numelm,2);
    for n = 1:numelm
        en(n,:) = [n n+1];
    end
    
    % set location matrix
    
    lm = zeros(numelm,2);
    for n = 1:numelm
        lm(n,:) = [n-1 n];
    end
   % lm(numelm,1) = numelm;
    
% assemble system matrix and vector

    sysmat = zeros(numeqs);
    sysvec = zeros(numeqs,1);
    
    for n = 1 : numelm
        
        % element coordinates
        
        xe(1) = x(en(n,1));
        xe(2) = x(en(n,2));
        
        % element quadrature coordinate and weight
        
        xq = (xe(1)+xe(2))/2;
        wq = xe(2)-xe(1);
        
        % values of shape functions and their derivatives at quadrature 
        
        shp(1) = 0.5;
        shp(2) = 0.5;
        dshp(1) =-1/(xe(2)-xe(1));
        dshp(2) = 1/(xe(2)-xe(1));
        
        % element matrix
        
        elmmat(1,1) = dshp(1)*dshp(1)*wq;
	    elmmat(2,2) = dshp(2)*dshp(2)*wq;
	    elmmat(1,2) = dshp(1)*dshp(2)*wq;
	    elmmat(2,1) = elmmat(1,2);

        % update system matrix
        
        for i = 1 : 2
            if(lm(n,i)~=0)
                for j = 1 : 2
                    if(lm(n,j)~=0)
                        sysmat(lm(n,i),lm(n,j)) = sysmat(lm(n,i),lm(n,j)) + elmmat(i,j);
                    end
                end
            end
        end
        
        % elment vector
        
        elmvec(1) =-shp(1)*f(xq)*wq; 
	    elmvec(2) =-shp(2)*f(xq)*wq;
        if(n == numelm)
            elmvec(1) = elmvec(1) + h;
        end
        if(n == 1)
            elmvec(1) = elmvec(1) - elmmat(1,2)*g;
            elmvec(2) = elmvec(2) - elmmat(2,2)*g;
        end

        % update system vector

        for i = 1 : 2
            if(lm(n,i) ~= 0)
                sysvec(lm(n,i)) = sysvec(lm(n,i)) + elmvec(i);
            end
        end
    end
    
    % solve system matrix equations
    
    sysvar = sysmat\sysvec;
    
    % assign nodal values
    
    u = [g; sysvar];
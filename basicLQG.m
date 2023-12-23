function [L,cost] = basicLQG(A,B,Q,R,x0,oZeta)

% [L,COST] = basicLQG(A,B,Q,R,x0,oZeta)
%
%   A, B: Model matrices, x[k+1] = Ax[x]+Bu[k]
%   Q, R: Cost parameters
%   X0: Initial position
%   OZETA: Covariance matrix of the system noise
%
%   L: Time series of optimal feedback control gains
%   COST: Total expected cost

n = size(Q,3);
St = Q(:,:,end);
scost = 0;
sizeState = size(A,1);
sizeControl = size(B,2);
L = zeros(sizeControl,sizeState,n-1);

for i = n-1:-1:1
    
    L(:,:,i) = (R(:,:,i)+B'*St*B)\B'*St*A;
    Sttemp = St;
    St = Q(:,:,i)+A'*Sttemp*(A-B*L(:,:,i));
    scost = scost+trace(Sttemp*oZeta);
    
end
    
cost = x0'*St*x0+scost;    




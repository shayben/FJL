function [Y,e]=FJL(X,k)
%FJL - Fast Johnson–Lindenstraus transform attempt
%X = original matrix.
%k = target dimensions.
err=10^-2;
C=10; %what is this, i dont even?
X=X'; %points need to be column vectors. Matlab usually keep them as row vectors.
[d,n]=size(X);
H=HadamardForuier(d);
D=diag(ones(1,d)-2*(randi(2,1,d)-1));
R=round(C*k^2/err^2);
S=SparseGausMat(k,d,R);
A=S*H*D;
Y=A*X;

end


function H=HadamardForuier(d)
    H=zeros(d);
    len=nextpow2(d);
    for i=1:d
        %should figure out if zero or one based
        in=cellfun(@str2num, cellstr(sprintf('%c',dec2bin(i,len))'));
        for j=1:d
           jn=cellfun(@str2num, cellstr(sprintf('%c',dec2bin(j,len))'));
           H(i,j)=(-1)^dot(in,jn);
        end
    end
    sqrt(1/d)*H;
end

function S=SparseGausMat(k,d)
Cr=d/(k*R);
S=sparse(zeros(k,d));
for i=1:k
    J=randi(d,R);
    for j=1:R
        S(i,J(j))=S(i,J(j))+normrnd(0,1);
    end
end
S=S.*Cr;
end
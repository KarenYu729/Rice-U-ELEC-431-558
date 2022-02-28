% coefficients of H(z)
a1 = 0.7;
a2 = -0.63;
a3 = -0.162;
den_H = [1,-a1,-a2,-a3];
% roots of H(z)
root_H = roots(den_H);
% coefficient of X(z)
A1 = 0.9^4*exp((pi/2)*1j);
b1 = 0.9*exp((pi/4)*1j);
b2 = exp((2*pi/5)*1j);
b2_con = conj(b2);
% numerator of Y(z) i.e.X(z)
% num_X1 = conv([0,2,-2*b1,0,A1],[1,-(b2+b2_con),b2*b2_con,0,0]);
% num_X1 = num_X1(1:7);
% num_X = num_X1-conv([5,-5*real(b2),0,0],[1,-b1,0,0]);
num_X1 = conv([0,2,-2*b1,0,A1],[1,-2*real(b2),abs(b2),0,0])-conv([5,-5*real(b2),0,0,0],[1,-b1,0,0,0]);
num_X = num_X1(1:7);
% num_X -> coefficient of the polynomial
% the nominator of Y is written in polynomial style
% denominator of Y(z)
Y_de1 = [1,-b1];
Y_de2 = [1,-b2];
Y_de3 = [1,-b2_con];
denom_Y = conv(conv(Y_de1,Y_de2),conv(Y_de3,den_H));
root_Y = roots(denom_Y);
% now the denominator of Y is written in factor style
root_R_ex = kron(root_Y,ones(1,6)).';
% num is a6/d6
num = deconv(num_X(7),denom_Y(7));
Nres = num_X - num*denom_Y;
Nres = Nres(1:6);
% N(z)|ri
num_R = root_R_ex.^(kron([0 -1 -2 -3 -4 -5],ones(6,1)).');
num_R = Nres*num_R;
%D(z)/(1-ri*z^-1)|ri
denum_f = 1-root_R_ex./(root_R_ex.');
denum_f = denum_f+diag(ones(1,6));
denum_f = prod(denum_f,2);
A = (num_R.')./denum_f;
seqn = [0:1999];
seqx = -5*cos(2*pi*seqn/5);
seqx(2) = seqx(2) + 2;
seqx(5:end) = seqx(5:end) + (exp(-j*pi/2)*((0.9*exp(j*pi/4)).^seqn(5:end)));
seqy = A(1)*root_Y(1).^seqn+A(2)*root_Y(2).^seqn+A(3)*root_Y(3).^seqn+A(4)*root_Y(4).^seqn+A(5)*root_Y(5).^seqn+A(6)*root_Y(6).^seqn;
seqy(1) = seqy(1)+num;
newxSeq = [0 0 0 seqx];
newySeq = [0 0 0 seqy];
testright = newxSeq(4:end);
testright(4:end) = testright(4:end)+a1*newySeq(3:(end-4));
testright(4:end) = testright(4:end)+a2*newySeq(2:(end-5));
testright(4:end) = testright(4:end)+a3*newySeq(1:(end-6));
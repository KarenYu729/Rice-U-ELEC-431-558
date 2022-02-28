clc;
clf;
filt3 = load("HW3_filt.mat");
% filt3.hw3_filt length 45
H_z = [filt3.hw3_filt];
H_len = length(H_z);

k = [0:999];
len_k = length(k);
zk = exp(1j*2*pi*k/1000);
% root of H(z)
N_1_Root = roots(H_z);
figure(1);
plot(N_1_Root, '*r');title('roots of polynomial H(z)')
% % plug zk in H(z)
% % 1-[roots].*[zk].^[-1]
% ext_Roots = kron(N_1_Root,ones(1,len_k));
% ext_power = kron([-1], ones(H_len-1,1));
% ext_z = kron(zk, ones(H_len-1,1));
% plug_R = 1-(ext_Roots.*(ext_z.^ext_power));
prod_col = H_z*(kron(ones(45,1),zk).^kron([0:44]',ones(1,1000)));

res = sqrt(poly(N_1_Root)/H_z);

% prod_col = prod(plug_R);
figure(2);
plot(abs(prod_col));

% % design sub-system random choose roots
% index = randperm(44);
% % index1 = index(1:22);
% % index1 = [37 20 10 39 41 25 27 28 8 17 12 9 13 15 5 4 1 34 14 44 18 31];
% % index1 = [3 12 24 9 39 28 21 29 18 7 41 15 1 32 8 36 44 27 11 20 30 40];
% index1 = [6 26 30 2 25 22 36 15 32 9 5 39 27 44 3 33 20 34 10 13 28 42];
% Root_1P = N_1_Root(index1);
% % index2 = index(23:44);
% % index2 = [22 43 33 6 19 3 21 42 23 30 16 36 24 11 40 26 38 35 2 7 29 32];
% % index2 = [5 14 31 19 26 42 22 37 38 2 13 35 4 16 33 17 34 25 23 6 43 10];
% index2 = [18 23 41 38 12 14 8 1 24 21 11 19 29 43 16 17 4 37 40 31 35 7];
% Root_2P = N_1_Root(index2);
% newHsub1 = poly(Root_1P)/(res*1.2);
% newHsub2 = poly(Root_2P)/(res/1.2);
% % figure(3);
% % subplot(1,2,1);plot(Root_1P, '*r');
% % subplot(1,2,2);plot(Root_2P, '*r');
% figure(4);
% prod_P1 = newHsub1*(kron(ones(23,1),zk).^kron([0:22]',ones(1,1000)));
% subplot(1,2,1);plot(abs(prod_P1)); title('The first sub-system')
% prod_P2 = newHsub2*(kron(ones(23,1),zk).^kron([0:22]',ones(1,1000)));
% subplot(1,2,2);plot(abs(prod_P2)); title('The second sub-system')
% figure(5);
% plot(abs(prod_P1.*prod_P2));
% figure(6);
% plot(abs(prod_P1),'r');hold on;
% plot(abs(prod_P2),'b');hold on;


% design sub-system
Root_1P = N_1_Root(1:22);
Root_2P = N_1_Root(23:44);
newHsub1 = poly(Root_1P)/(res*50);
newHsub2 = poly(Root_2P)/(res/50);
figure(3);
subplot(1,2,1);plot(Root_1P, '*r');
subplot(1,2,2);plot(Root_2P, '*r');
figure(4);
prod_P1 = newHsub1*(kron(ones(23,1),zk).^kron([0:22]',ones(1,1000)));
subplot(1,2,1);plot(abs(prod_P1)); title('The first sub-system')
prod_P2 = newHsub2*(kron(ones(23,1),zk).^kron([0:22]',ones(1,1000)));
subplot(1,2,2);plot(abs(prod_P2)); title('The second sub-system')
figure(5);
plot(abs(prod_P1.*prod_P2));
figure(6);
plot(abs(prod_P1),'r');hold on;
plot(abs(prod_P2),'b');hold on;




% match sub-system
% index1 = index(1:22);
% index1 = [37 20 10 39 41 25 27 28 8 17 12 9 13 15 5 4 1 34 14 44 18 31];
% index1 = [3 12 24 9 39 28 21 29 18 7 41 15 1 32 8 36 44 27 11 20 30 40];
index1 = [6 26 30 2 25 22 36 15 32 9 5 39 27 44 3 33 20 34 10 13 28 42];
Root_1P = N_1_Root(index1);
% index2 = index(23:44);
% index2 = [22 43 33 6 19 3 21 42 23 30 16 36 24 11 40 26 38 35 2 7 29 32];
% index2 = [5 14 31 19 26 42 22 37 38 2 13 35 4 16 33 17 34 25 23 6 43 10];
index2 = [18 23 41 38 12 14 8 1 24 21 11 19 29 43 16 17 4 37 40 31 35 7];
Root_2P = N_1_Root(index2);
newHsub1 = poly(Root_1P)/(res*1.2);
newHsub2 = poly(Root_2P)/(res/1.2);
% figure(3);
% subplot(1,2,1);plot(Root_1P, '*r');
% subplot(1,2,2);plot(Root_2P, '*r');
figure(7);
prod_P1 = newHsub1*(kron(ones(23,1),zk).^kron([0:22]',ones(1,1000)));
subplot(1,2,1);plot(abs(prod_P1)); title('The first sub-system')
prod_P2 = newHsub2*(kron(ones(23,1),zk).^kron([0:22]',ones(1,1000)));
subplot(1,2,2);plot(abs(prod_P2)); title('The second sub-system')
figure(8);
plot(abs(prod_P1.*prod_P2));
figure(9);
plot(abs(prod_P1),'r');hold on;
plot(abs(prod_P2),'b');hold on;
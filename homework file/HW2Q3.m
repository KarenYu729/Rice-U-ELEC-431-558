% define polynomial
my_poly1 = [1 2 3 4 3 2 1];
my_poly2 = [1 1 1 1 1 1 1];
% length of th polynomial
poly1_len = length(my_poly1);
poly2_len = length(my_poly2);
%roots of the polynomial
poly1_root = roots(my_poly1);
poly2_root = roots(my_poly2);
% plot roots(zeros)
figure(1);
subplot(2,1,1);plot(poly1_root, '*r');title('roots of polynomial 1');
subplot(2,1,2);plot(poly2_root, '*r');title('roots of polynomial 2');
% creat finite grid
grange = 2;
delta = 1/100;
oneD_grid = [-grange:delta:grange];
gridN = length(oneD_grid);
% creat gridN × gridN matrix
AA_real = kron(ones(gridN,1),oneD_grid);
AA_imag = kron(oneD_grid'*1j,ones(1,gridN));
AA = AA_real+AA_imag;
% 
BB = kron(ones(1,poly1_len),AA);
powerBB = kron([0:-1:-(poly1_len)+1],ones(gridN));
coefBB1 = kron(my_poly1,ones(gridN));

CC1 = (BB.^powerBB).*coefBB1;
DD1 = reshape(CC1,gridN,gridN,[]);
EE1 = sum(DD1,3);

coefBB2 = kron(ones(gridN),my_poly2);

CC2 = (BB.^powerBB).*coefBB2;
DD2 = reshape(CC2,gridN,gridN,[]);
EE2 = sum(DD2,3);

figure(2);
mesh(oneD_grid,oneD_grid,max(min(abs(EE1),5),0));
figure(3);
mesh(oneD_grid,oneD_grid,max(min(abs(EE2),5),0));

















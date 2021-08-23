syms x11 x12 x13 x14 x15 x16 x17 ...
     x21 x22 x23 x24 x25 x26 x27 ...
     x31 x32 x33 x34 x35 x36 x37 ...
     W W1 W2 W3 W4 W5 W6 W7 ... 
     
A =[x11 x12 x13 x14 x15 x16 x17;
    x21 x22 x23 x24 x25 x26 x27;
    x31 x32 x33 x34 x35 x36 x37 ];
W = [W1 W2 W3 W4 W5 W6 W7];
diag_W = diag([W1 W2 W3 W4 W5 W6 W7]);

ans1 = A*diag_W*A';
ans2 = 0;
for i = 1:7  
    ans2 = ans2 + W(i)*A(:,i)*A(:,i)';
end    

ans1 = x11;
ans2 = x32;

if(ans1 == ans2)
    disp('yes')
end    
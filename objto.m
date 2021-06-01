function [f] = objto(x)

x1 = x(1);
x2 = x(2);

temp1 = x1^2+x2-11;
temp2 = x1+x2^2-7;

f = temp1^2+temp2^2;


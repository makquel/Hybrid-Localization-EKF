function A = AWrap(A)

i = A>90;
A(i) = A(i) - 90;
A(~i) = 90 - A(~i);

function y=where_nearest(A,b)
[Asort,index]=sort(abs(A(:)-b));
y=index(1);
clear; clc;
A = hilb(10);
B = ifft2(A);
fileID = fopen('ifft2.data','w');
fprintf(fileID,'%d\n',size(B,2));
fprintf(fileID,[repmat('(%.6f,%.6f)\t', 1, size(B,2)) '\n'], [real(reshape(B.', 1, [])); imag(reshape(B.', 1, []))]);
type('ifft2.data')
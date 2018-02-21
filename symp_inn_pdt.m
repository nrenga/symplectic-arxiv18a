function inn = symp_inn_pdt(X, Y)

inn = mod(sum(X.*fftshift(Y,2), 2), 2);

end
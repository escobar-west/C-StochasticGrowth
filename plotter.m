fileID = fopen('ODEsolution.txt', 'r');
formatSpec = '%f';
A = fscanf(fileID, formatSpec);

N = A(1);
trials = A(2);

figure
hold on

for i = 1:trials
    plot( A(3:3 + N), A( i*(N+1)+3: i*(N+1)+3 + N ))
end

hold off

fclose(fileID)

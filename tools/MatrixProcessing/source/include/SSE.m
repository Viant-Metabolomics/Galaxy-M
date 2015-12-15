function L=SSE(lambda, alpha, y) % Calculate sum of squared errors
N=size(y,2);
len=size(y,1);
for i=1:N
    z(:,i)=jglog(y(:,i),alpha,lambda);
end
s = 0;
%grand_mean = mean(z(:)); 
mean_spec=mean(z,2);
for i=1:N
    %row_mean = mean(z(:,i));
    for j=1:len
        %col_mean = mean(z(j,:));
        s = s + (z(j,i)-mean_spec(j,1))^2;
    end
end
% Plot to show turning point
%plot(lambda,s,'bo', 'MarkerFaceColor','k','MarkerSize', 10)
%hold on
L=s;
%disp([alpha lambda L])
end

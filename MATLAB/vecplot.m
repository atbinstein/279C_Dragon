% Inputs: The time history of each component of a vector
% Plots a subplot of each component

function vecplot(t,X)

if size(X,2) ~= length(t)
    error('Nope')
end

k = length(t);
n = size(X,1);

figure;
hold on
subplot(n,1,1)
plot(t,X(1,:))
for ii = 2:n
    subplot(n,1,ii);
    plot(t,X(ii,:));
end

end
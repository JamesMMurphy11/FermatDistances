function X = SampleNonUniform_Dumbbell(N,tau)
%Sample non-uniformly from a density using rejection sampling.
%Designed for the bells of a dumbell

X=2*rand(10*N,2)-[1,1];
X=X(X(:,1).^2+X(:,2).^2<1,:);

Score=4*(X(:,1)).^2+tau;
Comparison=rand(size(X ,1),1);
X=X(Score>Comparison,:);

X=X(randsample(size(X,1),N),:);


end

%Toss out (1-alpha) proportion of things in the middle


%{
Idx=find(abs(X(:,1))<1/2);
X(Idx(rand(length(Idx),1)>alpha),:)=[];
%}
function X = SampleGaussianMixture(N,tau)

%%

X=zeros(N,2);

for i=1:N
    NoSampleYet=1;
    val=rand(1);
    if val<tau

        X(i,:)=rand(1,2);

    elseif val<tau+(1-tau)/2

        while NoSampleYet
            X(i,:)=.16*randn(1,2)+[.2,.2];
            if ~(X(i,1)<0 || X(i,1)>1 || X(i,2)<0 || X(i,2)>1)
                NoSampleYet=0;
            end
        end

    elseif val>tau+(1-tau)/2

        while NoSampleYet
            X(i,:)=.16*randn(1,2)+[.7,.7];
            if ~(X(i,1)<0 || X(i,1)>1 || X(i,2)<0 || X(i,2)>1)
                NoSampleYet=0;
            end
        end
    end
end

figure; scatter(X(:,1),X(:,2))

end
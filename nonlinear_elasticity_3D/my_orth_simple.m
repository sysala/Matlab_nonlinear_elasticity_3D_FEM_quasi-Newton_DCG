function [ W_orth ] = my_orth_simple( V,W,tol,maxsize)
n=size(W,2);
    
W_orth=[];

 if ~exist('maxsize','var')
     % third parameter does not exist, so default it to something
      maxsize = 1000;
 end

  if ~exist('tol','var')
     % third parameter does not exist, so default it to something
      tol = 1e-6;
  end
 
if size(V,2)>=maxsize
    return
end

for i=1:n
    temp=W(:,i)/norm(W(:,i));
    if ~isempty(V)
        temp=temp-V*(temp'*V)';
        temp=temp-V*(temp'*V)';
%         temp=temp-V*(temp'*V)';
%         temp=temp-V*(temp'*V)';
%         temp=temp-V*(temp'*V)';
%         temp=temp-V*(temp'*V)';
    end
    if ~isempty(W_orth)
        for j=1:size(W_orth,2)
            temp=temp-dot(temp,W_orth(:,j))*W_orth(:,j);
        end
        for j=1:size(W_orth,2)
            temp=temp-dot(temp,W_orth(:,j))*W_orth(:,j);
        end
%         for j=1:size(W_orth,2)
%             temp=temp-dot(temp,W_orth(:,j))*W_orth(:,j);
%         end
%         for j=1:size(W_orth,2)
%             temp=temp-dot(temp,W_orth(:,j))*W_orth(:,j);
%         end
%         for j=1:size(W_orth,2)
%             temp=temp-dot(temp,W_orth(:,j))*W_orth(:,j);
%         end
%         for j=1:size(W_orth,2)
%             temp=temp-dot(temp,W_orth(:,j))*W_orth(:,j);
%         end
    end
    n_temp = norm(temp);
    if n_temp>tol
        W_orth=[W_orth temp/n_temp];
%     else
%         fprintf('!');
    end
end
end
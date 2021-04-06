factors = factor(length(H));
assert(max(factors) <= 5, 'This code cannot accept a matrix whose dimensions contain a prime factor greater than 5.');

sigma{1}=[0 1; 1 0];
sigma{2}=[0 -1i; 1i 0];
sigma{3}=[1 0 ; 0 -1];
sigma{4}=eye(2);

m{1}=[0 1 0;1 0 0;0 0 0];
m{2}=[0 -1i 0;1i 0 0;0 0 0];
m{3}=[1  0 0;0 -1 0;0 0 0];
m{4}=[0 0 1;0 0 0;1 0 0];
m{5}=[0  0 -1i;0  0  0; 1i 0  0];
m{6}=[0 0 0;0 0 1;0 1 0];
m{7}=[0  0  0;0  0 -1i;0 1i  0];
m{8}=[1 0 0;0 1 0;0 0 -2]/sqrt(3);
m{9}=eye(3);

five{1} = [0 2 0 0 0 ; 2 0 sqrt(6) 0 0 ; 0 sqrt(6) 0 sqrt(6) 0 ; 0 0 sqrt(6) 0 2 ; 0 0 0 2 0];
five{2} = [0 -2i 0 0 0 ; 2i 0 -sqrt(6)*1i 0 0 ; 0 sqrt(6)*1i 0 -sqrt(6)*1i 0 ;  0 0 sqrt(6)*1i 0 -2i ; 0 0 0 2i 0];
five{3} = [2 0 0 0 0 ; 0 1 0 0 0 ; 0 0 0 0 0 ; 0 0 0 -1 0 ; 0 0 0 0 -2];
five{4} = eye(5);

IND = ones([1 length(factors)]);

max = [];
for i = 1:length(factors)
    if factors(i) == 2
        max = [max, 4];
    elseif factors(i) == 3
        max = [max, 9];
    elseif factors(i) == 5
        max = [max, 4];
    end
end

while IND(1) <= max(1)
    temp = eye(1);
    for i = 1:length(IND)
        if factors(i) == 2
            temp = kron(temp,sigma{IND(i)});
        elseif factors(i) == 3
            temp = kron(temp,m{IND(i)});
        elseif factors(i) == 5
            temp = kron(temp,five{IND(i)});
        end
    end
    
    basis(:,sub2ind(max,IND(1),IND(2),IND(3)))=temp(:); % Figuring out how to do this for an arbitrary amount of particles is the last thing I need to do before converting to a function
    
    IND(length(IND)) = IND(length(IND)) + 1;
    for i = length(IND):-1:1
        if IND(i) > max(i) && i ~= 1;
            IND(i) = 1;
            if i ~= 1
                IND(i - 1) = IND(i - 1) + 1;
            end
        end
    end
end

c=sparse(round(basis\H(:)));

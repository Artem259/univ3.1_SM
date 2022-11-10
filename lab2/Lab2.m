%#ok<*NOPTS>
clear % видалення попередніх змінних

eps = 10^(-6);

X = imread("x1.bmp");
figure("Name", "X");
imshow(X);
X = double(X);
X = cat(1, X, ones(1,size(X,2)));

Y = imread("y8.bmp");
figure("Name", "Y");
imshow(Y);
Y = double(Y);

calculate(X, Inverse_with_Grevil(X), Y, "Grevil")
calculate(X, Inverse_with_Moor(X), Y, "Moor")

function calculate(X, X_inv, Y, method_name)
    A = Y*X_inv + zeros(size(Y,1), size(X,1))*Z(X_inv,X);
    figure("Name", method_name);
    imshow(uint8(A*X));
end

% ------------------------------------------

function res = Z(A, A_inv)
    res = eye(size(A_inv, 1)) - A_inv*A;
end

function res = R(A_inv)
    res = A_inv * (A_inv');
end

function res = Grevil(A, a, A_inv)
    r = R(A_inv);
    z = Z(A,A_inv);
    if ((a'*z*a) < eps)
        res = ((r*a)*a'*A_inv) / (1+a'*r*a);
        res = A_inv - res;
        res = cat(2, res, (r*a)/(1+a'*r*a));
        
    else
        res = ((z*a)*a'*A_inv) / (a'*z*a);
        res = A_inv - res;
        res = cat(2, res, (z*a)/(a'*z*a));
    end
end

function res = Inverse_with_Grevil(A)
    scalar = A(1,:) * A(1,:)';
    if (scalar == 0)
        A_inv = A(1,:);
    else
        A_inv = A(1,:) / scalar;
    end
    A_inv = A_inv';
    
    M = A(1,:);
    for i = 2:size(A,1)
        A_inv = Grevil(M, A(i,:)', A_inv);
        M = cat(1, M, A(i,:));
    end
    res = A_inv;
end

% ------------------------------------------

function res = Approximation(A, d)
    res = A' * (A*A' + (d^2)*eye(size(A,1)))^(-1);
end

function res = Inverse_with_Moor(A)
    d = 10.0;
    A_curr = Approximation(A, d);
    while 1
        d = d/2;
        A_next = Approximation(A, d);
        if (norm(A_next-A_curr) < eps)
            break
        end
        A_curr = A_next;
    end
    res = A_next;
end

% ------------------------------------------

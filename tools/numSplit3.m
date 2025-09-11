function number = numSplit3(k)

switch k
    case 0
        number = [0; 0; 0];
    case 1
        number = [1 0 0;
                  0 1 0;
                  0 0 1];
    case 2
        number = [2 1 1 0 0 0;
                  0 1 0 2 1 0;
                  0 0 1 0 1 2];
    case 3
        number = [3 2 2 1 1 1 0 0 0 0;
                  0 1 0 2 1 0 3 2 1 0;
                  0 0 1 0 1 2 0 1 2 3];
        
    case 4
        number = [4 3 3 2 2 2 1 1 1 1 0 0 0 0 0;
                  0 1 0 2 1 0 3 2 1 0 4 3 2 1 0;
                  0 0 1 0 1 2 0 1 2 3 0 1 2 3 4];
    case 5
        number = [5 4 4 3 3 3 2 2 2 2 1 1 1 1 1 0 0 0 0 0 0;
                  0 1 0 2 1 0 3 2 1 0 4 3 2 1 0 5 4 3 2 1 0;
                  0 0 1 0 1 2 0 1 2 3 0 1 2 3 4 0 1 2 3 4 5];
    otherwise
        number = zeros(3, (k+1)*(k+2)/2);
        cnt = 0;
        for i = 1 : k+1
            number(1, cnt + (1:i)) = k-i+1;
            number(2, cnt + (1:i)) = i-1:-1:0;
            number(3, cnt + (1:i)) = 0:i-1;
            cnt = cnt + i;
        end
end

end
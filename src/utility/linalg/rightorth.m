function [R, Q] = rightorth(A, alg)




% TODO have a look at https://github.com/iwoodsawyer/factor

switch alg
    case 'svd'
        [U, S, Q] = svd(A, 0);
        R = U * S;
        Q = Q';
        return;
        
    case 'polar'
        newalg = 'polar';
    case 'lq'
        newalg = 'qr';
    case 'lqpos'
        newalg = 'qrpos';
    case 'rq'
        newalg = 'ql';
    case 'rqpos'
        newalg = 'qlpos';
end

[Q, R] = leftorth(A.', newalg);
Q = Q.';
R = R.';

end

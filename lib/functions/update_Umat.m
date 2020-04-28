function U = update_Umat(Y, Fs, modes, rank)
% U = update_Umat(Y, Fs, modes, rank)
% Updates U matrix by pseudoinverse
% 
[U, S, V] = svds(t2m(Y, modes), rank);
s = 1./diag(S);
U = t2m(Fs, modes)*(V*diag(s)*U');
end

% function U = update_Umat(Y, Fs, modes, rank)
% % U = update_Umat(Y, Fs, modes, rank)
% % Updates U matrix by pseudoinverse
% % 
% tau = 0.3;
% m = length(modes);
% sz = size(Y);
% [U, S, V] = svds(t2m(Y, modes), rank);
% s = 1./diag(S);
% U = t2m(Fs, modes)*(V*diag(s)*U');
% U = reshape(U, [sz(modes), size(U, 2)]);
% Ui = U2Ui_tau(U, tau);
% temp = U2Ui_tau(reshape(Ui{m+1}',[sz(modes), size(Ui{m+1},1)]), tau);
% if length(temp)~= m+1
%     temp_nrm = norm(temp{m},'fro');
%     temp{m} = temp{m}/temp_nrm;
%     temp{m+1} = temp_nrm;
% end
% Ui([m+2:2*m+1,m+1]) = temp;
% U = Lunfold(mergeFactors(Ui(1:m)),m+2)*Ui{m+1}'*Lunfold(mergeFactors(Ui(m+2:end)),m+2)';
% end

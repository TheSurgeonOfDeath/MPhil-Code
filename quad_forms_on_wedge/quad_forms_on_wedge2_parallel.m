% H = R^n, we want to determine possible signatures of symmetric bilinear
% forms on Wedge^2 H, which can be identified with q in Wedge^4 H^*.

% parpool; 


n = 7;

% fid = fopen("unique_sgns_" + n + "_parallel.csv", "w");
% % x = [1,2,3];
% % % fprintf(fid, "1,2,3\n")
% % fprintf(fid, "%d,%d,%d\n", x);
% fclose(fid);


wedge_basis_idx = nchoosek(uint8(1:n),2);
k = length(wedge_basis_idx);
quad_forms_basis_idx = nchoosek(uint8(1:n),4);
m = size(quad_forms_basis_idx,1);

% quad_forms = @(i, b1, b2) q_form(quad_forms_basis_idx(i,:), b1, b2)
% quad_forms(1,b1,b2)

quad_forms = cell(1,m);
mats = cell(1,m);
for i = 1:m
    quad_forms{i} = @(b1,b2) q_form(quad_forms_basis_idx(i,:),b1,b2);
    mats{i} = symm_matrix(wedge_basis_idx, quad_forms{i});
end

% comb_mats = nchoosek(mats,2)
% comb_mats(1,:)

% allmats = sparse(k,k,2^m - 1)
% allmats = cell(1,2^m - 1);
% idxs = zeros(2^m - 1,2);
% allmats = {}

filename = "unique_sgns_" + n + "_parallel_with_forms.csv";
% filename = "unique_sgns_" + n + "_parallel_with_forms_8+.csv";
fid = fopen(filename, "w");

current_idx = 0;
unique_sgns = zeros(1,3);
for i = 1:n % only sum up to n terms because after that the indices start repeating.
% for i = n+1:m
    p = nchoosek(m,i);
    comb_mats = nchoosek(1:m,i);
    for j = 1:p
        sum_mat = mats{comb_mats(j,1)};
        for r = 2:i
            sum_mat = sum_mat + mats{comb_mats(j,r)};
        end
        current_mat = sparse(sum_mat);
        current_sgn = signature_matrix(current_mat);
        if ~ismember(current_sgn, unique_sgns, "rows")
            current_idx = current_idx + 1;
            unique_sgns(current_idx,:) = current_sgn
            fprintf(fid, "%d,%d,%d,", current_sgn);

            forms = quad_forms_basis_idx(comb_mats(j,:),:);
            str = forms_idx_str(forms);
            fprintf(fid,"%s\n", str);
            % comb_mats{j,:}
            % sum_mat
            [i j]
        end
    end
end
fclose(fid);

% unique_sgns = full(unique_sgns(1:current_idx,:));

% filename = "unique_sgns_" + n + "_parallel.csv";
% writematrix(unique_sgns, filename);



% A = [1 2 3; 4 5 6; 7 8 9];
% x = [4 5 6];
% 
% if ismember(x,A,"rows")
%     "yes"
% end

% signatures = zeros(length(allmats),3);
% for i = 1:length(allmats)
%     signatures(i,:) = signature_matrix(allmats{i});
% end
% signatures
% unique_sgns = unique(signatures,"rows")

% filename = "unique_sgns_" + n + ".mat";
% save(filename, "unique_sgns");

% signature_quad_form(wedge_basis_idx, quad_forms{1})

%%
% c = find(signatures(:,3) == 0)
% crit = 17; % (5 choose 1) + (5 choose 1) = 15 so in (5 choose 3).
% nondegen_mats = nchoosek(mats,3)
% 
% nondegen_mats{2,:}

%%
% q1 = @(b1,b2) q_form([1 2 3 4],b1,b2);
% q2 = @(b1,b2) q_form([1 3 4 5],b1,b2);
% q3 = @(b1,b2) q_form([1 2 4 5],b1,b2);
% mat1 = symm_matrix(wedge_basis_idx, q1)
% mat2 = symm_matrix(wedge_basis_idx, q2)
% mat3 = symm_matrix(wedge_basis_idx, q3)
% eig(mat1)
% sgn1 = signature_matrix(mat1)
% eig(mat1 + mat2)
% sgn2 = signature_matrix(mat1 + mat2)
% eig(mat1 + mat2 + mat3)
% sgn2 = signature_matrix(mat1 + mat2 + mat3)

% mats = cell(m);
% for i = 1:m
%     mats{i} = symm_matrix(wedge_basis_idx, quad_forms{i});
% end


% signatures = zeros(m,3);
% for i = 1:m
%     signatures(i,:) = signature_quad_form(wedge_basis_idx, quad_forms{i});
% end
% signatures;
% 
% signature_quad_form(wedge_basis_idx, quad_forms{1})

%%
q1 = @(b1,b2) q_form([1 2 3 4],b1,b2);
mat1 = symm_matrix(wedge_basis_idx, q1)

q2 = @(b1,b2) q_form([1 2 3 5],b1,b2);
mat2 = symm_matrix(wedge_basis_idx, q2)

A = mat1 + mat2

eig(A)

%%

% quad_forms_basis_idx_cell = cell(1,4);
% for row = 1:m
%     quad_forms_basis_idx_cell{row} = quad_forms_basis_idx(row, :);
% end
% quad_forms_combs = nchoosek(quad_forms_basis_idx_cell,i)

comb_mats = nchoosek(1:m,2)
mats(comb_mats(6,:))

quad_forms_basis_idx(comb_mats(20,:),:)

forms = quad_forms_basis_idx(comb_mats(20,:),:)
forms(1,:)
string(mat2str(forms(1,:)))

str = "" + string(mat2str(forms(1,:)));
for i = 2:size(forms,1)
    str = str + "," + string(mat2str(forms(2,:)));
end
str

% fid = fopen("test.csv", "w")
% sgn = [1 2 3]
% fprintf(fid, "%d,%d,%d,", sgn)
% fprintf(fid, "%s", quad_forms_basis_idx(comb_mats(20,:),:))
% fclose(fid);



% comb_mats(j,:)
% 
% index = cellfun(@(x) x==comb_mats{j,1}, mats, 'UniformOutput', 1);
% 
% find(mats{:} == comb_mats{j,1})


function str = forms_idx_str(forms)
    str = "" + string(mat2str(forms(1,:)));
    for i = 2:size(forms,1)
        str = str + "," + string(mat2str(forms(2,:)));
    end
end



function mat = symm_matrix(wedge_basis_idx, quad_form)
    k = length(wedge_basis_idx);
    indices = zeros(3,2);
    vals = zeros(3,1);
    % mat = sparse(k,k);
    current_idx = 0;
    for i = 1:k
        for j = 1:i
            val = quad_form(wedge_basis_idx(i,:), wedge_basis_idx(j,:));
            if val
                current_idx = current_idx + 1;
                indices(current_idx,:) = [i,j];
                vals(current_idx) = val;
            end
        end
    end
    upper = sparse(indices(:,1), indices(:,2),vals, k,k);
    lower = sparse(indices(:,2), indices(:,1),vals,k,k);
    mat = upper + lower;
end

function sgn = signature_matrix(symm_matrix)
    eig_mat = eig(full(symm_matrix));
    p = length(eig_mat(eig_mat > 1e-8));
    q = length(eig_mat(eig_mat < -1e-8));
    r = size(symm_matrix,1) - (p + q);
    sgn = [p,q,r];
end

function sgn = signature_quad_form(wedge_basis_idx, quad_form)
    k = length(wedge_basis_idx);
    symm_matrix = zeros(k,k);
    for i = 1:k
        for j = 1:k
            symm_matrix(i,j) = quad_form(wedge_basis_idx(i,:), wedge_basis_idx(j,:));
        end
    end
    eig_mat = eig(symm_matrix);
    p = length(eig_mat(eig_mat > 0));
    q = length(eig_mat(eig_mat < 0));
    r = k - (p + q);

    sgn = [p,q,r];
end

function k = q_form(q,b1,b2)
    if sum(ismember(q,[b1,b2])) ~= 4
        k = 0;
        return
    end
    x = logical(alt_sum(ismember(q,b1)));
    k = (-1)^x;
end


function y = alt_sum(x)
    y = sum(x(1:2:end))- sum(x(2:2:end));
end


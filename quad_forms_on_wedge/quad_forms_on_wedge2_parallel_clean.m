% H = R^n, we want to determine possible signatures of symmetric bilinear
% forms on Wedge^2 H, which can be identified with q in Wedge^4 H^*.

% parpool; 


n = 7;


% basis for wedge^2 H
wedge_basis_idx = nchoosek(uint8(1:n),2);
k = length(wedge_basis_idx); % n choose 2

% basis for wedge^4 H
quad_forms_basis_idx = nchoosek(uint8(1:n),4);
m = size(quad_forms_basis_idx,1); % n choose 4


% all quadrilinear forms % Vectorize matrices
quad_forms = cell(1,m);
mats = cell(1,m);
for i = 1:m
    quad_forms{i} = @(b1,b2) q_form(quad_forms_basis_idx(i,:),b1,b2);
    mats{i} = symm_matrix(wedge_basis_idx, quad_forms{i});
end




% file to save data
filename = "unique_sgns_" + n + "_parallel_with_forms.csv";
fid = fopen(filename, "w");

% generate signatures by generating all linear combinations of basic quad_forms
current_idx = 0;
seen_sgns = containers.Map('KeyType', 'char', 'ValueType', 'logical');
for i = 1:n % only sum up to n terms because after that the indices start repeating.
    p = nchoosek(m,i);
    comb_mats = nchoosek(1:m,i);
    for j = 1:p
        sum_mat = mats{comb_mats(j,1)};
        for r = 2:i
            sum_mat = sum_mat + mats{comb_mats(j,r)};
        end
        current_mat = sparse(sum_mat);
        current_sgn = signature_matrix(current_mat);
        sgn_key = sprintf('%d_%d_%d', current_sgn);

        % if signature is not unique, skip
        if isKey(seen_sgns, sgn_key)
            continue;
        end
        seen_sgns(sgn_key) = true;

        % write signature to file
        current_idx = current_idx + 1;
        % unique_sgns(current_idx,:) = current_sgn
        % current_sgn
        fprintf(fid, "%d,%d,%d", current_sgn);
        
        % write forms to file
        forms = quad_forms_basis_idx(comb_mats(j,:),:);
        str = forms_idx_str(forms);
        fprintf(fid,",%s\n", str);
        [i j]
    end
end
fclose(fid);




function str = forms_idx_str(forms)
    str = string(mat2str(forms(1,:)));
    for i = 2:size(forms,1)
        str = str + "," + string(mat2str(forms(i,:)));
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


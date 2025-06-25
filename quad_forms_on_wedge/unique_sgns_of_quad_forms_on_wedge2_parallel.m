% H = R^n, we want to determine possible signatures of symmetric bilinear
% forms on Wedge^2 H, which can be identified with q in Wedge^4 H^*.
n = 5;

wedge_basis_idx = nchoosek(1:n,2);
k = length(wedge_basis_idx);
quad_forms_basis_idx = nchoosek(1:n,4);
m = size(quad_forms_basis_idx,1);

% quad_forms = @(i, b1, b2) q_form(quad_forms_basis_idx(i,:), b1, b2)
% quad_forms(1,b1,b2)

quad_forms = cell(1,m);
mats = cell(1,m);
for i = 1:m
    quad_forms{i} = @(b1,b2) q_form(quad_forms_basis_idx(i,:),b1,b2);
    mats{i} = symm_matrix(wedge_basis_idx, quad_forms{i});
end

comb_mats = nchoosek(mats,2)
comb_mats(1,:)

%%

allmats = cell(1,2^m - 1);
current_idx = 0;
for i = 1:m
    p = nchoosek(m,i);
    comb_mats = nchoosek(mats,i);
    for j = 1:p
        current_idx = current_idx + 1;
        c_mats = comb_mats(j,:);
        sum_mat = comb_mats{j,1};
        for r = 2:i
            sum_mat = sum_mat + comb_mats{j,r};
        end
        allmats{current_idx} = sum_mat;
    end
end

%%
signatures = zeros(length(allmats),3);
for i = 1:length(allmats)
    signatures(i,:) = signature_matrix(allmats{i});
end
signatures

unique_sgns = unique(signatures,"rows");

% filename = "unique_sgns_" + n + ".mat";
% save(filename, "unique_sgns");
filename = "unique_sgns.csv";
writematrix(unique_sgns, filename);

% signature_quad_form(wedge_basis_idx, quad_forms{1})

%%
q1 = @(b1,b2) q_form([1 2 3 4],b1,b2);
q2 = @(b1,b2) q_form([1 3 4 5],b1,b2);
q3 = @(b1,b2) q_form([1 2 4 5],b1,b2);
mat1 = symm_matrix(wedge_basis_idx, q1)
mat2 = symm_matrix(wedge_basis_idx, q2)
mat3 = symm_matrix(wedge_basis_idx, q3)
eig(mat1)
sgn1 = signature_matrix(mat1)
eig(mat1 + mat2)
sgn2 = signature_matrix(mat1 + mat2)
eig(mat1 + mat2 + mat3)
sgn2 = signature_matrix(mat1 + mat2 + mat3)

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

function mat = symm_matrix(wedge_basis_idx, quad_form)
    k = length(wedge_basis_idx);
    mat = zeros(k,k);
    for i = 1:k
        for j = 1:k
            mat(i,j) = quad_form(wedge_basis_idx(i,:), wedge_basis_idx(j,:));
        end
    end
end

function sgn = signature_matrix(symm_matrix)
    eig_mat = eig(symm_matrix);
    p = length(eig_mat(eig_mat > 0));
    q = length(eig_mat(eig_mat < 0));
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


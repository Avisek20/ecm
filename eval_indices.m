function [t_ari, t_nmi, t_wli] = eval_indices(data, centers, u, labels)
    [~, mem] = max(u,[],2);
    t_ari = adjrand(labels,mem);
    t_nmi = nmi(labels,mem);
    t_wli = WLI(data, centers, u);
end
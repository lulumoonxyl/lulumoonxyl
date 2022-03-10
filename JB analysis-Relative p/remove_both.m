function [dat_remove_120s_edge,fields_remove_120s_edge]=remove_both(dat_JB,fields_JB,edge,ns,both)
if ~isempty(both)&~isempty(edge)&~isempty(ns)
    [dat_edge_removed,fields_edge_removed]=eliminate_edge(dat_JB,fields_JB,edge);
    [dat_remove_120s_edge,fields_remove_120s_edge]=remove_the_first_ns(dat_edge_removed,fields_edge_removed,ns);
else
    dat_remove_120s_edge={};fields_remove_120s_edge={};
end
end 
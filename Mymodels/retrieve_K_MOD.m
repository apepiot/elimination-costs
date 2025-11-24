function [k,mod] = retrieve_K_MOD(k_mod)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if iscell(k_mod)
    kmod = k_mod{:};
end

if ischar(k_mod)
    kmod = k_mod;
end
idx = strfind(kmod,'_');
k   = kmod(1:(idx-1));
mod = kmod(idx+1:end);
end


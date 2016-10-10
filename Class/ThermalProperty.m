classdef Property < handle
    properties
        temperature;
        cp_def;rho_def;mu_def;k_def;
    end
    
    methods
        function y=cp(varargin)
            y=[];
        end
        
        function y=rho(varargin)
            y=[];
        end
        
        function y=mu(varargin)
            y=[];
        end
        
        function y=k(varargin)
            y=[];
        end
        
        function y=setdefault(obj,T)
            obj.temperature=T;
            obj.cp_def=obj.cp(T);
            obj.rho_def=obj.rho(T);
            obj.mu_def=obj.mu(T);
            obj.k_def=obj.k(T);
            y=obj;
        end
        
        function [cp_T,rho_T,mu_T,k_T]=getall(obj,T)
            cp_T=obj.cp(T);
            rho_T=obj.rho(T);
            mu_T=obj.mu(T);
            k_T=obj.k(T);
        end
    end
end
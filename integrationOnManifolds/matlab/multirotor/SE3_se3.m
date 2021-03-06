classdef SE3_se3_ < handle
    %This class represents the semi direct product of SE(3) and little
    % se(3)
    
    properties
        g;     % Element of SE(3)
        v;     % Element of se(3)
    end
    
    methods(Static)
        
        % Performs group adjoint action on the 
        % lie algebra
        function v_new = Ad_SE3_se3(g,v)
            
            R2 = g(1:3,1:3);
            t2 = g(1:3,4);
            g_inv = [R2' -R2'*t2;zeros(1,3), 1];
            
            
            v_new = g*v*g_inv;
%             v_new = g_inv*v*g;
            
        end
        
        % Performs the group co adjoint action on the 
        % lie algebra
        function v_new = CO_Ad_SE3_se3(g,v)
            
            R2 = g(1:3,1:3);
            t2 = g(1:3,4);
%             M = [R2 zeros(3,3); -SE3_se3.SO3_wedge(t2)*R2, R2];
        M = [R2' zeros(3,3); -R2'*SE3_se3.SO3_wedge(t2), R2'];
            V = M*SE3_se3.SE3_vee(v);
            
%             g_inv = [R2' -R2'*t2;zeros(1,3), 1];
            
            
            v_new = SE3_se3.SE3_wedge(V);
%             v_new = g_inv*v*g;
            
        end
        
        % Returns the identity element
        function obj_new = I()
            g = eye(4);
            v = zeros(4,4);
            obj_new = SE3_se3(g,v);
        end
        
        % SSM
        function w_wedge = SO3_wedge(w)
            w1 = w(1);
            w2 = w(2);
            w3 = w(3);

            w_wedge = [0 -w3 w2;...
                w3  0  -w1;...
                -w2 w1 0];
 
        end
        
        function w_vee = SO3_vee(w_wedge)
            x = w_wedge(3,2);
            y = w_wedge(1,3);
            z = w_wedge(2,1);

            w_vee = [x;y;z];
        end
        
        % SE3 Vee
        function W_vee = SE3_vee(W_wedge)
            omega = W_wedge(1:3,1:3);
            omega_vee = SE3_se3.SO3_vee(omega);
            rho = W_wedge(1:3,4);
            
            W_vee = [rho;omega_vee];
        end
        
        % SE3 Wedge
        function W_wedge = SE3_wedge(W_vee)
           W_wedge=[SE3_se3.SO3_wedge(W_vee(4:6)),W_vee(1:3);zeros(1,4)]; 
        end
        
        % se3 Matrix Adjoint
        % v is the 4x4 matrix
        % m is the 6x6 matrix representation
        function m = se3_m_ad(v)
            omega = v(1:3,1:3);
            rho = v(1:3,4);
            
            m = [omega, SE3_se3.SO3_wedge(rho);zeros(3,3), omega];
        end
        
        
        % Note that v is the norm 4x4 matrix
        % but s is the 6x6 representation 
        function dexp_inv = se3_dexp_inv(v)
%             B = [1, -0.5, 1/6, 0, -1/30, 0, 1/42, 7, -1/30, 0, 5/66, 0, -691/2730,0,7/6]; % Bernouli numbers
%             s = 0;
%             k = length(B);
%             ad = SE3_se3.se3_m_ad(v);
% 
%             for ii = 0:k-1
%               s = s + (-1)^ii*B(ii+1)/factorial(ii)*ad^ii;  
% 
%             end

            w = v(1:3,1:3);
            th = norm(w);
            dexp_inv_w = eye(3)+0.5*w-(th*cot(th/2)-2)/(2*th^2)*w^2;
            
            p = v(1:3,4);
            px = SE3_se3.SO3_wedge(p);
            a_th = (cos(th)-1)/th^2;
            b_th=(th-sin(th))/th^3;
            c_th = -sin(th)/th^3+2*((1-cos(th))/th^4);
            d_th = -2/th^4+3/th^5*sin(th)-1/th^4*cos(th);
            q = (SE3_se3.SO3_vee(w)'*p)*(c_th*w+d_th*w^2);
            B = (a_th*px+b_th*(w*px+px*w)+q);
            
            dexp_inv = [dexp_inv_w -dexp_inv_w*B*dexp_inv_w;...
                        zeros(3,3), dexp_inv_w];
        end
        
    end
    
    methods
        function obj = SE3_se3(g,v)
            
            obj.g = g;
            obj.v = v;             
            
        end
        
        % Performs left group multiplication
        % Using the left semi-direct product
        function obj3 = mtimes(obj1,obj2)
            
            g3 = obj1.g*obj2.g;
            R2 = obj2.g(1:3,1:3);
            t2 = obj2.g(1:3,4);
            g2_inv = [R2' -R2'*t2;zeros(1,3), 1];
            v3 = obj2.v + obj1.Ad_SE3_se3(g2_inv,obj1.v);
            
%             v3 = obj2.v + obj1.CO_Ad_SE3_se3(obj2.g,obj1.v);
            obj3 = SE3_se3(g3,v3);
            
        end
        
        % Computes the invese element
        function obj_inv = inv(obj)
           
            g2 = inv(obj.g);
            v2 = -obj.Ad_SE3_se3(obj.g,obj.v);
            v2 = round(v2,14);
            obj_inv = SE3_se3(g2,v2);
            
        end
        
        function [v,u] = log(obj)
            
            v = logm(obj.g);
            dv = SE3_se3.se3_dexp_inv(v);
            u = SE3_se3.SE3_wedge(dv*SE3_se3.SE3_vee(obj.v));
            
%             w = v(1:3,1:3);
%             th = norm(w);
%             dexp_inv_w = eye(3)+0.5*w-(th*cot(th/2)-2)/(2*th^2)*w^2;
%             
%             p = v(1:3,4);
%             px = SE3_se3.SO3_wedge(p);
%             a_th = (cos(th)-1)/th^2;
%             b_th=(th-sin(th))/th^3;
%             c_th = -sin(th)/th^3+2*((1-cos(th))/th^4);
%             d_th = -2/th^4+3/th^5*sin(th)-1/th^4*cos(th);
%             q = (SE3_se3.SO3_vee(w)'*p)*(c_th*w+d_th*w^2);
%             B = (a_th*px+b_th*(w*px+px*w)+q);
%             
%             dexp_inv = [dexp_inv_w -dexp_inv_w*B*dexp_inv_w;...
%                         zeros(3,3), dexp_inv_w];
                    
    
                    
%             u2 = SE3_se3.SE3_wedge(dexp_inv*SE3_se3.SE3_vee(obj.v))
           
            
        end
        

        
        
    end
    

end


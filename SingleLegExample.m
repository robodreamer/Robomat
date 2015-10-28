%% 06-28-14 Single Leg Floating-base Example

% written by Andy Park at Purdue

classdef SingleLegYan < ConstructRobotSpatial
    
    properties (SetAccess=protected)
        name = 'Single Leg - Yan'
        NJ = 6;
        NB = 6 + 1; % number of bodies: 11 + 1
        N_joint = 6 + 6; % number of joints: 11 + 6(FB)
        N_ext = 7 + 1; % total number of bodies + dummy bodies (end-effector frames)
        ag_0 = [0 0 0 0 0 -9.81]';
        display_mass = 1;
        display_inertia = 1;
    end
    
    properties (Constant = true)
        J_LHY = 1;
        J_LHR = 2;
        J_LHP = 3;
        J_LKP = 4;
        J_LAP = 5;
        J_LAR = 6;
    end
    
    properties (Access=public) 

    end
    
    % useful functions
    methods
           
        %%% construct a simple humanoid model
        function model = SingleLegYan()
            
            R_base = eye(3);
            
            L = 10;
            m_link = 4;
            scale = L;
            
            q = zeros(1,model.N_joint);
%             qd = zeros(1,model.N_joint);
            
            %=== left leg
            % base body (floating)
            KBody = 1;
            model.bodyname{KBody} = 'base';
            model.parent{KBody} = 0; % parent body (0: none)
            %             model.phi{KBody}    = eye(6); % motion subspace (joint motion that moves this body)
            model.i_p_com{KBody}  = [0 0 -L/8]';
            model.m{KBody}      = m_link;
            model.I_com{KBody}  = model.Inertia([3 3 1 0 0 0])*scale;
            model.phi{KBody}    = [zeros(3) eye(3); eye(3) zeros(3)];
            model.pi_p_i{KBody}   = q(1:3)'; % position of the joint w.r.t inertial frame
            model.pi_R_i{KBody}   = R_base*0 + eye(3); % rotation of the joint
            model.joints{KBody}  = [1:KBody + 5]; % joints (moving this body)
            model.q{KBody}      = q(1:KBody + 5)'; % initial joint angles
            
            KBody = KBody + 1;
            model.bodyname{KBody} = 'bodyLHY';
            model.parent{KBody} = model.FindBody('base');
            model.i_p_com{KBody}  = [0 0 0]';
            model.m{KBody}      = m_link/5;
            model.I_com{KBody}  = model.Inertia([1 1 1 0 0 0])*scale;
            model.phi{KBody}    = [0 0 1 0 0 0]';
            model.pi_p_i{KBody}   = [0 0 -L/4]';
            model.pi_R_i{KBody}   = eye(3);
            model.joints{KBody} = [KBody + 5];
            model.q{KBody}      = q(KBody + 5)';
            
            KBody = KBody + 1;
            model.bodyname{KBody} = 'bodyLHR';
            model.parent{KBody} = model.FindBody('bodyLHY');
            model.i_p_com{KBody}  = [0 0 0]';
            model.m{KBody}      = m_link/5;
            model.I_com{KBody}  = model.Inertia([1 1 1 0 0 0]);
            model.phi{KBody}    = [1 0 0 0 0 0]';
            model.pi_p_i{KBody}   = [0 0 -L/20]';
            model.pi_R_i{KBody}   = eye(3);
            model.joints{KBody} = [KBody + 5];
            model.q{KBody}      = q(KBody + 5)';
            
            KBody = KBody + 1;
            model.bodyname{KBody} = 'bodyLHP';
            model.parent{KBody} = model.FindBody('bodyLHR');
            model.i_p_com{KBody}  = [0 0 -L/2]';
            model.m{KBody}      = m_link;
            model.I_com{KBody}  = model.Inertia([20 20 1 0 0 0])*scale;
            model.phi{KBody}    = [0 1 0 0 0 0]';
            model.pi_p_i{KBody}   = [0 0 -L/20]';
            model.pi_R_i{KBody}   = eye(3);
            model.joints{KBody} = [KBody + 5];
            model.q{KBody}      = q(KBody + 5)';
            
            KBody = KBody + 1;
            model.bodyname{KBody} = 'bodyLKP';
            model.parent{KBody} = model.FindBody('bodyLHP');
            model.i_p_com{KBody}  = [0 0 -L/2]';
            model.m{KBody}      = m_link;
            model.I_com{KBody}  = model.Inertia([20 20 1 0 0 0])*scale;
            model.phi{KBody}    = [0 1 0 0 0 0]';
            model.pi_p_i{KBody}   = [0 0 -L]';
            model.pi_R_i{KBody}   = eye(3);
            model.joints{KBody} = [KBody + 5];
            model.q{KBody}      = q(KBody + 5)';
            
            KBody = KBody + 1;
            model.bodyname{KBody} = 'bodyLAP';
            model.parent{KBody} = model.FindBody('bodyLKP');
            model.i_p_com{KBody}  = [0 0 -L/4]';
            model.m{KBody}      = m_link;
            model.I_com{KBody}  = model.Inertia([5 5 5 0 0 0])*scale;
            model.phi{KBody}    = [0 1 0 0 0 0]';
            model.pi_p_i{KBody}   = [0 0 -L]';
            model.pi_R_i{KBody}   = eye(3);
            model.joints{KBody} = [KBody + 5];
            model.q{KBody}      = q(KBody + 5)';
            
            KBody = KBody + 1;
            model.bodyname{KBody} = 'bodyLAR';
            model.parent{KBody} = model.FindBody('bodyLAP');
            model.i_p_com{KBody}  = [0 0 0]';
            model.m{KBody}      = m_link/5;
            model.I_com{KBody}  = model.Inertia([1 1 1 0 0 0]);
            model.phi{KBody}    = [1 0 0 0 0 0]';
            model.pi_p_i{KBody}   = [0 0 -L/20]';
            model.pi_R_i{KBody}   = eye(3);
            model.joints{KBody} = [KBody + 5];
            model.q{KBody}      = q(KBody + 5)';
            
            % LF
            KBody = KBody + 1;
            model.bodyname{KBody} = 'LF';
            model.parent{KBody} = model.FindBody('bodyLAR');
            model.phi{KBody}    = [];
            model.i_p_com{KBody}  = [0 0 0]';
            model.m{KBody}      = 0;
            model.I_com{KBody}  = diag([0,0,0]);
            model.pi_p_i{KBody}   = [0, 0, -L/2]';
            model.pi_R_i{KBody}   = eye(3);
            model.joints{KBody} = [];
            model.q{KBody}       = [];
            
        end
        
        % add apperance to the featherstone (v2) model
        function model = AddAppearance(model, floor_z)
            m = model.featherstone;
            % floor
            % tiles profile: ( X, Y, Z, tilesize, colour )
            if ~exist('floor_z','var')
                floor_z = -28; % set floor height by default
            end
                
            m.appearance.base = { 'tiles', [-50 50; -50 50; floor_z floor_z], 5 };
            
            % add cylindrical bodies
            scale = 1;
            for I=model.N_ext:-1:1
                
                % find child links
                child_links = find([model.parent{:}] == I);
                if(isempty(child_links)) % end-bodies
                    %                         m.appearance.body{I-5} = { 'cyl', [0 0 0;0 0 -scale*10], scale*1};
                else
                    for k=1:length(child_links)
                        i_child = child_links(k);
                        % get the relative position parameter between bodies
                        % the child body (pi_p_i) parameter will be the
                        % parameter for the ith body here
                        p_xyz = model.pi_p_i{i_child}';
                        m.appearance.body{I}((k-1)*3+1:k*3) = {'cyl', [0 0 0; p_xyz], scale*1};
                    end
                end
            end
            
            % add the foot shape (box)
            m.appearance.body{model.FindBody('bodyLAR')}(4:5) = { 'box', [-5 -3 -4; 5 3 -5]};

            m.camera.body = 0;% fixed camera
            m.camera.direction = [1 0 0.1];
            m.camera.locus = [0 0.5];
            
            %=== add parameters for floating-base
            m.appearance.body = {{}, {}, {}, {}, {}, m.appearance.body{:}};
%             if m.camera.body > 0
%                 m.camera.body = m.camera.body + 5;
%             end
                        
            model.featherstone = m;
        end
        
        %add contact points to the featherstone model
        function model = AddContactPoints(model)
            m = model.featherstone;
            X = 5;
            Y = 3;
            Z = -5;
            m.gc.point = [-X X -X X;-Y -Y Y -Y;Z Z Z Z];
            m.gc.body = 4*ones(1,4);
            model.featherstone = m;
        end
    end
    
end



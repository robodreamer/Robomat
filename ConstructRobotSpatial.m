%% Superclass of Spatial Rigid Body (FB) Robot

%@written on 06-21-14
%@author: andypark at purdue . edu

classdef ConstructRobotSpatial
    
    properties (Abstract = true, SetAccess=protected)
        name; % robot name
        NB; % number of bodies:
        N_joint; % number of joints (including floating-body)
        N_ext; % total number of bodies + dummy bodies (end-effector frames)
        NJ;
        ag_0;
        display_mass;
        display_inertia;
    end
    
    properties (Constant = true)
        base_motion = 1;
    end
    
    properties (Access=public)
        parent;
        phi;
        i_p_com;
        m;
        I_com;
        pi_p_i;
        pi_R_i;
        joints;
        q;
        X_J;
        i_XL_pi;
        i_X_pi;
        T0;
        R0;
        I;
        I_C;
        i_X_0;
        pi_T_i;
        R0_6;
        v;
        a;
        chi;
        f;
        a0;
        q_a;
        q_p;
        bodyname;
        featherstone = []; % featherstone v2 library compatible
        %         featherstoneV1 = []; % featherstone v1 library compatible
        
    end
    
    %     properties (Acce
    
    % global A q qd CandG
    
    methods (Static = true)
        function M = crossSkew(x)
            M = [0 -x(3) x(2) ; x(3) 0 -x(1) ; -x(2) x(1) 0];
        end
        
        % from rotation matrix to 6x6 matrix
        function X = r2X(R)
            X = [R eye(3);eye(3) R];
        end
        
        % cross product operation for spatial motion
        function M = crossSpv(spV)
            w = spV(1:3);
            v = spV(4:6);
            
            M = [ConstructRobotSpatial.crossSkew(w) zeros(3);ConstructRobotSpatial.crossSkew(v) ConstructRobotSpatial.crossSkew(w)];
        end
        
        % cross product operation for spatial force
        function M = crossSpf(spV)
            M = -ConstructRobotSpatial.crossSpv(spV)';
        end
        
        function plotRows(A)
            plot3(A(1,:), A(2,:), A(3,:));
        end
        
        function scatterRows(A)
            scatter3(A(1,:), A(2,:), A(3,:));
        end
        
        function inertia = Inertia(I)
            % Inertia_robot: generate an inertia matrix from 1x6 row vector definition.
            % written on 02-10-13 by Andy Park
            
            Ixx = I(1);
            Iyy = I(2);
            Izz = I(3);
            Ixy = I(4);
            Iyz = I(5);
            Ixz = I(6);
            
            %Inertia = [Ixx -Ixy -Ixz; -Ixy Iyy -Iyz; -Ixz -Iyz Izz];
            inertia = [Ixx Ixy Ixz; Ixy Iyy Iyz; Ixz Iyz Izz];
            
        end
        
        function  dXj = djcalcV2( jtyp, q )
            % djcalc  joint transform and motion subspace differential matrices.
            % [Xj,S]=jcalc(type,q)  returns the joint transform and motion subspace
            % matrices for a joint of the given type.  jtyp can be either a string or a
            % structure containing a string-valued field called 'code'.  Either way,
            % the string contains the joint type code.  For joints that take
            % parameters (e.g. helical), jtyp must be a structure, and it must contain
            % a field called 'pars', which in turn is a structure containing one or
            % more parameters.  (For a helical joint, pars contains a parameter called
            % 'pitch'.)  q is the joint's position variable.
            
            % modified on 06-24-14 by andy park
            
            switch jtyp
                case 'FB'				% floating-base joint
                    
                case 'Rx'				% revolute X axis
                    dXj = dXrotx(q);
                case 'Ry'				% revolute Y axis
                    dXj = dXroty(q);
                case 'Rz'			% revolute Z axis
                    dXj = dXrotz(q);
                case 'Px'				% prismatic X axis
                    dXj = dXtrans([q 0 0]);
                    dXj = reshape(dXj(:,1),6,6);
                case 'Py'				% prismatic Y axis
                    dXj = dXtrans([0 q 0]);
                    dXj = reshape(dXj(:,2),6,6);
                case 'Pz'			% prismatic Z axis
                    dXj = dXtrans([0 0 q]);
                    dXj = reshape(dXj(:,2),6,6);
                otherwise
                    error( 'unrecognised joint code ''%s''', code );
            end
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%% functions added after July-2014 %%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % 07-09-14 compute the orientation error given two quaternions
        function err_o = compute_error_quat(qt_curr, qt_des);
            err_o = -Quat(qt_des)'*qt_curr;
        end
        
        % 07-08-14
        % Convert homogeneous/rotation transform to a unit-quaternion
        function qt = r2q(t, method_convert)
            if ~exist('method_convert','var')
                method_convert = 3; %by default
            end
            
            method_convert = 3;
            % Peter Corke's method
            if(method_convert == 1)
                qs = sqrt(trace(t)+1)/2.0;
                kx = t(3,2) - t(2,3);	% Oz - Ay
                ky = t(1,3) - t(3,1);	% Ax - Nz
                kz = t(2,1) - t(1,2);	% Ny - Ox
                
                if (t(1,1) >= t(2,2)) & (t(1,1) >= t(3,3))
                    kx1 = t(1,1) - t(2,2) - t(3,3) + 1;	% Nx - Oy - Az + 1
                    ky1 = t(2,1) + t(1,2);			% Ny + Ox
                    kz1 = t(3,1) + t(1,3);			% Nz + Ax
                    add = (kx >= 0);
                elseif (t(2,2) >= t(3,3))
                    kx1 = t(2,1) + t(1,2);			% Ny + Ox
                    ky1 = t(2,2) - t(1,1) - t(3,3) + 1;	% Oy - Nx - Az + 1
                    kz1 = t(3,2) + t(2,3);			% Oz + Ay
                    add = (ky >= 0);
                else
                    kx1 = t(3,1) + t(1,3);			% Nz + Ax
                    ky1 = t(3,2) + t(2,3);			% Oz + Ay
                    kz1 = t(3,3) - t(1,1) - t(2,2) + 1;	% Az - Nx - Oy + 1
                    add = (kz >= 0);
                end
                
                if add
                    kx = kx + kx1;
                    ky = ky + ky1;
                    kz = kz + kz1;
                else
                    kx = kx - kx1;
                    ky = ky - ky1;
                    kz = kz - kz1;
                end
                nm = norm([kx ky kz]);
                
                if nm == 0,
                    qt = [1 0 0 0]';
                else
                    s = (sqrt(1 - qs^2) / nm);
                    qv = s*[kx ky kz];
                    
                    qt = [qs qv]';
                end
                
                % text book
            elseif(method_convert == 2)
                R = t;
                
                Rxx = R(1,1); Rxy = R(1,2); Rxz = R(1,3);
                Ryx = R(2,1); Ryy = R(2,2); Ryz = R(2,3);
                Rzx = R(3,1); Rzy = R(3,2); Rzz = R(3,3);
                
                w = sqrt( trace( R ) + 1 ) / 2;
                
                % check if w is real. Otherwise, zero it.
                if( imag( w ) > 0 )
                    w = 0;
                end
                
                x = sqrt( 1 + Rxx - Ryy - Rzz ) / 2;
                y = sqrt( 1 + Ryy - Rxx - Rzz ) / 2;
                z = sqrt( 1 + Rzz - Ryy - Rxx ) / 2;
                
                [element, ind_tmp ] = max( [w,x,y,z] );
                
                if( ind_tmp == 1 )
                    x = ( Rzy - Ryz ) / (4*w);
                    y = ( Rxz - Rzx ) / (4*w);
                    z = ( Ryx - Rxy ) / (4*w);
                end
                
                if( ind_tmp == 2 )
                    w = ( Rzy - Ryz ) / (4*x);
                    y = ( Rxy + Ryx ) / (4*x);
                    z = ( Rzx + Rxz ) / (4*x);
                end
                
                if( ind_tmp == 3 )
                    w = ( Rxz - Rzx ) / (4*y);
                    x = ( Rxy + Ryx ) / (4*y);
                    z = ( Ryz + Rzy ) / (4*y);
                end
                
                if( ind_tmp == 4 )
                    w = ( Ryx - Rxy ) / (4*z);
                    x = ( Rzx + Rxz ) / (4*z);
                    y = ( Ryz + Rzy ) / (4*z);
                end
                
                Q = [ w; x; y; z ];
                qt = Q;
                
                % Roy Featherstone's method
            elseif(method_convert == 3)
                % for sufficiently large q0, this function formulates 2*q0 times the
                % correct return value; otherwise, it formulates 4*|q1| or 4*|q2| or 4*|q3|
                % times the correct value.  The final normalization step yields the correct
                % return value.
                E = t(1:3,1:3);
                tr = trace(E);				% trace is 4*q0^2-1
                v_ = -skew(E);				% v is 2*q0 * [q1;q2;q3]
                
                if tr > 0
                    qt = [ (tr+1)/2; v_ ];
                else
                    E = E - (tr-1)/2 * eye(3);
                    E = E + E';
                    if E(1,1) >= E(2,2) && E(1,1) >= E(3,3)
                        qt = [ 2*v_(1); E(:,1) ];
                    elseif E(2,2) >= E(3,3)
                        qt = [ 2*v_(2); E(:,2) ];
                    else
                        qt = [ 2*v_(3); E(:,3) ];
                    end
                    if qt(1) < 0
                        qt = -qt;
                    end
                end
                
                qt = qt / norm(qt);
                
                % reverse the sign
                qt(2:4) = -qt(2:4);
            end
        end
        
        % 07-08-14 compute the conversion matrix (q->w)
        function Q = Quat(qt)
            q0 = qt(1);
            q1 = qt(2:4);
            Q = [-q1'; q0*eye(3) + skew(q1)];
        end
        
        % 06-28-14 reverse of the tr2diff (w comes first)
        % this only works for small difference for rotational part
        function d = T2diff(t1,t2)
            if nargin == 1,
                d = [0.5*[t1(3,2)-t1(2,3); t1(1,3)-t1(3,1); t1(2,1)-t1(1,2)];t1(1:3,4)];
            else
                d = [0.5*(cross(t1(1:3,1), t2(1:3,1)) + ...
                    cross(t1(1:3,2), t2(1:3,2)) + ...
                    cross(t1(1:3,3), t2(1:3,3)) ...
                    ); t2(1:3,4)-t1(1:3,4)];
            end
        end
        
        % 06-28-14 take the passive part of the vector/matrix
        function out = passive(in)
            if(size(in,2) == 1) % only one column (e.g. q,qd,qdd)
                out = in(1:6,1);
            else % if no. column > 1 (e.g. J)
                out = in(:,1:6);
            end
        end
        
        % 06-28-14 take the actuated part of the vector/matrix
        function out = active(in)
            if(size(in,2) == 1) % only one column (e.g. q,qd,qdd)
                out = in(7:end,1);
            else % if no. column > 1 (e.g. J)
                out = in(:,7:end);
            end
        end
        
        % 06-28-14 take the positional part of the vector/matrix
        function out = TakePos(in)
            out = in(4,6,:);
        end
        
        % 06-28-14 take the rotational part of the vector/matrix
        function out = TakeOri(in)
            out = in(1,3,:);
        end
        
        % 06-27-14 cross product with multiple vectors
        % a and b have to be column vectors (a=3x1, b=3xn)
        function out = crossmany(a, b)
            n = size(b,2);
            out = zeros(3,n);
            for i=1:n
                out(:,i) = cross(a,b(:,i));
            end
        end
        
        
        % 06-27-14 compute the support vectors for the friction cone constraints
        function Betas = computeSupportVectorsFC(nSV, mu)
            if ~exist('nSV','var')
                nSV = 4; %by default (polyhedral cone linear approximation)
            end
            if ~exist('mu','var')
                mu = 0.6; %by default (alpha approx. 30 deg, thus 60 deg cone)
            end
            
            % angle of friction cone
            alpha = atan(mu);
            % radius of the top part
            r = sin(alpha);
            % height of the cone
            h = cos(alpha);
            
            % compute linear approximation (polyhedra) of the cone
            Betas = zeros(3,nSV);
            for k=1:nSV
                theta = 2*pi*(k-1)/nSV;
                Betas(:,k) = [r*cos(theta), r*sin(theta), h]';
            end
        end
        
        % % 06-27-14 compute the support vectors for the friction cone constraints
        %     function Betas = computeSupportVectorsFC_old(pts, nSV, mu)
        %         if ~exist('nSV','var')
        %             nSV = 4; %by default (polyhedral cone linear approximation)
        %         end
        %         if ~exist('mu','var')
        %             mu = 0.6; %by default (alpha approx. 30 deg, thus 60 deg cone)
        %         end
        %
        %         alpha = atan(mu);
        %         r = sin(alpha);
        %         h = cos(alpha);
        %
        %         if(size(pts,2) ~= 3) % if the pts dimension incorrect
        %             pts = pts'; % try transpose
        %         end
        %         nPts = size(pts,1);
        %
        %         if(nPts > 1) % multiple contact points
        %             Betas = cell(nPts,1);
        %
        %             for i=1:nPts
        %                 % initialize support vectors
        %                 Betas{i} = zeros(3,nSV);
        %                 pt_i = pts(i,:);
        %                 for k=1:nSV
        %                     theta = 2*pi*(k-1)/nSV;
        %                     %                     disp(theta)
        %                     B_i = pt_i + [r*cos(theta), r*sin(theta), h];
        %                     Betas{i}(:,k) = B_i';
        %                 end
        %             end
        %         elseif(nPts == 1)% one contact points
        %             Betas = zeros(3,nSV);
        %
        %             pt_i = pts;
        %             for k=1:nSV
        %                 theta = 2*pi*(k-1)/nSV;
        %                 B_i = pt_i + [r*cos(theta), r*sin(theta), h];
        %                 Betas(:,k) = B_i';
        %             end
        %         end
        %     end
        
        
        % 06-27-14 draw a linear approximation of the friction cone
        % Betas are support vectors around the cone so this translates to the point
        function varargout = meshFC(Betas, pt)
            %         figure(1); clf;
            %         [v f] = createTetrahedron;
            %         drawMesh(v, f);
            %         view(3); axis('vis3d'); axis off;
            %         title('Tetrahedron');
            % get number of nodes (nSV+1)
            nSV = size(Betas,2);
            if(nSV > 1)
                nodes = zeros(nSV+1, 3);
                % first node : contact point
                nodes(1,:) = pt;
                edges = [];
                faces = [];
                % build nodes
                for k=1:nSV
                    nodes(k+1,:) = pt + Betas(:,k)';
                end
                
                % build edges, faces
                edges = zeros(nSV,2);
                faces = zeros(nSV,3);
                for k=1:nSV
                    edges(k,:) = [1 k+1];
                    if(k < nSV)
                        faces(k,:) = [1 k+1 k+2];
                    else % k=nSV
                        faces(k,:) = [1 faces(1,2) faces(k-1,3)];
                    end
                end
                
                % format output
                varargout = formatMeshOutput(nargout, nodes, edges, faces);
            end
        end
        
        
        function varargout = drawPolyhedron_andy(parent, nodes, faces, varargin)
            %DRAWPOLYHEDRON Draw polyhedron defined by vertices and faces
            %
            %   drawPolyhedron(NODES, FACES)
            %   Draws the polyhedron defined by vertices NODES and the faces FACES.
            %   NODES is a NV-by-3 array containing coordinates of vertices, and FACES
            %   is either a NF-by3 or NF-by-4 array containing indices of vertices of
            %   the triangular or rectangular faces.
            %   FACES can also be a cell array, in the content of each cell is an array
            %   of indices to the nodes of the current face. Faces can have different
            %   number of vertices.
            %
            %   H = drawPolyhedron(...);
            %   Also returns a handle to the created patche.
            %
            %   Example:
            %   [n f] = createSoccerBall;
            %   drawPolyhedron(n, f);
            %
            %   See also:
            %   polyhedra, drawMesh, drawPolygon
            %
            %   ---------
            %
            %   author : David Legland
            %   INRA - TPV URPOI - BIA IMASTE
            %   created the 10/02/2005.
            %
            
            %   HISTORY
            %   07/11/2005 update doc.
            %   04/01/2007 typo
            %   18/01/2007 add support for 2D polyhedra ("nodes" is N-by-2 array), and
            %       make 'cnodes' a list of points instead of a list of indices
            %   14/08/2007 add comment, add support for NaN in faces (complex polygons)
            %   14/09/2007 rename as drawPolyhedron
            %   16/10/2008 better support for colors
            %   27/07/2010 copy to 'drawMesh'
            
            
            %% Initialisations
            
            
            % process input arguments
            switch length(varargin)
                case 0
                    % default color is red
                    varargin = {'facecolor', [1 0 0]};
                case 1
                    % use argument as color for faces
                    varargin = {'facecolor', varargin{1}};
                otherwise
                    % otherwise do nothing
            end
            
            % overwrites on current figure
            hold on;
            
            % if nodes are 2D points, add a z=0 coordinate
            if size(nodes, 2) == 2
                nodes(1,3) = 0;
            end
            
            
            %% main loop : for each face
            
            if iscell(faces)
                % array FACES is a cell array
                h = zeros(length(faces(:)), 1);
                
                for f = 1:length(faces(:))
                    % get nodes of the cell
                    face = faces{f};
                    
                    if sum(isnan(face))~=0
                        % Special processing in case of multiple polygonal face.
                        % each polygonal loop is separated by a NaN.
                        
                        % find indices of loops breaks
                        inds = find(isnan(face));
                        
                        % replace NaNs by index of first vertex of each polygon
                        face(inds(2:end))   = face(inds(1:end-1)+1);
                        face(inds(1))       = face(1);
                        face(length(face)+1)= face(inds(end)+1);
                    end
                    
                    % draw current face
                    cnodes  = nodes(face, :);
                    h(f)    = patch(cnodes(:, 1), cnodes(:, 2), cnodes(:, 3), [1 0 0]);
                    %                 h(f)    = patch('Parent',parent,cnodes(:, 1), cnodes(:, 2), cnodes(:, 3), [1 0 0]);
                end
                
            else
                % array FACES is a NC*NV indices array, with NV : number of vertices of
                % each face, and NC number of faces
                h = zeros(size(faces, 1), 1);
                for f = 1:size(faces, 1)
                    % get nodes of the cell
                    cnodes = nodes(faces(f,:)', :);
                    h(f) = patch(cnodes(:, 1), cnodes(:, 2), cnodes(:, 3), [1 0 0]);
                    %                 h(f) = patch('Parent',parent, cnodes(:, 1), cnodes(:, 2), cnodes(:, 3), [1 0 0]);
                end
            end
            
            % set up drawing options
            if ~isempty(varargin)
                set(h, varargin{:});
            end
            
            % format output parameters
            if nargout > 0
                varargout = {h};
            end
            
            % % 06-27-14 draw a linear approximation of the friction cone
            %     function varargout = meshFC_old(Betas, pt)
            %         %         figure(1); clf;
            %         %         [v f] = createTetrahedron;
            %         %         drawMesh(v, f);
            %         %         view(3); axis('vis3d'); axis off;
            %         %         title('Tetrahedron');
            %         % get number of nodes (nSV+1)
            %         nSV = size(Betas,2);
            %         if(nSV > 1)
            %             nodes = zeros(nSV+1, 3);
            %             % first node : contact point
            %             nodes(1,:) = pt;
            %             edges = [];
            %             faces = [];
            %             % build nodes
            %             for k=1:nSV
            %                 nodes(k+1,:) = Betas(:,k)';
            %             end
            %
            %             % build edges, faces
            %             edges = zeros(nSV,2);
            %             faces = zeros(nSV,3);
            %             for k=1:nSV
            %                 edges(k,:) = [1 k+1];
            %                 if(k < nSV)
            %                     faces(k,:) = [1 k+1 k+2];
            %                 else % k=nSV
            %                     faces(k,:) = [1 faces(1,2) faces(k-1,3)];
            %                 end
            %             end
            %
            %             % format output
            %             varargout = formatMeshOutput(nargout, nodes, edges, faces);
            %         end
            %     end
            
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%% useful functions for class object %%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    methods
        
        function result = isInSubtree(model,leaf,candidate)
            parent = model.parent{leaf};
            while parent ~= 0
                if parent == candidate
                    result = 1;
                    return
                else
                    parent = model.parent{parent};
                end
            end
            result = 0;
        end
        
        % compute jacobian
        function jac = jacob(model,link,coord)
            if ~exist('coord','var')
                coord = 0; %world coordinate by default
            end
            jac = jacob_pt(model,link,[0 0 0]',coord);
            %             if ~exist('coord','var')
            %                 coord = 0; %world coordinate by default
            %             end
            %
            %             % Initalize Jacobian
            %             jac = zeros(6,model.N_joint);
            %
            %             % Recursively build jacobian
            %             if link > model.NB % for virtual body (frames)
            %                 parent = model.parent{link};
            %                 c_X_j = model.i_X_pi{link};
            %             else
            %                 c_X_j = eye(6);
            %                 parent = link;
            %             end
            %
            %             while parent ~= 0
            %                 j = parent;
            %                 %                 disp (j);
            %                 % c_X_j: transform from the frame of interest to the each joint frame
            %                 % since phi is defined w.r.t each joint frame
            %                 if(j==1)
            %                     % only for floating base the world frame is where the
            %                     % joints are defined
            %                     joints_0 = model.joints{j};
            %                     phi_0 = model.phi{j};
            %                     jac(:,joints_0(1:3)) = model.i_X_pi{j}*phi_0(:,1:3);
            %                     jac(:,joints_0(4:6)) = model.i_X_pi{j}*phi_0(:,4:6);
            %                 else
            %                     jac(:,model.joints{j}) = c_X_j*model.phi{j};
            %                     c_X_j = c_X_j*model.i_X_pi{j};
            %                 end
            %                 %                 cxj =[model.R0{link} zeros(3,3); zeros(3,3) model.R0{link}]*c_X_j;
            %                 %                 j0 = [model.R0{link} zeros(3,3); zeros(3,3) model.R0{link}]*jac;
            %
            %                 parent = model.parent{j};
            %             end
            %
            %             if(coord == 0)
            %                 %Transform to global
            %                 jac = [model.R0{link} zeros(3,3); zeros(3,3) model.R0{link}]*jac;
            %             else
            %             end
        end
        
        % 06-25-14 compute jacobian dot
        function jacdot = jacob_dot(model,link,qd,coord)
            %             if ~exist('coord','var')
            %                 coord = 0; %world coordinate by default
            %             end
            %             jacdot = jacob_pt_dot(model,link,qd,[0 0 0]',coord);
            if ~exist('coord','var')
                coord = 0; %world coordinate by default
            end
            % Initalize Jacobian
            jacdot = zeros(6,model.N_joint);
            
            % compute the velocity of each body
            %% RNEA outward recursion
            for i=1:model.NB
                p = model.parent{i};
                if p == 0
                    v{i} = model.phi{i}*qd(model.joints{i});
                else
                    v{i} = model.i_X_pi{i}*v{p} + model.phi{i}*qd(model.joints{i});
                end
            end
            
            % Recursively build jacobian
            if link > model.NB % for virtual body (frames)
                parent = model.parent{link};
                c_X_j = model.i_X_pi{link};
            else
                c_X_j = eye(6);
                parent = link;
            end
            
            while parent ~= 0
                j = parent;
                % v_i x si*qd_i
                if(j == 1)
                    %                     joints_0 = model.joints{j};
                    %                     phi_0 = model.phi{j};
                    %                     jacdot(:,joints_0(4:6)) = c_X_j*model.crossSpv(v{j})*phi_0(:,4:6);
                    c_X_j = c_X_j*model.i_X_pi{j};
                    %                     jacdot(:,joints_0(1:3)) = c_X_j*model.crossSpv(v{j})*phi_0(:,1:3);
                    %                     c_X_j = c_X_j*model.i_X_pi{j};
                    jacdot(:,model.joints{j}) = c_X_j*model.crossSpv(v{j})*model.phi{j};
                else
                    jacdot(:,model.joints{j}) = c_X_j*model.crossSpv(v{j})*model.phi{j};
                    c_X_j = c_X_j*model.i_X_pi{j};
                end
                parent = model.parent{j};
            end
            
            if(coord == 0)
                %Transform to global
                jacdot = [model.R0{link} zeros(3,3); zeros(3,3) model.R0{link}]*jacdot;
            else
            end
        end
        
        % compute jacobian of a certain position of the body i
        function jac = jacob_pt(model,link,pt,coord)
            if ~exist('coord','var')
                coord = 0; %world coordinate by default
            end
            % Initalize Jacobian
            jac = zeros(6,model.N_joint);
            
            % Recursively build jacobian
            if link > model.NB % for virtual body (frames)
                parent = model.parent{link};
                c_X_j = model.i_X_pi{link};
            else
                c_X_j = xlt(pt);
                parent = link;
            end
            
            while parent ~= 0
                j = parent;
                if(j==1)
                    %                     joints_0 = model.joints{j};
                    %                     phi_0 = model.phi{j};
                    c_X_j = c_X_j*model.i_X_pi{j};
                    jac(:,model.joints{j}) = c_X_j*model.phi{j};
                    %                     jac(:,joints_0(1:3)) = c_X_j*phi_0(:,1:3);
                    %                     jac(:,joints_0(4:6)) = c_X_j*phi_0(:,4:6);
                    %                     c_X_j = c_X_j*model.i_X_pi{j};
                    %                     jac(:,model.joints{j}) = c_X_j*model.phi{j};
                else
                    jac(:,model.joints{j}) = c_X_j*model.phi{j};
                    c_X_j = c_X_j*model.i_X_pi{j};
                end
                parent = model.parent{j};
            end
            
            if(coord == 0)
                %Transform to global
                jac = [model.R0{link} zeros(3,3); zeros(3,3) model.R0{link}]*jac;
            else
            end
        end
        
        % 06-25-14 compute jacobian dot for a point in the body
        function jacdot = jacob_pt_dot(model,link,pt,qd,coord)
            if ~exist('coord','var')
                coord = 0; %world coordinate by default
            end
            % Initalize Jacobian
            jacdot = zeros(6,model.N_joint);
            
            % compute the velocity of each body
            %% RNEA outward recursion
            for i=1:model.NB
                p = model.parent{i};
                if p == 0
                    v{i} = model.phi{i}*qd(model.joints{i});
                else
                    v{i} = model.i_X_pi{i}*v{p} + model.phi{i}*qd(model.joints{i});
                end
            end
            
            % Recursively build jacobian
            if link > model.NB % for virtual body (frames)
                parent = model.parent{link};
                c_X_j = model.i_X_pi{link};
            else
                %                 c_X_j = eye(6);
                c_X_j = xlt(pt);
                parent = link;
            end
            
            while parent ~= 0
                j = parent;
                % v_i x si*qd_i
                if(j==1)
                    %                     joints_0 = model.joints{j};
                    %                     phi_0 = model.phi{j};
                    c_X_j = c_X_j*model.i_X_pi{j};
                    jacdot(:,model.joints{j}) = c_X_j*model.crossSpv(v{j})*model.phi{j};
                    %                     jacdot(:,joints_0(1:3)) = c_X_j*model.crossSpv(v{j})*phi_0(:,1:3);
                    %                     jacdot(:,joints_0(4:6)) = c_X_j*model.crossSpv(v{j})*phi_0(:,4:6);
                    %                     c_X_j = c_X_j*model.i_X_pi{j};
                    %                     jacdot(:,model.joints{j}) = c_X_j*model.crossSpv(v{j})*model.phi{j};
                else
                    jacdot(:,model.joints{j}) = c_X_j*model.crossSpv(v{j})*model.phi{j};
                    c_X_j = c_X_j*model.i_X_pi{j};
                end
                parent = model.parent{j};
            end
            
            if(coord == 0)
                %Transform to global
                jacdot = [model.R0{link} zeros(3,3); zeros(3,3) model.R0{link}]*jacdot;
            else
            end
        end
        
        %06-21-14 compute the CoM
        function Pcom = getCoM(model)
            c = 0;
            mb = 0;
            for i = 1:model.NB
                mb = mb + model.m{i};
                c_i = model.T0{i}*[model.i_p_com{i} ;1];
                c = c + model.m{i}*c_i;
            end
            Pcom = c(1:3)/mb;
        end
        
        % compute COM jacobian
        function Jcom = getJcom(model)
            % Initalize Jacobian
            
            jac_com = zeros(6,model.N_joint);
            
            if(1) % more concise way to compute Jmi
                
                mb = 0;
                for link=1:model.NB
                    jac_mi = model.jacob_pt(link, model.i_p_com{link});
                    mb = mb + model.m{link};
                    jac_com = jac_com + model.m{link}*jac_mi;
                end
                Jcom = jac_com/mb;
                
                % J CoM w.r.t world frame
                %                 if(1)
                %                     Jcom = model.R0_6{1}'*Jcom;
                %                 end
                
            else % original implementation
                mb = 0;
                for link=1:model.NB
                    
                    jac_mi = zeros(6,model.N_joint);
                    mb = mb + model.m{link};
                    
                    % Recursively build jacobian
                    if link > model.NB %for dummy body
                        parent = model.parent{link};
                        c_X_j = model.i_X_pi{link};
                    else
                        %                     c_X_j = eye(6);
                        c_X_j = xlt(model.i_p_com{link});
                        parent = link;
                    end
                    
                    while parent ~= 0
                        j = parent;
                        %                 disp (j);
                        jac_mi(:,model.joints{j}) = c_X_j*model.phi{j};
                        % brings the transform w.r.t the end-frame
                        c_X_j = c_X_j*model.i_X_pi{j};
                        parent = model.parent{j};
                    end
                    
                    %Transform to global
                    jac_mi = [model.R0{link} zeros(3,3); zeros(3,3) model.R0{link}]*jac_mi;
                    %                 % verify the jacob_pt result
                    %                 valuecheck(model.jacob_pt(link, model.i_p_com{link}),jac_mi);
                    %                 jac_mi(4:end,7:end)
                    jac_com = jac_com + model.m{link}*jac_mi;
                end
                %             jac_com
                %             mb
                Jcom = jac_com/mb;
            end
        end
        
        % compute COM jacobian
        function Jcom_dot = getJcom_dot(model,qd)
            % Initalize Jacobian
            jac_com_dot = zeros(6,model.N_joint);
            
            mb = 0;
            for link=1:model.NB
                jac_mi_dot = model.jacob_pt_dot(link, model.i_p_com{link}, qd);
                mb = mb + model.m{link};
                jac_com_dot = jac_com_dot + model.m{link}*jac_mi_dot;
            end
            Jcom_dot = jac_com_dot/mb;
        end
        
        
        % find the index of the body by name
        function index = FindBody(model,str)
            a = model.bodyname;
            ind = ind2sub(size(a),find(cellfun(@(x)strcmp(x,str),a)));
            % error if the body name is incorrect
            if(isempty(ind))
                cprintf('Errors','link name is incorrect\n');
            else
                index = ind;
            end
        end
        
        % 07-10-14 compute the position of a point in a body
        function pos = posBody(model,link,pt)
            if ~exist('pt','var')
                pt = zeros(3,1); %world coordinate by default
            end
            pt = [pt;1];
            T0 = model.T0{link};
            pos = T0*pt;
            pos = pos(1:3);
        end
        
        % 06-26-14 compute dXJ for joints
        function di_X_pi = dXup(model,i,qd)
            if(i==1) % at floating-base
                %                 jtype = 'FB';
                q = model.q{i};
                %                 dRx = ConstructRobotSpatial.djcalcV2('Rx', q(4));
                %                 dRy = ConstructRobotSpatial.djcalcV2('Ry', q(5));
                %                 dRz = ConstructRobotSpatial.djcalcV2('Rz', q(6));
                dRx = djcalcV2('Rx', q(4));
                dRy = djcalcV2('Ry', q(5));
                dRz = djcalcV2('Rz', q(6));
                Rx = jcalcV2('Rx', q(4));
                Ry = jcalcV2('Ry', q(5));
                Rz = jcalcV2('Rz', q(6));
                dRx = dRx(1:3,1:3);
                dRy = dRy(1:3,1:3);
                dRz = dRz(1:3,1:3);
                Rx = Rx(1:3,1:3);
                Ry = Ry(1:3,1:3);
                Rz = Rz(1:3,1:3);
                
                P_base = skew(q(1:3));
                dP_base = skew(qd(1:3));
                
                %                 Rup = Rz*Ry*Rx;
                %                 dRup = qd(6)*dRz*Ry*Rx + qd(5)*Rz*dRy*Rx + qd(4)*Rz*Ry*dRx;
                Rup = Rx*Ry*Rz;
                dRup = qd(4)*dRx*Ry*Rz + qd(5)*Rx*dRy*Rz + qd(6)*Rx*Ry*dRz;
                
                %                 di_X_pi = [dRup, zeros(3); -dP_base*Rup-P_base*dRup, dRup];
                di_X_pi = [dRup, zeros(3); -Rup*dP_base-dRup*P_base, dRup];
            else
                jointaxis = find(model.phi{i}(1:3) == 1);
                if(jointaxis == 1) % Rx
                    jtype = 'Rx';
                elseif(jointaxis == 2) % Rx
                    jtype = 'Ry';
                elseif(jointaxis == 3) % Rx
                    jtype = 'Rz';
                end
                dXidq = model.djcalcV2(jtype, model.q{i}); % modified by Andy Park (06-24-14)
                di_X_pi = dXidq * model.i_XL_pi{i} * qd(model.joints{i});
            end
        end
        
        
        %%% update the model with the joint configurations
        function model = UpdateModel(model,q,display)
            if ~exist('display','var')
                display = 0; %off by default
            end
            
            % actuated joints
            model.q_a = q(7:end);
            model.q_p = q(1:6);
            
            % base body (floating)
            K = 6;
            KBody = 1;
            if(model.base_motion == 1)
                model.pi_p_i{KBody}   = q(1:3); % position of the joint
                q_ = q(4:6);
                %                 R_base = expm(-ConstructRobotSpatial.crossSkew(ones(3,1).*q_));
                R_base = rotz(q_(3))*roty(q_(2))*rotx(q_(1));
                %                 R_base = rotx(q_(1))*roty(q_(2))*rotz(q_(3));
                model.pi_R_i{KBody} = R_base; % rotation of the joint
            end
            model.q{KBody}      = q(1:K); % initial joint angles
            
            for KBody=2:model.NB
                model.q{KBody}      = q(KBody+5);
            end
            
            if(display == 1)
                figure(1);
                clf;
                %                 plot3(0,0,0);
                % view(178,2)
                hold on
                axis square
                %                 axis vis3d
                x = model.pi_p_i{1}(1);
                y = model.pi_p_i{1}(2);
                z = model.pi_p_i{1}(3);
                % axis([x-1 x+1 y-1 y+1 z-1 z+1])
                daspect([1 1 1])
                xlabel('x');
                ylabel('y');
                zlabel('z');
                grid on;
            end
            
            
            %% Setup all neccecaries for model structs
            for i=1:model.N_ext
                %                 j  = model;
                pred = model.parent{i};
                
                % Setup Coordinates
                
                if i==1 % base bodies
                    % translation and rotation
                    %                     model.i_XL_pi{i} = [model.pi_R_i{i}' zeros(3); -skew(model.pi_p_i{i}) model.pi_R_i{i}'];
                    model.i_XL_pi{i} = eye(6);
                else % non-base body
                    model.i_XL_pi{i} = [model.pi_R_i{i}' zeros(3); -model.pi_R_i{i}'*skew(model.pi_p_i{i}) model.pi_R_i{i}'];
                end
                
                % rotation first translation next
                %                 model.i_XL_pi{i} = [model.pi_R_i{i}' zeros(3); -model.pi_R_i{i}'*skew(model.pi_p_i{i}) model.pi_R_i{i}'];
                % translation first rotation next
                %                 model.i_XL_pi{i} = [model.pi_R_i{i}' zeros(3); -skew(model.pi_p_i{i})*model.pi_R_i{i}' model.pi_R_i{i}'];
                %                 model.i_XL_pi{i} = xlt(model.pi_p_i{i})*model.r2X(model.pi_R_i{i}'); % translation -> rotation
                
                if i == 1 % floating-base
                    R = eye(3);
                    %                     model.X_J{i} = eye(6);
                    %                     model.X_J{i} = [model.pi_R_i{i}' zeros(3); -skew(model.pi_p_i{i}) model.pi_R_i{i}'];
                    model.X_J{i} = [model.pi_R_i{i}' zeros(3); -model.pi_R_i{i}'*skew(model.pi_p_i{i}) model.pi_R_i{i}'];
                    
                elseif i > model.NB %end-effector frame (no transform due to joint rotation)
                    R = eye(3);
                    model.X_J{i} = eye(6);
                    
                else %1 < i < model.NB
                    q_ = model.q{i};
                    
                    %                     R = expm(-ConstructRobotSpatial.crossSkew(model.phi{i}*q_));
                    %                     tmp = model.phi{i};
                    R = angvec2r(-q_, model.phi{i}(1:3)); % equivalent expression
                    %         sum(sum(R - R2))
                    model.X_J{i} = [R zeros(3); zeros(3) R];
                end
                
                % 6x6 Adj transform
                model.i_X_pi{i} = model.X_J{i} * model.i_XL_pi{i};
                
                % 4x4 HT
                if i == 1 % at floating-joint
                    model.pi_T_i{i} =  [ model.pi_R_i{i} model.pi_p_i{i} ; 0 0 0 1];
                else
                    model.pi_T_i{i} =  [ model.pi_R_i{i} model.pi_p_i{i} ; 0 0 0 1]*[ R' [0;0;0] ; 0 0 0 1];
                end
                %                 model.pi_T_i{1}
                
                % Initialize inertia
                I_com_bar = model.I_com{i} + model.m{i} * ConstructRobotSpatial.crossSkew(model.i_p_com{i}) * ConstructRobotSpatial.crossSkew(model.i_p_com{i})';
                % inertia w.r.t the i-th link coordinate
                model.I{i} = [I_com_bar model.m{i}*ConstructRobotSpatial.crossSkew(model.i_p_com{i}) ; model.m{i} * ConstructRobotSpatial.crossSkew(model.i_p_com{i})' model.m{i} * eye(3);];
                
                % Initialize global transforms
                if model.parent{i} == 0
                    model.T0{i} = model.pi_T_i{i};
                    model.R0{i} = model.pi_T_i{i}(1:3,1:3);
                    %                     model.R0{i}
                    model.i_X_0{i} = model.i_X_pi{i};
                else
                    model.T0{i} = model.T0{pred} * model.pi_T_i{i};
                    model.R0{i} = model.R0{pred} * model.pi_T_i{i}(1:3,1:3);
                    model.i_X_0{i} = model.i_X_pi{i}*model.i_X_0{pred};
                    
                    if(display == 1)
                        ConstructRobotSpatial.plotRows([model.T0{i}*[0 0 0 1]' model.T0{pred}*[0 0 0 1]']);
                        
                    end
                end
                % for jacobian transform
                model.R0_6{i} = [model.R0{i} zeros(3) ; zeros(3) model.R0{i}];
                
                if(display == 1 && model.display_mass == 1)
                    if i<=model.NB
                        % Plot Masses (only for the bodies that are moved by joints)
                        ConstructRobotSpatial.scatterRows(model.T0{i}*[model.i_p_com{i} ;1]);
                        
                        if(model.display_inertia == 1)
                            
                            [V D] = eig(model.I_com{i});
                            %D=sqrt(diag(D));
                            %D=D/10;
                            S = (model.m{i}/5*(ones(3)-eye(3)))\diag(D);
                            S = sqrt(S)/5*(model.m{i})^.25;
                            v1 = V(:,1)*S(1); v2 = V(:,2)*S(2); v3 = V(:,3)*S(3);
                            
                            theta = [0:.1:2*pi 0];
                            ct = cos(theta); st = sin(theta); onest = 0*ct+1;
                            
                            
                            % Plot inertias
                            ConstructRobotSpatial.plotRows(model.T0{i}*[(v1*ct+v2*st)+ model.i_p_com{i}*onest;onest])
                            ConstructRobotSpatial.plotRows(model.T0{i}*[(v1*ct+v3*st)+ model.i_p_com{i}*onest;onest])
                            ConstructRobotSpatial.plotRows(model.T0{i}*[(v3*ct+v2*st)+ model.i_p_com{i}*onest;onest])
                        end
                    end
                end
            end
            
        end
        
        % fixed-base inverse-dynamics
        function tauTotal = inverseDynamics(model,qdd,qd,f_ext,grav0)
            tauTotal = zeros(model.N_joint,1);
            
            %% RNEA outward recursion
            for i=1:model.NB
                %                 r = model;
                p = model.parent{i};
                if p == 0 % at base body
                    model.v{i}   = model.phi{i}*qd(model.joints{i});
                    model.chi{i} = ConstructRobotSpatial.crossSpv(model.v{i})*model.phi{i}*qd(model.joints{i});
                    model.a{i}   = model.i_X_pi{i}*(-grav0) + model.chi{i} + model.phi{i}*qdd(model.joints{i});
                    model.f{i}   = model.I{i}*model.a{i} +...
                        ConstructRobotSpatial.crossSpf(model.v{i})*model.I{i}*model.v{i} - f_ext{i};
                else
                    model.v{i} = model.i_X_pi{i}*model.v{p} + model.phi{i}*qd(model.joints{i});
                    model.chi{i} = ConstructRobotSpatial.crossSpv(model.v{i})*model.phi{i}*qd(model.joints{i});
                    model.a{i}   = model.i_X_pi{i}*model.a{p} + model.chi{i} + model.phi{i}*qdd(model.joints{i});
                    model.f{i}   = model.I{i}*model.a{i} + ...
                        ConstructRobotSpatial.crossSpf(model.v{i})*model.I{i}*model.v{i} - f_ext{i};
                end
            end
            %% rnea inward recursion
            for i=model.NB:-1:1
                p = model.parent{i};
                tauTotal(model.joints{i}) =  model.phi{i}'*model.f{i};
                if p>0
                    model.f{p} = model.f{p}+model.i_X_pi{i}'*model.f{i};
                end
            end
            %             model;
            
        end
        
        % restructured floating-base inverse dynamics
        function [a0_method3 tau_method3] = inverseDynamicsFB(model,qdd,qd,f_ext,grav0)
            % initialization for the floating base
            model.a{1} = -model.ag_0;
            model.I_C{1} = model.I{1};
            model.v{1} = model.phi{1}*qd(1:6);
            %             model.v{1} = qd(1:6); % because of this multiplication qd is reversed again!
            model.f{1} = model.I{1}* model.a{1} + ...
                ConstructRobotSpatial.crossSpf(model.v{1}) * model.I{1}*model.v{1} -...
                f_ext{1};
            
            % inverse dynamics with floating base accel ==0
            % outward rnea
            %             qdd = [zeros(6,1) ; qddJoints];
            for i = 2:model.NB
                p = model.parent{i};
                model.v{i} = model.i_X_pi{i}*model.v{p}+model.phi{i}*qd(model.joints{i});
                model.chi{i} = model.crossSpv(model.v{i})*model.phi{i}*qd(model.joints{i});
                model.a{i}   = model.i_X_pi{i}*model.a{p} + model.chi{i}+model.phi{i}*qdd(model.joints{i});
                model.f{i}   = model.I{i}*model.a{i} + ...
                    ConstructRobotSpatial.crossSpf(model.v{i})*model.I{i}*model.v{i} - f_ext{i};
            end
            % inward rnea
            tauTotal = zeros(model.N_joint,1);
            for i=model.NB:-1:2
                p = model.parent{i};
                tauTotal(model.joints{i}) =  model.phi{i}'*model.f{i};
                model.f{p} = model.f{p}+model.i_X_pi{i}'*model.f{i};
            end
            
            % Composite rigid body recursion for I_C of floating base
            for i = 1:model.NB
                model.I_C{i} = model.I{i};
            end
            for i = model.NB:-1:2
                p = model.parent{i};
                %             r = model(i);
                model.I_C{p} = model.I_C{p} + model.i_X_pi{i}' * model.I_C{i} * model.i_X_pi{i};
            end
            
            a0_method3 = (-model.I_C{1})\model.f{1};
            model.a0{1} = a0_method3;
            % Add in torques due to floating base acceleration
            for i = 2:model.NB
                p = model.parent{i};
                model.a0{i} = model.i_X_pi{i} * model.a0{p};
                % Add in additional torques
                tauTotal(model.joints{i}) = tauTotal(model.joints{i}) + model.phi{i}'*model.I_C{i}*model.a0{i};
            end
            tau_method3 = tauTotal(7:end);
        end
        
        % update spatial vector quantities in the system
        function model = UpdateSpatialVectors(model,qdd,qd,f_ext,grav0)
            
            %% RNEA outward recursion
            for i=1:model.NB
                %                 r = model;
                p = model.parent{i};
                
                if p == 0
                    model.v{i}   = model.phi{i}*qd(model.joints{i});
                    model.chi{i} = ConstructRobotSpatial.crossSpv(model.v{i})*model.phi{i}*qd(model.joints{i});
                    model.a{i}   = model.i_X_pi{i}*(-grav0) + model.chi{i} + model.phi{i}*qdd(model.joints{i});
                    
                    model.f{i}   = model.I{i}*model.a{i} +...
                        ConstructRobotSpatial.crossSpf(model.v{i})*model.I{i}*model.v{i} - f_ext{i};
                else
                    model.v{i} = model.i_X_pi{i}*model.v{p} + model.phi{i}*qd(model.joints{i});
                    model.chi{i} = ConstructRobotSpatial.crossSpv(model.v{i})*model.phi{i}*qd(model.joints{i});
                    model.a{i}   = model.i_X_pi{i}*model.a{p} + model.chi{i} + model.phi{i}*qdd(model.joints{i});
                    model.f{i}   = model.I{i}*model.a{i} + ...
                        ConstructRobotSpatial.crossSpf(model.v{i})*model.I{i}*model.v{i} - f_ext{i};
                end
            end
            
        end
        
        % compute the joint space mass matrix
        function H = computeH(model)
            %% Calculate H
            H = zeros(model.NB,model.NB);
            
            for i = 1:model.NB
                I_C{i} = model.I{i};
            end
            
            for i = model.NB:-1:1
                F = I_C{i}*model.phi{i};
                %model(i).I_C
                %fprintf('Link (%d,%d), (%f,%f,%f,%f,%f,%f)\NB',i,i,F(1),F(2),F(3),F(4),F(5),F(6));
                H(model.joints{i},model.joints{i}) = model.phi{i}'*F;
                if model.parent{i} ~= 0
                    pred = model.parent{i};
                    I_C{pred} = I_C{pred} + model.i_X_pi{i}'*I_C{i}*model.i_X_pi{i};
                end
                j = i;
                while model.parent{j} ~= 0
                    F = model.i_X_pi{j}'*F;
                    j = model.parent{j};
                    %fprintf('Link (%d,%d), (%f,%f,%f,%f,%f,%f)\NB',i,j,F(1),F(2),F(3),F(4),F(5),F(6));
                    H(model.joints{i},model.joints{j}) = F' * model.phi{j};
                    H(model.joints{j},model.joints{i}) = H(model.joints{i},model.joints{j})';
                end
            end
            
        end
        
        % 06-26-14 compute C and G component in the rigid body system
        function out = computeCandG(model,qd,method,grav)
            if ~exist('method','var')
                method = 0; % compute C and G both by default
            end
            if ~exist('grav','var')
                grav = model.ag_0; % compute C and G both by default
            end
            f_ext_zero = cell(model.NB,1);
            for i = 1:model.NB
                f_ext_zero{i} = zeros(6,1);
            end
            
            qdd = zeros(model.N_joint,1);
            
            if(method == 0) % compute C and G
                %                 grav = model.ag_0;
            elseif(method == 1) % compute C
                grav = zeros(6,1);
            end
            
            out = model.inverseDynamics(qdd,qd,f_ext_zero,model.ag_0);
        end
        
        % 06-26-14 compute Gravity torque component in the rigid body system
        function out = computeG(model,grav)
            if ~exist('grav','var')
                grav = model.ag_0; % use gravity from the model definition
            end
            f_ext_zero = cell(model.NB,1);
            qdd = zeros(model.N_joint,1);
            qd = zeros(model.N_joint,1);
            out = model.inverseDynamics(qdd,qd,f_ext_zero,grav);
        end
        
        % 06-20-14 compute the composite rigid body inertia for
        % floating-base
        function I_C_0 = CompositeInertia(model)
            % Composite rigid body recursion for I_C of floating base
            for i = 1:model.NB
                I_C{i} = model.I{i};
            end
            
            I_C_0 = 0;
            for i = model.NB:-1:1
                p = model.parent{i};
                %             r = model;
                if p > 0
                    I_C{p} = I_C{p} + model.i_X_pi{i}' * I_C{i} * model.i_X_pi{i};
                else %p=0 at floating base
                    I_C_0 = I_C_0 + model.i_X_pi{i}' * I_C{i} * model.i_X_pi{i};
                end
            end
        end
        
        % 06-21-14 compute the whole-body mass
        function M_total = getMass(model)
            % Composite rigid body recursion for I_C of floating base
            mb = 0;
            for i = 1:model.NB
                mb = mb + model.m{i};
            end
            
            M_total = mb;
        end
        
        
        
        %         % compute COM jacobian
        %         function Pcom = getCoMQ(model,q)
        %             model.UpdateModel(q);
        %
        %             c = 0;
        %             mb = 0;
        %             for i = 1:model.NB
        %                 mb = mb + model.m{i};
        %                 c_i = model.T0{i}*[model.i_p_com{i} ;1];
        %                 c = c + model.m{i}*c_i;
        %             end
        %             Pcom = c(1:3)/mb;
        %
        %         end
        
        % 06-25-14 compute Centroidal Momentum matrix
        function [A_G1, h_G1, I_G1, v_G1] = computeCMM(model,qd, method)
            if ~exist('method','var')
                method = 1; % use recursive method by default
            end
            
            if(method == 1) % recursive
                I_C = cell(model.NB,1);
                I_C_0 = 0;
                for i = 1:model.NB
                    I_C{i} = model.I{i};
                end
                for i = model.NB:-1:1
                    p = model.parent{i};
                    if (p~=0)
                        I_C{p} = I_C{p} + model.i_X_pi{i}' * I_C{i} * model.i_X_pi{i};
                    else % floating-base
                        I_C_0 = I_C_0 + model.i_X_pi{i}' * I_C{i} * model.i_X_pi{i};
                    end
                end
                
                h_G1 = 0;
                % we need minus sign
                X_G = xlt(-skew(I_C_0(1:3,4:6))/I_C_0(6,6)); % 0_X_G
                i_X_G = cell(model.NB,1);
                i_X_G_ = cell(model.NB,1);
                A_G_i1 = cell(model.NB,1);
                A_G1 = [];
                
                for i = 1:model.NB
                    % compute pi_X_G (recursive) use the property that
                    % parent body index < children body indices
                    p = model.parent{i};
                    if(p==0) % at floating base
                        i_X_G{i} = model.i_X_pi{i} * X_G;
                    else
                        i_X_G{i} = model.i_X_pi{i} * i_X_G{p};
                    end
                    % %                     % another way to compute i_X_G (non-recursive way)
                    % %                     i_X_G_{i} = model.i_X_0{i}*X_G;
                    % %                     valuecheck(i_X_G{i},i_X_G_{i});
                    
                    %                     A_G_i1{i} = i_X_G{i}'*I_C{i}*model.phi{i};
                    if(i==1) % at floating base joint
                        A_G_i1{i} = i_X_G{i}'*I_C{i}*model.i_X_pi{i}*model.phi{i};
                    else
                        A_G_i1{i} = i_X_G{i}'*I_C{i}*model.phi{i};
                    end
                    h_G1 = h_G1 + A_G_i1{i}*qd(model.joints{i});
                    A_G1 = [A_G1 A_G_i1{i}];
                end
                I_G1 = X_G'*I_C_0*X_G;
                v_G1 = I_G1\h_G1;
            end
            
        end
        
        % 06-26-14 compute Centroidal Momentum matrix Derivative
        function [A_G, A_G_dot] = computeCMM_dot(model,qd, method)
            
            if ~exist('method','var')
                method = 1; % use recursive method by default
            end
            
            if(method == 1) % recursive
                I_C = cell(model.NB,1);
                dI_C = cell(model.NB,1); % derivative of body spatial inertias
                I_C_0 = 0;
                dI_C_0 = 0;
                
                for i = 1:model.NB
                    I_C{i} = model.I{i};
                    dI_C{i} = zeros(6,6);
                end
                
                for i = model.NB:-1:1
                    p = model.parent{i};
                    dXup{i} = model.dXup(i,qd);
                    if (p~=0)
                        I_C{p} = I_C{p} + model.i_X_pi{i}' * I_C{i} * model.i_X_pi{i};
                        dI_C{p} = dI_C{p} + (dXup{i}'*I_C{i} + model.i_X_pi{i}'*dI_C{i})*model.i_X_pi{i} + model.i_X_pi{i}'*I_C{i}*dXup{i};
                    else % floating-base
                        I_C_0 = I_C_0 + model.i_X_pi{i}' * I_C{i} * model.i_X_pi{i};
                        dI_C_0 = dI_C_0 + (dXup{i}'*I_C{i} + model.i_X_pi{i}'*dI_C{i})*model.i_X_pi{i} + model.i_X_pi{i}'*I_C{i}*dXup{i};
                    end
                end
                
                h_G1 = 0;
                % we need minus sign
                M = I_C_0(6,6);
                p_G = skew(I_C_0(1:3,4:6))/M;
                X_G = xlt(-p_G); % 0_X_G
                % compute x_G_dot from Jcom
                Jcom = model.getJcom();
                xdot_com = Jcom(4:6,:)*qd;
                %                 dX_G = xlt(-xcom_dot);
                dXtr = dXtrans(-p_G);
                dX_G = -(reshape(dXtr(:,1),[6 6])*xdot_com(1) + ...
                    reshape(dXtr(:,2),[6 6])*xdot_com(2) + ...
                    reshape(dXtr(:,3),[6 6])*xdot_com(3));
                
                i_X_G = cell(model.NB,1);
                di_X_G = cell(model.NB,1);
                A_G_i1 = cell(model.NB,1);
                dA_G_i1 = cell(model.NB,1);
                A_G1 = [];
                dA_G1 = [];
                
                for i = 1:model.NB
                    % compute pi_X_G (recursive) use the property that
                    % parent body index < children body indices
                    p = model.parent{i};
                    
                    if(p==0) % at floating base
                        i_X_G{i} = model.i_X_pi{i} * X_G;
                        di_X_G{i} = dXup{i} * X_G + model.i_X_pi{i}*dX_G;
                    else
                        i_X_G{i} = model.i_X_pi{i} * i_X_G{p};
                        di_X_G{i} = dXup{i} * i_X_G{p} + model.i_X_pi{i}*di_X_G{p};
                    end
                    
                    %                     A_G_i1{i} = i_X_G{i}'*I_C{i}*model.phi{i};
                    %                     dA_G_i1{i} = di_X_G{i}'*I_C{i}*model.phi{i} + i_X_G{i}'*dI_C{i}*model.phi{i};
                    if(i==1) % at floating base joint
                        A_G_i1{i} = i_X_G{i}'*I_C{i}*model.i_X_pi{i}*model.phi{i};
                        dA_G_i1{i} = di_X_G{i}'*I_C{i}*model.i_X_pi{i}*model.phi{i} + i_X_G{i}'*dI_C{i}*model.i_X_pi{i}*model.phi{i};
                    else
                        A_G_i1{i} = i_X_G{i}'*I_C{i}*model.phi{i};
                        dA_G_i1{i} = di_X_G{i}'*I_C{i}*model.phi{i} + i_X_G{i}'*dI_C{i}*model.phi{i};
                    end
                    %                     h_G1 = h_G1 + A_G_i1{i}*qd(model.joints{i});
                    A_G1 = [A_G1 A_G_i1{i}];
                    dA_G1 = [dA_G1 dA_G_i1{i}];
                    %TODO: fix the size of A, dA
                end
                %                 I_G1 = X_G'*I_C_0*X_G;
                %                 v_G1 = I_G1\h_G1;
                A_G = A_G1;
                A_G_dot = dA_G1;
            end
            
            
            
        end
        
        
        %% 06-23-14 translate a model to featherstone's model
        % this function will allow us to translate our model to
        % featherstone's model definition so that we can use featherstone's
        % library + drakes's library for computations
        function model = model2FeatherStoneV2(model)
            m = [];
            m.NB = model.NB;
            m.gravity = model.ag_0(4:end);
            
            for i=1:model.NB
                m.parent(i) = model.parent{i};
                m.Xtree{i} = model.i_XL_pi{i};
                %                 m.Xup{i} = model.i_X_pi{i};
                %                 m.XJ{i} = model.X_J{i};
                if(i==1)
                    m.jtype{i} = 'R';
                else
                    jointaxis = find(model.phi{i}(1:3) == 1);
                    if(jointaxis == 1) % Rx
                        m.jtype{i} = 'Rx';
                    elseif(jointaxis == 2) % Rx
                        m.jtype{i} = 'Ry';
                    elseif(jointaxis == 3) % Rx
                        m.jtype{i} = 'Rz';
                    end
                end
                m.I{i} = model.I{i}; % use mcI is an option
            end
            
            % if the base is floating-base
            if(model.joints{1}(end) == 6)
                m = floatbaseV2(m);
            end
            
            %             m = AddAppearance(m)
            
            model.featherstone = m;
        end
        
        %         % add apperance to the featherstone (v2) model
        %         function m = AddAppearance(m)
        %
        %             % floor
        %             % tiles profile: ( X, Y, Z, tilesize, colour )
        %             m.appearance.base = { 'tiles', [-0.3 0.3; -0.3 0.3; 0 0], 0.1 };
        %             %             { 'tiles', [-0.3 0.3; -0.3 0.3; 0 0], 0.1 };
        %             %               { 'colour', [0.9 0 0], 'line', [0 0 0; 2 0 0], ...
        %             %     'colour', [0 0.9 0], 'line', [0 0 0; 0 2 0], ...
        %             %     'colour', [0 0.3 0.9], 'line', [0 0 0; 0 0 2] };
        %             % body
        %             for i=6:m.NB
        %                 %                 NB = i - 5;
        %                 m.appearance.body{i-5} = { 'box', [0 0 0; 1 1 1] };
        %                 %                 { 'box', [0 0 0; 3 2 1] }
        %                 % box: {coords,colour) creates a box having the specified parent,
        %                 % coordinates and colour, and returns its handle.  Coords is a 2x3 matrix
        %                 % containing  any two diametrically opposite vertices of the box.  The box
        %                 % itself is aligned with the coordinate axes.  Colour is an RGB colour
        %                 % vector.
        %                 %                 { 'facets', npt, 'cyl', [0 -T/2 0; 0 T/2 0], R };
        %             end
        %
        %             % Camera settings
        %
        %             m.camera.body = 1;
        %             m.camera.direction = [1 0 0.1];
        %             m.camera.locus = [0 0.5];
        % %             m.featherstone = m;
        %         end
    end
    
    
end



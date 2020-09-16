classdef TO_main < handle
    
    properties
        FT = 1;
        rm = 0.04;
        pen = 3;      
        
    end
    
    methods
        function this = TO_main()
        end
        
        %% PRE---------------------------------------------------
        function C = PRE(this,res,input,num,fn)
            filename = fn + num; 
            volfrac = input.volfrac;
            len = res(2)-1; 
            wid = res(1)-1; 
            
            % Defining the cases: CS 1 for passive, 0 for active, shp 1 for
            % circle 2 for square-----------------------------
%             shp = input.shp;
            % Define the passive and active
            % elements---------------------------
%             inp = zeros(wid,len);
%             active = input.active;
%             passive= input.passive;
              elem = input.element;
%             x = input.xpos;
%             y = input.ypos;
%             radius = input.radius;
%             cs = input.cs;
%             if cs == 1
%                 passive = element(x,y,radius,inp,1,shp);
%             elseif cs ==2
%                 active = element(x,y,radius,inp,0,shp);
%             end
            
            % Loading Condition-----------------------------------------
            force_pos = input.force_pos;
            force_dir = input.force_dir;
            %force_pos = 2*(wid+1)*(len+1); % For multiple use = [N1,N2];
            %force_dir = -1; %% +1 or -1 that's it: THERE IS NO EFFECT of force on shape : Multiple Case, use array = [-1 -1]
            if size(force_dir,2) ==1
                F = sparse(force_pos,1,force_dir, 2*(wid+1)*(len+1),1);% (Point number,1,direction,max number of points,1)
            elseif size(force_dir,2) ==2
                F = sparse(force_pos,[1 2],force_dir,2*(wid+1)*(len+1),2);
            else 
                F = sparse(force_pos,[1 2,3,4],force_dir,2*(wid+1)*(len+1),4);
            end
            
            fixeddofs = input.fixeddofs;% Boundary Condition----------------------

            % Print the initial Conditions------------------------------
            initial_model(volfrac,elem,force_pos,force_dir,fixeddofs,len,wid,filename);

            % TOPOP OPTIMIZATION-------------------------------
            [topopt,C] = this.RUN(len,wid,volfrac,this.pen,this.rm*wid,this.FT,F,fixeddofs,elem,size(force_dir,2));
            final_output = ones(wid+1,len+1);
            final_output(2:wid+1,2:len+1) = topopt;
            imwrite(final_output,"./data/OUT/" + "OUT_"+ filename + ".png")
        end
        
%% -----------------------------------------MAIN ALGORITHM: SIMP with Helmotz Filter -----------------------------------------------------------------
        
        function [topopt,c] = RUN(this,nelx,nely,volfrac,penal,rmin,ft,F,fixeddofs,elem,no)
            
            % MATERIAL PROPERTIES------------------------------------
            E0 = 1;
            Emin = 1e-9;
            nu = 0.3;
            
            % PREPARE FINITE ELEMENT ANALYSIS------------------------------
            A11 = [12  3 -6 -3;  3 12  3  0; -6  3 12 -3; -3  0 -3 12];
            A12 = [-6 -3  0  3; -3 -6 -3 -6;  0 -3 -6  3;  3 -6  3 -6];
            B11 = [-4  3 -2  9;  3 -4 -9  4; -2 -9 -4 -3;  9  4 -3 -4];
            B12 = [ 2 -3  4 -9; -3  2  9 -2;  4  9  2  3; -9 -2  3  2];
            KE = 1/(1-nu^2)/24*([A11 A12;A12' A11]+nu*[B11 B12;B12' B11]);
            nodenrs = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);
            edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1,nelx*nely,1);
            edofMat = repmat(edofVec,1,8)+repmat([0 1 2*nely+[2 3 0 1] -2 -1],nelx*nely,1);
            iK = reshape(kron(edofMat,ones(8,1))',64*nelx*nely,1);
            jK = reshape(kron(edofMat,ones(1,8))',64*nelx*nely,1);
            
            % DEFINE LOADS AND
            % SUPPORTS--------------------------------------
            U = zeros(2*(nely+1)*(nelx+1),2);
            alldofs = [1:2*(nely+1)*(nelx+1)];
            freedofs = setdiff(alldofs,fixeddofs);
            
            % PREPARE FILTER -----------------------------------
            iH = ones(nelx*nely*(2*(ceil(rmin)-1)+1)^2,1);
            jH = ones(size(iH));
            sH = zeros(size(iH));
            k = 0;
            for i1 = 1:nelx
              for j1 = 1:nely
                e1 = (i1-1)*nely+j1;
                for i2 = max(i1-(ceil(rmin)-1),1):min(i1+(ceil(rmin)-1),nelx)
                  for j2 = max(j1-(ceil(rmin)-1),1):min(j1+(ceil(rmin)-1),nely)
                    e2 = (i2-1)*nely+j2;
                    k = k+1;
                    iH(k) = e1;
                    jH(k) = e2;
                    sH(k) = max(0,rmin-sqrt((i1-i2)^2+(j1-j2)^2));
                  end
                end
              end
            end
            H = sparse(iH,jH,sH);
            Hs = sum(H,2);
            
            % INITIALIZE ITERATION
            % -------------------------------------------
            x = repmat(volfrac,nely,nelx);
            xPhys = x;
            loop = 0;
            change = 1;

            % START ITERATION-------------------------------
            while change > 0.02
              loop = loop + 1;
              % BREAKING THE LOOP
              if loop>400
                  break;
              end

             % FE-ANALYSIS -------------------------------------
              sK = reshape(KE(:)*(Emin+xPhys(:)'.^penal*(E0-Emin)),64*nelx*nely,1);
              K = sparse(iK,jK,sK); K = (K+K')/2;
              if no == 1
                U(freedofs) = K(freedofs,freedofs)\F(freedofs);
              elseif no ==2
                U(freedofs,:) = K(freedofs,freedofs)\F(freedofs,:);
               elseif no ==4
                U(freedofs,:) = K(freedofs,freedofs)\F(freedofs,:);
              end
              
             % OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
             % ---------------------------
              %ce = reshape(sum((U(edofMat)*KE).*U(edofMat),2),nely,nelx);
              %c = sum(sum((Emin+xPhys.^penal*(E0-Emin)).*ce));
             % dc = -penal*(E0-Emin)*xPhys.^(penal-1).*ce;
                dc=0;
                c = 0;
                for i = 1:size(F,2)
                    Ui = U(:,i);
                    ce = reshape(sum((Ui(edofMat)*KE).*Ui(edofMat),2),nely,nelx);
                    c = c + sum(sum((Emin+xPhys.^penal*(E0-Emin)).*ce));
                    dc = dc - penal*(E0-Emin)*xPhys.^(penal-1).*ce;
                end
              dv = ones(nely,nelx);

              % FILTERING/MODIFICATION OF SENSITIVITIES
              % -----------------------------
              if ft == 1
                dc(:) = H*(x(:).*dc(:))./Hs./max(1e-3,x(:));
              elseif ft == 2
                dc(:) = H*(dc(:)./Hs);
                dv(:) = H*(dv(:)./Hs);
              end
              
              %  OPTIMALITY CRITERIA UPDATE OF DESIGN VARIABLES AND
              %  PHYSICAL DENSITIES----------------------------
              l1 = 0; l2 = 1e9; move = 0.2; poop = 1;
              while (l2-l1)/(l1+l2) > 0.02
                if poop>300
                    break;
                end
                lmid = 0.5*(l2+l1);
                xnew = max(0,max(x-move,min(1,min(x+move,x.*sqrt(-dc./dv/lmid)))));
                if ft == 1
                  xPhys = xnew;
                elseif ft == 2
                  xPhys(:) = (H*xnew(:))./Hs;
                end
                xPhys(elem ==1) = 0;
                xPhys(elem ==2) = 1;
                if sum(xPhys(:)) > volfrac*nelx*nely, l1 = lmid; else l2 = lmid; end
                poop = poop + 1;
              end
              change = max(abs(xnew(:)-x(:)));
              x = xnew;
              
              %PRINT RESULTS --------------------
              if rem(loop,10)==0
                fprintf('-');
              end
              %fprintf(' It.:%5i Obj.:%11.4f Vol.:%7.3f ch.:%7.3f\n',loop,c,mean(xPhys(:)),change);
              % PLOT DENSITIES
              topopt = 1- xPhys;
              %colormap(gray); 
              %imshow(1-xPhys); caxis([0 1]); axis equal; axis off; drawnow;
            end  
            
        end
        %% ----------GET COMP-----------------
        
        function C = GET_COMP(input)
            
            volfrac = input.volfrac;
            len = res(2)-1; 
            wid = res(1)-1; 
            elem = input.element;
            force_pos = input.force_pos;
            force_dir = input.force_dir;
            %force_pos = 2*(wid+1)*(len+1); % For multiple use = [N1,N2];
            %force_dir = -1; %% +1 or -1 that's it: THERE IS NO EFFECT of force on shape : Multiple Case, use array = [-1 -1]
            if size(force_dir,2) ==1
                F = sparse(force_pos,1,force_dir, 2*(wid+1)*(len+1),1);% (Point number,1,direction,max number of points,1)
            elseif size(force_dir,2) ==2
                F = sparse(force_pos,[1 2],force_dir,2*(wid+1)*(len+1),2);
            else 
                F = sparse(force_pos,[1 2,3,4],force_dir,2*(wid+1)*(len+1),4);
            end
            
            fixeddofs = input.fixeddofs;% Boundary Condition
        end

            
    end
    
end
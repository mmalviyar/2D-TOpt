classdef TO < handle
    
    properties
        VF = [0.25,0.3,0.35,0.4,0.5];% Volume fraction
        BC = {};
        BCHD = {};% Boundary cases
        FC = {}
        FCHD = {}% Force cases
        Fmagx =  round(cos(0:pi/6:pi/2),4);
        Fmagy =  round(sin(0:pi/6:pi/2),4);
        Fdir = [-1 1];
        Fnodes = NaN; 
        res = [];
        HD = [];
        pt = 5; 
        lr = [1.5,2,2.5,3];
        tolerance = 5;
        
    end
    
    methods
        
        function this = TO()
                this.res = [20,60]; % wid,len
                this.HD = this.pt*(this.res);                
                this.init(this.res,1);
                this.init(this.HD,round(this.pt));
        end        
       
        function init(this,res,p)
            
                len = res(2)-1; wid = res(1)-1;
                xdof = 1:2:2*(len+1)*(wid+1);
                xdof = reshape(xdof,wid+1,len+1);
                ydof = xdof + 1;

                % CASE 1: CANTILEVER
                BC_1 = [xdof(:,1);ydof(:,1)]; % CANTILEVER
                FC_1 = [xdof(1:p:end,end), ydof(1:p:end,end)]; % CANTILEVER FX FY END
                
                % CASE 2: MBB with ROLLER SUPPORTS                               
                BC_2 = union(xdof(:,1),ydof(end,end));
                FC_2 = [xdof(1,1:p:end/2)',ydof(1,1:p:end/2)']; % LEFT TOP TO MID
                              
                                
                % CASE 3: %MBB with fixed on one end and roller on the other
                BC_3 =  [xdof(end,1),ydof(end,1),ydof(end,end)];
                FC_3 = [xdof(1,end/4 + 1 :p: 3*end/4)',ydof(1,end/4 + 1 :p: 3*end/4)']; % TOP MID FY
                                
                
                % CASE 4: MBB both end fixed, Force on Bottom
                BC_4 = [xdof(end,1),ydof(end,1),ydof(end,end), xdof(end,end)]; 
                FC_4 = [xdof(end,end/4 + 1 :p: 3*end/4)',ydof(end,end/4 + 1 :p: 3*end/4)']; % BOT MID FY
                
                
                % CASE 5: MBB both end fixed, Force on TOP
                BC_5 = [xdof(end,1),ydof(end,1),ydof(end,end), xdof(end,end)]; 
                FC_5 = [xdof(1,end/4 + 1 :p: 3*end/4)',ydof(1,end/4 + 1 :p: 3*end/4)']; % TOP MID FY

                              
                % SAVING THE FORCES AND BC NODES
                if p ==1
                    this.BC = {BC_1,BC_2,BC_3,BC_4,BC_5};
                    this.FC = {FC_1, FC_2, FC_3,FC_4,FC_5};
                else
                    this.BCHD = {BC_1,BC_2,BC_3,BC_4,BC_5};
                    this.FCHD = {FC_1, FC_2, FC_3,FC_4,FC_5};
                end              
                
                this.Fnodes = length(FC_1);
           
        end
        
               
        function [gncases,gncasesHD] =  gencases(this,total) % GENERATE CASES
           
           num_max = total;
           DATA = zeros(num_max,10);
           gncases = {};   
           gncasesHD = {};   
           
           for num = 1: num_max
              
              % SELECTING BOUNDARY CONDITION--------
              B = randi([1,length(this.BC)],1,1); % BC
              input.fixeddofs = this.BC{:,B};
              inputHD.fixeddofs = this.BCHD{:,B};
              
              % SELECTING FORCE LOCATIONS---------------
              Floc = randi([1,this.Fnodes],1,1); % LOAD NODE 
              input.force_pos = [this.FC{B}(Floc,1),this.FC{B}(Floc,2)];
              inputHD.force_pos = [this.FCHD{B}(Floc,1),this.FCHD{B}(Floc,2)];                

              % SELECTING FORCE MAG AND DIRECTION---------------
              Fd = this.Fdir(randi([1,2],1,1)); % LOAD DIRECTION
              Fm = randi([1,length(this.Fmagx)],1,1); % LOAD MAG 
              input.force_dir = Fd'.*[this.Fmagx(Fm),this.Fmagy(Fm)];  % LOAD MAG*Direction
              inputHD.force_dir = Fd'.*[this.Fmagx(Fm),this.Fmagy(Fm)];
               
              % VOLUME FRACTION--------------
              vf = this.VF(randi([1,length(this.VF)],1,1)); % VOLUME FRAC
              inputHD.volfrac =  vf;
              input.volfrac = vf;
              
              
              % ELEMENT CASE-------------
              el_case = randi([1,3],1,1); % Element Generator: 1: NO, 2: ACTIVE ,3: PASSIVE, 4: AP        
              [elem,elemHD,xloc,yloc,radii, ratio] = this.element_generator(el_case,input.fixeddofs,vf,input.force_pos);
              input.element = elem;
              inputHD.element = elemHD;
                           
              % FINAL STRUCTURE--------------------------
              gncases{num} = input; 
              gncasesHD{num} = inputHD; 
              % DATA FILE AS STORAGE
              DATA(num,:) = [vf,B,Floc,Fd,Fm,el_case,xloc,yloc,radii,ratio];   
              
           end
           
           %Uniqueness------------
%            [data,K] = unique(DATA,'rows');
%            data = data(1:total,:);
%            gncases = gncases{K};
%            gncasesHD = gncasesHD{K};
%            P = 1:total;
%            gncases = gncases{P};
%            gncasesHD = gncasesHD{P};           
           
           %
           save('training_data.mat','gncases','gncasesHD','DATA','-v7.3');
           
        end
        
        function [elem,elemHD,xloc,yloc,radii,a] = element_generator(this,condition,Bnodes,VF,fnodes)
            
            % INTIALIZING
            len = this.res(2)-1;
            wid = this.res(1)-1;            
            inp = zeros(wid,len);
            
            widHD = this.HD(1)-1;
            lenHD = this.HD(2)-1;
            inpHD = zeros(widHD,lenHD);            

            
            % LOCATIONS RESTRICTED BY BOUNDARY CONDITIONS --ADATABLE
            loc = zeros(wid+1,len+1);
            xdof = 1:2:2*(len+1)*(wid+1);
            xdof = reshape(xdof,wid+1,len+1);
            ydof = xdof + 1;
            BCx = ismember(xdof,Bnodes);
            BCy = ismember(ydof,Bnodes);
            loc(BCx==1) = 1;
            loc(BCy==1) = 1;
            Fx = ismember(xdof,fnodes);
            Fy = ismember(ydof,fnodes);
            loc(Fx==1) = 1;
            loc(Fy==1) = 1;
            tol = this.tolerance; % Tolerance
            
            switch condition
                
                                
                case 1 % NOTHING 
                    elem = inp;
                    elemHD = inpHD;
                    xloc = 0; yloc = 0; radii = 0; a = 0;
                    
                    
                case 2 % ACTIVE ELEMENTS
                    
                    % GENERATING RADII ARRAY
                    act_maxradii = sqrt(0.5*VF*(len*wid));
                    active_radius = round(linspace(len/10,act_maxradii,3)); 
                    
                    % SELECTING ONE RADII
                    R = randi([1,length(active_radius)],1,1);
                    radii = active_radius(R);
                    
                    % RECTANGLE DIMENSIONS            
                    a = this.lr(randi([1,length(this.lr)],1,1)); 
                    w = radii/sqrt(a); l = a*w;
                    
                    % GENERATING ELEMENT                    
                    xloc = randi([1,len],1,1);
                    yloc = randi([1,wid],1,1);
                    [X,Y] = this.high_Node(xloc,yloc);
                    elem = this.square_generator(inp,0,xloc,yloc,w,l); 
                    elemHD = this.square_generator(inpHD,0,X,Y,this.pt*w,this.pt*l);
                    
                case 3 % JUST PASSIVE
                    
                    while 1
                        
                       % GENERATE ARRAY OF RADIUS           
                        pas_maxradii = sqrt((1-2*VF)*(len*wid)); % Add /pi for circle
                        passive_radius = round(linspace(len/8,pas_maxradii,3));

                        % SELECT RADIUS
                        R = randi([1,length(passive_radius)],1,1);
                        radii = passive_radius(R);
                        a = this.lr(randi([1,length(this.lr)],1,1));
                        w = radii/sqrt(a); l = a*w;                        
                        
                        locb = imgaussfilt(loc,0.5,'FilterSize',11);
                        %locb(:,1:tol) = 1;
                        %locb(:,end-tol:end) = 1; 
                        
                        A = find(locb==0); 
                        K = A(randi([1,length(A)],1,1));                    
                        xloc =  ceil(K/this.res(1));                    
                        yloc = rem(K,this.res(1));
                        
                        
                        [X,Y] = this.high_Node(xloc,yloc);
                        if yloc<round(wid-w/2-tol+1) || w/2+tol<yloc
                            if 1<Y && Y<round(widHD-this.pt*w/2-tol*this.pt) || this.pt*w/2+tol*this.pt<Y
                                
                                elem = this.square_generator(inp,1,xloc,yloc,w,l); 
                                elemHD = this.square_generator(inpHD,1,X,Y,this.pt*w,this.pt*l);  
                                break;
                            end
                        end
                        
                    end
                    
            end
                    
        end
                
        
        
        function control(this,number)
            
                 T = TO_main();
                 if length(number) == 1
                    [train_cases,train_casesHD] = this.gencases(no_of_cases);
                    a = 1; b = number;
                 elseif length(number) == 2
                    load('training_data.mat');
                    train_cases = gncases;
                    train_casesHD = gncasesHD;
                    a = number(1); b = number(2);
                 end
                 this.gen_directory();
                 
                 comp_time = zeros(b-a+1,5);
                 
                 for i = a:b
                     try 
                         fprintf("CASE NUMBER: " + i);
                         inp = train_cases{i};
                         inpHD = train_casesHD{i};
                         
                         time_start = now();
                         fprintf("|L-");
                         CL = T.PRE(this.res,inp,i,"LOW");
                         tl = now()-time_start;
                         time_start = now();
                         fprintf("|H-");                     
                         CH = T.PRE(this.HD,inpHD,i,"HIGH");
                         th = now()-time_start;
                         comp_time(i,:) = [i CL 3600*24*tl CH 3600*24*th];
                         fprintf("\n");  
                     catch
                         print("ERROR")     
                     end
                                       
                     
                 end
                 
                 writematrix(comp_time,"comp_time_" + a +"_" + b + "_.csv");                
                 
                 

        end
        
        function mat = square_generator(this,input,cs,x,y,w,l)
                mat = input;                
                                
                for i = 1:size(input,2)
                    for j = 1:size(input,1)
                        if abs(j-y) <= round(w/2) && abs(i-x) <= round(l/2)
                            if cs == 0
                                mat(j,i) = 2; % Active
                            else
                                mat(j,i) = 1; % Passive
                            end
                        end
                    end
                end
            
        end
        
        function [xhd,yhd] = high_Node(this,x,y)
            
            xhd = this.pt*(x-1) + (this.pt+1)/2;
            yhd = this.pt*(y-1) + (this.pt+1)/2;
            
        end

        function gen_directory(this)
             mkdir data
             mkdir data/VF
             mkdir data/F
             mkdir data/B
             mkdir data/OUT
        end
        function S = round_odd(this,S)
            % round to nearest odd integer.
            idx = mod(S,2)<1;
            S = floor(S);
            S(idx) = S(idx)+1;
        end
        
        
        
    end
end

    
  
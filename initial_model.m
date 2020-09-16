function initial_model(vf,elem,force_pos,force_dir,fix_pos,len,wid,filename)

%% Plotting volume
% volume fraction as base -- SPITS OUT GRAY IMAGE WITH INTENSITY == GRAYSCALE
ch_vf = vf*ones(wid+1,len+1);
ch_bx = ones(wid+1,len+1);
ch_by = ones(wid+1,len+1);
ch_fx = 0.5*ones(wid+1,len+1);
ch_fy = 0.5*ones(wid+1,len+1);
fmax = 1;
%imwrite(ch_vf,"VF_" + filename + ".png")
%% Plotting active and passive elements
    % GRAY  SCALE IMAGE WITH ACTIVE AND PASSIVE ELEMENTS
    elem2 = zeros(wid+1,len+1);
    elem2(2:wid+1,2:len+1) = elem;
%     passive2 = zeros(wid+1,len+1);
%     passive2(2:wid+1,2:len+1) = passive;
    if norm(elem2) > 0
        ch_vf(elem2 == 1) = 1; % passive
        ch_vf(elem2 == 2) = 0; % active
    end
%imshow(ch_vf)   
imwrite(ch_vf,"./data/VF/"+ "VF_" + filename + ".png")
%% DOF To PIXEL REPRESENTATION
    maxdof = 2*(len+1)*(wid+1);
    xdof = 1:2:maxdof-1;
    xdof = reshape(xdof,[wid + 1,len + 1]);
    ydof = 2:2:maxdof;
    ydof = reshape(ydof,[wid + 1,len + 1]);
    
%% Plotting BC       
    BCx = ismember(xdof,fix_pos);
    BCy = ismember(ydof,fix_pos);
    ch_bx(BCx == 1) = 0; %y
    ch_by(BCy == 1) = 0; %x
imwrite(ch_bx,"./data/B/"+ "BFX_" + filename + ".png")
imwrite(ch_by,"./data/B/"+ "BFY_" + filename + ".png")
%% Ploting Force 
       for i = 1:size(force_dir,2)
           Px = ismember(xdof,force_pos(i));
           Py = ismember(ydof,force_pos(i)); 
           Val = 0.5*(1 + force_dir(i)/fmax);
           ch_fy(Py==1) = Val; ch_fx(Px==1) = Val;  
       end  

imwrite(ch_fx,"./data/F/"+ "FX_" + filename + ".png")
imwrite(ch_fy,"./data/F/"+ "FY_" + filename + ".png")
end
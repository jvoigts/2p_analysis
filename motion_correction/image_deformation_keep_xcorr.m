%
% Demo of a method for deforming and re-sampling 
%
%
%
% Mar 29, 2015
% Jakob Voigts (jvoigts@mit.edu)


Imsize=[80 80];
Nframes=5000;
Nregions=2; % vertically stacked regions to use for deformations - just do 2 for this test

% assume that the deformation is known
% so we'll just deform noise and check the spatial cross correlation

I=abs(rand(Imsize(1),Imsize(2),Nframes)).*1; 

if 0 % add grid to check re-sampling?
    for j=1:Nframes
        I(1:20:end,:)=I(1:20:end,:)+.1;
        I(:,1:20:end)=I(:,1:20:end)+.1;
    end;
end;

translations=randn(Nregions,2,Nframes).*25;

%%

[x y]=meshgrid([1:Imsize(1)],[1:Imsize(2)]);

regionboundaries=round(linspace(1,Imsize(1),Nregions+1));
hregheight = round(mean(diff(regionboundaries))/2);


for cond=1:3 % condition 1: bilinear 2: nearest neighbor(NN) 3: NN with replacement
    
    last_Iout=I(:,:,1); % keep copy of last aligned frame for doubled pixel replacement, init with 1st input frame
    disp(cond);
    
    for j=1:Nframes
        
        Iin=I(:,:,j);
        
        Iout=zeros(size(Iin));
        
        %deform input image based on translations
        Dx=Iin.*0;
        Dy=Iin.*0;
        
        % fill in deformation grid, then onterpolate to give pixel-wide
        % displacement for entire image
        sample_y = [ regionboundaries(1:end-1)'+hregheight  regionboundaries(1:end-1)'+hregheight ];
        sample_y(1,:)=1;
        sample_y(end,:)=size(Iin,1);
        sample_x= [ones(1,Nregions)' ones(1,Nregions)'.*size(Iin,2)];
        
        Vx= [translations(:,1,j) translations(:,1,j)]; %translation values (these would be estimated for each region)
        Vy= [translations(:,2,j) translations(:,2,j)];
        
        Dx = interp2(sample_x,sample_y,Vx,x,y,'spline'); % calculate displacement per pixel
        Dy = interp2(sample_x,sample_y,Vy,x,y,'spline');
        
        yd=x+Dx; % calc complete pixel coordinates
        xd=y+Dy;
        
        m=1-(  yd>size(Iin,2) | xd>size(Iin,1) |  xd<1 | yd<1 ); % which pixels dont have data assigned?
        
        if cond==1
            % interpolate
            xd_l=min(size(Iin,2),max(1,floor(xd))); % round and clip
            yd_l=min(size(Iin,1),max(1,floor(yd)));
            xd_h=min(size(Iin,2),max(1,ceil(xd))); % round and clip
            yd_h=min(size(Iin,1),max(1,ceil(yd)));
            
            %  get sub-pixel position component
            py=xd-floor(xd);
            px=yd-floor(yd);
            
            ii=sub2ind(size(Iin),yd_l,xd_l); % do 4x NN once for each neighboring pixel
            I_ll=Iin(ii);
            ii=sub2ind(size(Iin),yd_h,xd_l);
            I_hl=Iin(ii);
            ii=sub2ind(size(Iin),yd_l,xd_h);
            I_lh=Iin(ii);
            ii=sub2ind(size(Iin),yd_h,xd_h);
            I_hh=Iin(ii);
            
            Iout= ((1-px).*I_ll + (1-px).*I_lh + px.*I_hl + px.*I_hh +... % bilinear interpolation
                (1-py).*I_ll  + (1-py).*I_hl + py.*I_lh + py.*I_hh)./4;
           
            Iout=Iout.*m; % zero extrapolated data
        else
            % nearest neigbor:
            xdnn=min(size(Iin,1),max(1,round(xd))); % round and clip for NN
            ydnn=min(size(Iin,2),max(1,round(yd)));
            ii_nn=sub2ind(size(Iin),xdnn,ydnn);
            
            Iout=Iin(ii_nn);
            Iout=Iout.*m; % zero extrapolated data
            
            if cond==3
                
                %start by identifying all pixels that have fan-out > 1
                [C,IA,IC] = unique(ii_nn,'stable'); % get all unique ones (first occurrences in IA)
                fanouts=1:numel(Iout); fanouts(IA)=[]; % this leaves only the 2nd and further copies of pixels from the original to the deformed image
                
                if 0 % display which ones were doubled?
                    disp_fanouts=Iin.*0;
                    disp_fanouts(fanouts)=1;
                    imagesc(Iout+disp_fanouts);
                end;
                Iout(fanouts)=last_Iout(fanouts); % set doubles to independent samples from prev. image to maintain spatial xcorr
                
                last_Iout=Iout;
            end;
            
        end;
        
        registeredStack{cond}(:,:,j)=Iout;
        
        if 0
            clf;
            subplot(121); imagesc(Iin); daspect([1 1 1]);
            subplot(122); imagesc(Iout); daspect([1 1 1]);
            drawnow
        end;
        
    end;
    
    %calculate xcorrs
    px=50;
    py=50;
    
    mm=[-8:8];
    ref= (squeeze( registeredStack{cond}(ceil(py),ceil(px),:) ));
    xc{cond}=Iin.*0;
    for i=mm+px % first subsample every 2nd
        
        for j=mm+py
            c=corrcoef(squeeze(registeredStack{cond}(i,j,:)),ref);
            xc{cond}(i,j)=c(1,2);
        end
    end;
end;

clf;
imagesc([xc{1}(mm+px,mm+py) xc{2}(mm+px,mm+py) xc{3}(mm+px,mm+py)])



clear all; close all
warning off
showLiveTracking = 0; % if you want to see live tracking results
showCurvAxis = 0; % to display curvature axis of larva
   
h = fspecial('average',20);

areaThresh = 400; % fly larva area in pixels, maybe try 300-600 range, increasing might help more (default is 400 for fly larva)
weightMetric = 0.5; % value from 0 to 0.9999, higher will bias the sum pixel intensity, lower will bias the overlap
smallObjThresh = 25; % removes small objects, probably won't change much
% if you don't see anything it doesn't mean it's not working necessarily, it could
% means that the fish isn't moving. Make the 'filen' number of the loop
% equal to the folder you want to check

H = fspecial('disk',5);
filtParam = 1; % 1-3; 3 is more smoothed. for movies with less contrast, a value closer to zero (i.e. 0.01) works better
lastArea = 1;
frameSkip = 10; % use smaller number to get better background estimate (default is 30 for 900 frame movie)

foldern{1} = 'I:\plexb and sema2b single'; % thresh = 8



%% 

for fileNo = 1:numel(foldern)
    fileNo
    clear BG
    filen = foldern{fileNo};
    cd(filen);
    aviFiles = dir('*.avi');
    
     for aviNum = 1:length(aviFiles)
%     for aviNum = 7

    thisAviFile = (aviFiles(aviNum).name);
    thisMaskName = thisAviFile;
    thisMaskName(strfind(thisMaskName,'.'):end)=[];
    thisMaskName = [thisMaskName '_mask.mat'];
    load(thisMaskName);

    obj = VideoReader(thisAviFile);
    nframes = obj.FrameRate*obj.Duration;

    nWells = double(max(mask(:)));
    % make movie into matrix
    bgf = [1:frameSkip:nframes];
    BG=uint8(zeros(obj.Height,obj.Width,length(bgf)));
    for ff=1:length(bgf)
        temp  = read(obj,bgf(ff));
        BG(:,:,ff) = temp(:,:,1);
%         hb(:,:,ff) = im2single(BG(:,:,ff)); % need to convert to single precision to use 'mode' function
    end
    
    for wn = 1:nWells
%     for wn = 14;
        wellmask = mask==wn;
        [xvals,yvals]=find(wellmask);
        if ~isempty(xvals)
        bottom(wn) = max(xvals);
        rightside(wn) = max(yvals);
        top(wn) = min(xvals);
        leftside(wn) = min(yvals);
        else
            bottom(wn) = 0;
            rightside(wn)=0;
            top(wn)=0;
            leftside(wn) =0;
        end
    end

       %% mask edits
       clear wellBackgrounds
% for wellNum = [40]
    for wellNum = 1:nWells

    wellmask = mask==wellNum;

    clear bgDiffVals
    for bgFrameNo = 1:size(BG,3)-1
        bg1 = BG(:,:,bgFrameNo);
        bg2 = BG(:,:,bgFrameNo+1);
        bg1(~wellmask)=0;
        bg2(~wellmask)=0;
        bgDiff = imabsdiff(bg1,bg2);
        bgDiffVals(bgFrameNo,1)=sum(bgDiff(:));
    end

    [~,bgDiffSort] = sort(bgDiffVals,'descend');
    
    frames2use = bgDiffSort(1:round(numel(bgDiffSort)*.2));

    clear bgWellMov
    frameCnt = 1;
    for ii=1:numel(frames2use)
        bgWellMov(:,:,frameCnt) =  im2single(imcrop(BG(:,:,frames2use(ii)),[leftside(wellNum) top(wellNum) rightside(wellNum)-leftside(wellNum) abs(top(wellNum)-bottom(wellNum))]));
        frameCnt = frameCnt+1;
    end
    
%     wellBgMode = mode(bgWellMov,3);
%     wellBgMax = max(bgWellMov,[],3);
%     wellBg = (wellBgMode+wellBgMax)./2;
    wellBg = mode(bgWellMov,3);
%     imshow(wellBgMode);figure;imshow(wellBgMax);figure;imshow(wellBg)
    
    wellBg = imgaussfilt(wellBg,filtParam);

    wellBackgrounds{wellNum,1}=im2uint8(wellBg);
        
    end
    save(thisMaskName,'wellBackgrounds','-append');
    
%% tracking

 
    lastArea = 0; 
    clear cents areas curveMov
    cents = cell(nframes,nWells);
    areas = cell(nframes,nWells);
    circularity = cell(nframes,nWells);
    meanCurvatures = zeros(nframes,nWells);
    maxCurvatures = zeros(nframes,nWells);
    sumCurvatures = zeros(nframes,nWells);
    gofFits= zeros(nframes,nWells);
    
    threshVals = cell(nframes,nWells);

    cents(:,:)={[0,0]};

    tic
    for fn = 1:nframes
%     for fn = 1:6000
    disp(['Movie number: ' num2str(aviNum) ', Frame number: ' num2str(fn)])
        fn
        %tic
        f=read(obj,fn);
        f=rgb2gray(f);
%         f = f(:, :, 1);
%         fFilter = medfilt2(f,[filtParam filtParam]);  % get rid of RGB dimension (convert to grayscale)
%         fFilter = imgaussfilt(f,filtParam);  % get rid of RGB dimension (convert to grayscale)
% 
        for wellNum = 1:nWells
%           for wellNum = [1]
  
            try
                wellmask = mask==wellNum;
                [iCoord,jCoord] = ind2sub(size(wellmask),find(wellmask));

                fFilterCrop = imgaussfilt(f,filtParam);
                
%                 fFilterCrop = f;

                fFilterCrop=imcrop(fFilterCrop,[min(jCoord) min(iCoord) max(jCoord)-min(jCoord) max(iCoord)-min(iCoord)]);

    %             bgFilt = medfilt2(wellBackgrounds{wellNum,1},[filtParam filtParam]);
    %             frameDiff=mat2gray(imabsdiff(fFilterCrop,bgFilt));
                bgImg = wellBackgrounds{wellNum,1};
                frameDiff=mat2gray(imabsdiff(fFilterCrop,bgImg));
                
%                 imshow(bgImg)
%                 figure;
%                 imshow(fFilterCrop)   
%                 
%                 imshow((frameDiff))
                
%                 frameDiffFilt = imgaussfilt(frameDiff,0.01);
%                 imshow(frameDiffFilt)
                
                fishArea=1;
                thresh = 1;
                threshStep = 0.01;
                while fishArea <= areaThresh 
                    thresh=thresh-threshStep;
%                     thresh
                    bw = (frameDiff >= thresh);
                    bw2 = bwareaopen(bw, smallObjThresh, 4); % merge pixels into area if connectivity is 8 - this keeps areas of at least 200px and connectivity of 8
                    bw2props = regionprops(bw2,'Area');
                    fishArea=max([bw2props(:).Area]);
                    if isempty(fishArea)
                        fishArea=0;
                    end
%                     imshow(bw2)    
%                     threshVals{fn,wellNum}=thresh;
                end
          
                thresh = thresh + threshStep;
                
                bw = (frameDiff >= thresh);
                bw2 = bwareaopen(bw, smallObjThresh, 4); % merge pixels
                
                L = bwlabel(bw2); % label each region that passed the bwareaopen critera
                s = regionprops(L,frameDiff,'PixelIdx','area','centroid','PixelIdxList','WeightedCentroid','Solidity','Circularity'); % calculate the area, centroid, and axis length of each region
                area_vector = [s.Area]; 

                clear regionIntensities percentOverlap distanceFromLast regionPixels regionSolidity centDiffSq percentOverlap
                for thisRegion = 1:length(area_vector)
                    regionIntensities(thisRegion,1) = sum(frameDiff(s(thisRegion).PixelIdxList));
                    regionPixels{thisRegion,1} = s(thisRegion).PixelIdxList;
                    regionSolidity{thisRegion,1} = s(thisRegion).Solidity;
                    if fn>1
                    centDiffSq = ([cents{fn-1,wellNum}]-[s(thisRegion).WeightedCentroid]).^2;
                    distanceFromLast(thisRegion,1) = sqrt(centDiffSq(1)+centDiffSq(2));
                    percentOverlap(thisRegion,1) = numel(intersect(areas{fn-1,wellNum},s(thisRegion).PixelIdxList))/numel(s(thisRegion).PixelIdxList); 
                    end
                end

                if fn>1
                [~,regIntMax]=max(area_vector'.*regionIntensities);
                [~,distMin]=min(distanceFromLast);
                [~,overlapMax]=max(percentOverlap);
                solidity =[s(:).Solidity];
                [~,maxSolidity]=max(solidity);
                idx = mode([regIntMax distMin overlapMax]);
    %             idx = mode([regIntMax distMin maxSolidity]);
                else
                [~,regIntMax]=max(area_vector'.*regionIntensities);
                idx = regIntMax; 
                end

                centroids(1,:) = [s(idx).WeightedCentroid];

                lastArea = s(idx).PixelIdxList;
                thisCircularity = s(idx).Circularity;
                cents{fn,wellNum}=centroids;
                areas{fn,wellNum}=lastArea;
                circularity{fn,wellNum}=thisCircularity;
                % calculate long-axis polynomial fit and curvature
                
%                     if thisCircularity>0.7 % for cases when larva is rolled too tightly into ball, perform image erosion
% 
%                         se = strel('disk',4);
%                         bw3 = imerode(bw2,se);
%                         bw3 = bwskel(bw2);
%                         imshow(bw3)
%                         bwMov(:,:,fn)=bw3;
% 
%                         L = bwlabel(bw3); % label each region that passed the bwareaopen critera
%                         s = regionprops(bw3,frameDiff,'PixelIdxList','WeightedCentroid'); % calculate the area, centroid, and axis length of each region
% %                         ft = fittype( 'poly2' );
%                     else
%                         bwMov(:,:,fn)=bw2;
%                     end
%                     
                    
%                 bwMov(:,:,fn)=bw2;
            
                larvaFrame = zeros(size(fFilterCrop));

                roundedCent = round(s(idx).WeightedCentroid);
                larvaFrame(s(idx).PixelIdxList)=1;
                larvInds = find(larvaFrame);
              
                [rowData,colData]=ind2sub(size(larvaFrame),larvInds);

                [vert,theta,a]=vertRot(rowData,colData);
                rotIm = imrotate(larvaFrame,theta+90);
%                 rotIm = imrotate(frameDiff,theta+90);

                rotFrame = imrotate(fFilterCrop,theta+90);

                larvIndsRot = find(rotIm);

                [rowData,colData]=ind2sub(size(rotIm),larvIndsRot);
                [xData, yData] = prepareCurveData(colData,rowData);

                ft = fittype( 'poly3' );
                [fitresult, gof] = fit( xData, yData, ft );
                fitData{fn,1}=fitresult;
                fitData{fn,2}=xData;
                fitData{fn,3}=yData;
                
                xcoords = round(xData);
                ycoords = round(fitresult(xData));
                
                xCurve = (xData);
                yCurve = (fitresult(xData));
                
                Vertices = [xCurve yCurve];
                Vertices = unique(Vertices,'rows');

                Lines=[(1:size(Vertices,1))' (2:size(Vertices,1)+1)']; Lines(end,2)=1;
                k=abs(LineCurvature2D(Vertices,Lines));

                meanCurvatures(fn,wellNum)=mean(k);
                maxCurvatures(fn,wellNum)=max(k);
                sumCurvatures(fn,wellNum)=sum(k);
                gofFits(fn,wellNum)=gof.rsquare;

%                 figure,  hold on;
%                 N=LineNormals2D(Vertices,Lines);
%                 k=k*10;
%                 plot([Vertices(:,1) Vertices(:,1)+k.*N(:,1)]',[Vertices(:,2) Vertices(:,2)+k.*N(:,2)]','g');
%                 plot([Vertices(Lines(:,1),1) Vertices(Lines(:,2),1)]',[Vertices(Lines(:,1),2) Vertices(Lines(:,2),2)]','b');
%                 axis equal;
                

                    rot2 = rotIm;
                    rotFrame2 = rotFrame;

                    fitMat = zeros(size(rotIm));
                    lineInds = sub2ind(size(fitMat),ycoords,xcoords);
                    fitMat(lineInds)=1;
                    fitMat = bwmorph(fitMat,'bridge');

                    rotFrame2(fitMat)=0;

                    rot2(lineInds)=0;
                    rotFrame(lineInds)=255;
        %             imshow(rotIm);hold on
                    im3d = cat(3,rotIm,rot2,rot2);
                    im3dframe = cat(3,rotFrame,rotFrame2,rotFrame2);

                    im3drot = imrotate(im3d,-1*(theta+90));
                    im3drotFrame = imrotate(im3dframe,-1*(theta+90));
        %             rotFrame = imrotate(thisFrame,theta+90);

                    rprops = regionprops(im3drot(:,:,3),'Centroid');
                    rotCent = round(rprops.Centroid);
        %             im3drotCrop = imcrop(im3drot,[rotCent(1)-20 rotCent(2)-20 40 40]);
                    im3drotCrop = imcrop(im3drotFrame,[rotCent(1)-40 rotCent(2)-40 80 80]);
%                     curvMov(:,:,:,fn)=im3drotCrop; 

            if showLiveTracking
                xcoor=centroids(:,1);
                xcoor(xcoor==0)=[];
                ycoor=centroids(:,2);
                ycoor(ycoor==0)=[];
%                 imshow(fFilterCrop)
                imshow(bw2)
                hold on
                scatter(xcoor,ycoor,10,'+','r')
                drawnow;refresh
                hold off
                thisFrame = getframe;
                trackMov(:,:,:,fn)=thisFrame.cdata;
                fn
            end
            
            if showCurvAxis

                curveMov(:,:,:,fn)=im3drotCrop;
                
                imshow(im3drotCrop)
                drawnow;refresh;

    %                         waitforbuttonpress

    %             plot(xData,fitresult(xData))

            end
              

            catch % usually throws if no fish in well
                  
                cents{fn,wellNum}=[0 0];
                areas{fn,wellNum}=0;
                meanCurvatures(fn,wellNum)=NaN;
                maxCurvatures(fn,wellNum)=NaN;
                sumCurvatures(fn,wellNum)=NaN;
                circularity{fn,wellNum}=NaN;
                gofFits(fn,wellNum)=NaN;
            end
                  
        end
%         toc
        
        %   toc
        %   cents{wn} = centroids;
    end  
    toc

    centsMult=[cents];
    
    thisAviFile = (aviFiles(aviNum).name);
    thisCentFileName = thisAviFile;
    thisCentFileName(strfind(thisCentFileName,'.'):end)=[];
    thisCentFileName = [thisCentFileName '_centroids.mat'];
    save(thisCentFileName,'centsMult','wellBackgrounds','meanCurvatures','maxCurvatures','sumCurvatures','circularity');
    
    clear cents mask pixes centsMult wellBackgrounds
    end
    
end

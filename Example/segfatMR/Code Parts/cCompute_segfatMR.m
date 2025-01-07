classdef cCompute_segfatMR<cCompute
    properties
        %% these properties must exist
        pVersion_cCompute_segfatMR = '1.2';
        
        %% these properties are just for application specific use
        oBoundaries = boundary.empty;
        pSelectedBound = '';
        pSliceDone = logical.empty;    % 0 means not semented, 1 means segemented
        pUseSlice = false;
        pLoc1 = '';
        pLoc2 = '';
        pFatThresh = 0;
    end
    
    methods(Static)
        
    end
    
    methods
        % % % data loading % % %
        function oComp = mInit_oCompApp(oComp, imgs, oCont)
            pFcfg = oCont.pFcfg;
            pAcfg = oCont.pAcfg;
            %-------------------------------%
            % sort images and init dat struct
            sliceLocs = arrayfun(@(x) x.sliceLocation, imgs, 'un', false);
                % determine SpacingBetweenSlices
            sliceLocsSort = sort(unique(cellfun(@str2num, sliceLocs)), 'ascend');
            SpacingBetweenSlices = diff(sliceLocsSort);
            SpacingBetweenSlices(end+1) = SpacingBetweenSlices(end);
            SpacingBetweenSlices = round(SpacingBetweenSlices*1000)/1000;   % round to 0.000
            if sum(SpacingBetweenSlices-min(SpacingBetweenSlices))~=0
                uiwait(warndlg('Spacings between slices are not consistent! Do not ingnore this message!'));
            else
                SpacingBetweenSlices = SpacingBetweenSlices(1);
            end
            for i =1:numel(imgs)
                if ~isfield(imgs(i).dicomInfo, 'SpacingBetweenSlices')
                    imgs(i).dicomInfo.SpacingBetweenSlices = SpacingBetweenSlices;
                end
            end
                % END determine SpacingBetweenSlices
            i=0;
            for sliceLoc = unique(sliceLocs)
                i=i+1;
                indImg = find(ismember(sliceLocs, sliceLoc(1)));
                oComp(i) = feval(str2func(pAcfg.cComputeFcn));
                oComp(i).oImgs = imgs(indImg);
                oComp(i).pPatientName = imgs(indImg(1)).patientName;
                oComp(i).pDataName = oComp(i).pPatientName;
                oComp(i).pSliceLocation = imgs(indImg(1)).sliceLocation;
                oComp(i).pSliceDone = false;
            end
            [a sortInd] = sort(cellfun(@(x) str2num(x) ,{oComp.pSliceLocation}), 'ascend');
            oComp = oComp(sortInd);
            oComp = setStructArrayField(oComp, 'pLoc1', {pAcfg.table.columnFormat{5}{1}});
            oComp = setStructArrayField(oComp, 'pLoc2', {pAcfg.table.columnFormat{6}{1}});
        end
        
        % % % GUI update % % %magn
        function imgDisplay = mGetImg2Display(oComp, oCont)
            pFcfg = oCont.pFcfg;
            pAcfg = oCont.pAcfg;
            oComp = oCont.oComp(oCont.pTable.row);
            %-------------------------------%
            %% determine image to be shown
            stdImgType = pAcfg.standardImgType;
            switch pAcfg.imageDisplayMode
                case 'Water only'
                    imgDisplay = oComp.mWaterImg;
                case 'Fat only'
                    imgDisplay = oComp.mFatImg;
                case 'Out Phase only'
                    imgDisplay = oComp.mGetImgOfType('OutPhase');
                case 'In Phase only'
                    imgDisplay = oComp.mGetImgOfType('InPhase');
                case 'All Four'
                    msgbox('All Four: not implementet!');
                case stdImgType
                    imgDisplay = oComp.mGetImgOfType(stdImgType);
                    
            end
            
            %% convert
            imgDisplay = imgDisplay.scale2([0 255]);
        end
                
        function mPostPlot(oComp, oCont)
            pFcfg = oCont.pFcfg;
            pAcfg = oCont.pAcfg;
            oComp = oCont.oComp(oCont.pTable.row);
            imgAxis = oCont.pHandles.imgAxis;
            %-------------------------------%
            %% plot segmentation
            try 
                disp(['RS_BodyBounds used fit method: ' num2str(oComp.pVarious.AutoSegment_ind)]);
            end
            contCoord = {};
            colors = {};
            for i=1:numel(oComp.oBoundaries)
                cBound = oComp.oBoundaries(i);
                contCoord{i} = {cBound.coord};
                colors(i) = pAcfg.contour.colors(find(ismember(pAcfg.contour.names, cBound.name)));
            end
            % plot contours
            oCont.pHandles.contour = oCont.mDrawContour(imgAxis, contCoord, colors);
            %% plot line for miror axis
            if pAcfg.contour.showHemiLine
                %search for hemi line
                hemiLineInd = find(arrayfun(@(x) isfield(x.pVarious, 'hemiLine'), oCont.oComp));
                if isempty(hemiLineInd)
                        oComp.mShowHemiLine(oCont);
                else
                    coord = oCont.oComp(hemiLineInd).pVarious.hemiLine.pos;
                    oCont.mDrawContour(imgAxis, {fliplr(coord)}, {'white'});
                end
                
            end
            %% plot rect
            if pAcfg.contour.showFemurBox
                femurBoxInd = find(arrayfun(@(x) isfield(x.pVarious, 'femurBox'), oCont.oComp));
                if isempty(femurBoxInd)
                        oComp.mShowFemurBox(oCont);
                else
                    coord = bwboundaries(oCont.oComp(femurBoxInd).pVarious.femurBox);
                    pos = [min(coord{1}(:,2)) min(coord{1}(:,1)) max(coord{1}(:,2))-min(coord{1}(:,2)) max(coord{1}(:,1))-min(coord{1}(:,1))];
                    
                    oCont.mDrawContour(imgAxis, {coord}, {[1 1 1]});
                    
                    T1 = text(0, 0, '', 'Parent', oCont.pHandles.imgAxis, 'Color', 'white');
                    T2 = text(0, 0, '', 'Parent', oCont.pHandles.imgAxis, 'Color', 'white');
                    T1.String = num2str(round(pos(3)*oComp(1).oImgs(1).dicomInfo.PixelSpacing(2)));
                    T1.Position = [pos(1)+pos(3)/2-5, pos(2)-6];
                    T2.String = num2str(round(pos(4)*oComp(1).oImgs(1).dicomInfo.PixelSpacing(2)));
                    T2.Position = [pos(1)+pos(3)+2, pos(2)+pos(4)/2];
                end
                
%                 try 
%                     coord = bwboundaries(oCont.oComp(oCont.pVarious.femurBoxInd).pVarious.femurBox);
%                     oCont.mDrawContour(imgAxis, {coord}, {[1 1 1]});
%                     
%                     pos = oCont.pVarious.femurBoxSize;
%                     T1 = text(0, 0, '', 'Parent', oCont.pHandles.imgAxis, 'Color', 'white');
%                     T2 = text(0, 0, '', 'Parent', oCont.pHandles.imgAxis, 'Color', 'white');
%                     T1.String = num2str(round(pos(3)*oComp(1).oImgs(1).dicomInfo.PixelSpacing(2)));
%                     T1.Position = [pos(1)+pos(3)/2-5, pos(2)-6];
%                     T2.String = num2str(round(pos(4)*oComp(1).oImgs(1).dicomInfo.PixelSpacing(2)));
%                     T2.Position = [pos(1)+pos(3)+2, pos(2)+pos(4)/2];
%                 catch
%                     uiwait(msgbox('could not draw box'));
%                     oComp.mShowFemurBox(oCont);
%                 end
            end
            %% show VAT?
            if any(ismember(oCont.pActiveKey, 'v'))
                oComp.mVatOverlay(oCont);
            end
            %% TopLeft text in Image
            pos = [3, 6];
            gap = 8;
            letterSize = 10;
            letterStyle = 'bold';
            imgText = [pAcfg.imageDisplayMode ' - Slice ' num2str(oCont.pTable.row)];
            oCont.pHandles.imgText = text(pos(1),pos(2), imgText, 'Parent', imgAxis, 'Color', 'white'); pos(2) = pos(2)+gap;
            oCont.pHandles.imgText.FontSize = letterSize;
            oCont.pHandles.imgText.FontWeight = letterStyle;
            oCont.pHandles.imgText.HitTest = 'off';
            % boundary text im image
            for i = 1:numel(oComp.oBoundaries)
                if ~isempty(oComp.oBoundaries(i).coord)
                    % plot annotation
                    txt = oComp.oBoundaries(i).name;
                    t = text(pos(1),pos(2), txt, 'Parent', oCont.pHandles.imgAxis, 'Color', colors{i});
                    t.FontSize = letterSize;
                    t.FontWeight = letterStyle;
                    pos(2) = pos(2)+gap;
                end
            end
        end
        
        function oCont = mSurveyViewer(oComp, oCont)
            
            if isfield(oCont.pHandles, 'SV') && isvalid(oCont.pHandles.SV) && ishandle(oCont.pHandles.SV.handles.hline)
                oComp.mUpdateSurveyViewer(oCont);    
            else
                oComp.mStartSurveyViewer(oCont);    
            end

        end
        
        function oCont = mStartSurveyViewer(oComp, oCont)
            oCont.pHandles.SV = SurveyTool(oCont, oComp);   
        end
        
        function oCont = mUpdateSurveyViewer(oComp, oCont)
            disp('SV is open');
            ReflineH = oCont.pHandles.SV.handles.hline;
            
            NewSlice = oComp(oCont.pTable.row);
            SlicePos = NewSlice.oImgs(1).dicomInfo.ImagePositionPatient(3)+NewSlice.oImgs(1).dicomInfo.SpacingBetweenSlices/2;
            
            
            SlicePos = ((oCont.pHandles.SV.oImgs.dicomInfo.ImagePositionPatient(3)-oCont.pHandles.SV.oImgs.dicomInfo.PixelSpacing(2)/2)-SlicePos)/oCont.pHandles.SV.oImgs.dicomInfo.PixelSpacing(1);
            
            zPos = SlicePos;
            YDataNew = set(ReflineH, 'YData',[zPos, zPos]);
        end
        
        function oCont = mDrawGraph(oComp, oCont)
            pFcfg = oCont.pFcfg;
            pAcfg = oCont.pAcfg;
            datInd = oCont.pTable.row;
            oComp = oCont.oComp(datInd);
            graphAxis = oCont.pHandles.graphAxis;
            axes(graphAxis);
            visInd = oComp.oBoundaries.getBoundInd('visceralBound');
            %-------------------------------%
            % determine graphData to be plotted in graph
            if visInd==0
                graphData = oComp.mWaterImg.data;
            else
                visBounoCont = oComp.oBoundaries(visInd);
                VatMask = logical(oComp.mGetBoundMask(oComp.oImgs(1).data, oComp.oBoundaries(visInd).coord));
                graphData = oComp.mWaterImg.data(VatMask);
            end
            
            % plot hist with as many bars as the data has possible values
            if isfield(oCont.pHandles, 'graphPlot') && ishandle(oCont.pHandles.graphPlot)   % use current histogram
                %[histDat histInd] = histcounts(graphData, 'NumBins', max(graphData)/10)
                oCont.pHandles.graphPlot.Data = graphData;
                oCont.pHandles.graphPlot.NumBins = max(max(graphData))/10;
                oCont.pHandles.graphPlot.BinLimitsMode = 'manual';
                oCont.pHandles.graphPlot.BinLimits = [min(min(graphData)) max(max(graphData))];
                oCont.pHandles.graphAxis.YLim = [0 max(oCont.pHandles.graphPlot.Values)];
                
            else    % create histogram because not created, yet
                hold all;
                oCont.pHandles.graphPlot = histogram(graphData);
                oCont.pHandles.graphPlot.NumBins = max(max(graphData))/10;
                oCont.pHandles.graphPlot.BinLimitsMode = 'manual';
                oCont.pHandles.graphPlot.BinLimits = [min(min(graphData)) max(max(graphData))];
                oCont.pHandles.graphPlot.EdgeColor = oCont.pFcfg.apperance.designColor1;
                oCont.pHandles.graphPlot.EdgeAlpha = 0.66;
                oCont.pHandles.graphPlot.FaceColor = oCont.pFcfg.apperance.designColor1;
                oCont.pHandles.graphPlot.FaceAlpha = 0.2;
                oCont.pHandles.graphAxis.YTick = [];
                oCont.pHandles.graphAxis.XTickMode = 'auto';
                oCont.pHandles.graphAxis.FontWeight = 'bold';
                oCont.pHandles.graphAxis.Box = 'on';
                oCont.pHandles.graphPlot.HitTest = 'off';
                oCont.pHandles.graphAxis.YLabel.String = 'frequency';
                oCont.pHandles.graphAxis.XLabel.String = 'pixel intensity';
                hold off;
            end
            
            if visInd==0   % set hist line to fatThresh value
                try delete(oCont.pHandles.histLine); end
            else
                % update histogram threshold line
                oCont = oComp.mUpdateHistLine(oCont);
            end
            drawnow
            
            %% check if slice is done
            oComp = oComp.mSet_pSliceDone;
            %-------------------------------%
            oCont.oComp(datInd) = oComp;
        end
        
        function oCont = mUpdateHistLine(oComp, oCont)
            if isfield(oCont.pHandles, 'histLine') && isvalid(oCont.pHandles.histLine)    % use existing line
                oCont.pHandles.histLine.XData = [oComp.pFatThresh oComp.pFatThresh];
                oCont.pHandles.histLine.YData = [0 oCont.pHandles.graphAxis.YLim(2)];
            else    % create new line
                oCont.pHandles.histLine = line([oComp.pFatThresh oComp.pFatThresh], [0 oCont.pHandles.graphAxis.YLim(2)], 'parent', oCont.pHandles.graphAxis);
                oCont.pHandles.histLine.LineWidth = 3;
                oCont.pHandles.histLine.LineStyle = ':';
                oCont.pHandles.histLine.Color = 'green';
                oCont.pHandles.histLine.HitTest = 'off';
            end
        end
        
        function lines = mGetTextBoxLines(oComp, oCont)
            lines = string('');
        end
        
        % % % GUI Interaction % % %
        function oComp = mTableEdit(oComp, select)
            switch select.Source.ColumnName{select.Indices(2)}
                case 'UseIt'
                    oComp.pUseSlice = select.NewData;
            end
        end
        
        function oComp = mKeyPress(oComp, oCont, key)
            pAcfg = oCont.pAcfg;
            key = key.Key;
            %-------------------------------%
            oCont.pHandles.figure.WindowKeyReleaseFcn = {@oCont.mKeyRelease};
            %oCont.pHandles.figure.WindowKeyPressFcn = '';
            switch key
                case pAcfg.contour.keyAssociation  % normal grayed out contour display
                    oComp.pSelectedBound = pAcfg.contour.names(find(ismember(pAcfg.contour.keyAssociation, key)));
                    imgAxis = oCont.pHandles.imgAxis;
                    % delete axis childs
                    if exist('oCont.pHandles.contour')
                        for a = oCont.pHandles.contour
                            delete(a);
                        end
                    end
                    
                    %% Plot Contours grayed out
                    contCoord = {};
                    colors = {};
                    for i=1:numel(oComp.oBoundaries)
                        cBound = oComp.oBoundaries(i);
                        contCoord{i} = {cBound.coord};
                        if isequal({cBound.name}, oComp.pSelectedBound)
                            colors(i) = pAcfg.contour.colors(find(ismember(pAcfg.contour.names, oComp.pSelectedBound)));
                        else
                            colors{i} = [0.1 0.1 0.1];
                        end
                    end
                    % plot contours
                    oCont.pHandles.contour = oCont.mDrawContour(imgAxis, contCoord, colors);
                    
                    pause(0);
                    
                case pAcfg.key.deleteContour  % delete selected contour
                    disp('keyDelCont')
                    % delete contour if one contour is selected
                    otherKeys = oCont.pActiveKey(~ismember(oCont.pActiveKey, {key}))
                    if numel(otherKeys) > 1
                        disp('to many keys pressed');
                    else
                        switch otherKeys{1}
                            case pAcfg.contour.keyAssociation
                                contourName = pAcfg.contour.names(find(ismember(oCont.pAcfg.contour.keyAssociation, otherKeys{1})))
                                
                                boundInd = oComp.oBoundaries.getBoundInd(oComp.pSelectedBound);
                                if ~isempty(boundInd)
                                    oComp.oBoundaries(boundInd) = [];
                                    oComp.pSliceDone = false;
                                    oComp.pUseSlice = false;
                                end
                                
                        end
                    end

                case pAcfg.key.showVat  % show vat
                    oComp.mVatOverlay(oCont);
                case pAcfg.key.contourTracking
                    pAcfg.contour.contourTracking.enable = ~pAcfg.contour.contourTracking.enable;
            end
            %-------------------------------%
            oCont.pAcfg = pAcfg;
        end
        
        function mKeyRelease(oComp, oCont, keys)
            oCont.pHandles.figure.WindowKeyPressFcn = @oCont.mKeyPress;
            
            
            oCont.pHandles.imgDisplay.AlphaData = 1;
            oCont.mTableCellSelect('Caller', 'mKeyRelease');
        end
        
        function mImgAxisButtonDown(oComp, oCont, hit)
            switch oCont.pActiveKey{1}
                case ''
                otherwise
                oCont.pLineCoord = [oCont.pHandles.imgAxis.CurrentPoint(1,2) oCont.pHandles.imgAxis.CurrentPoint(1,1)];
                oCont.pActiveMouse = oCont.pHandles.figure.SelectionType;
                
                oCont.pHandles.draw = oCont.mDrawContour(oCont.pHandles.imgAxis, {{oCont.pLineCoord}}, {'green'});
                oCont.pHandles.draw = oCont.pHandles.draw{1};
                if isempty(oCont.pHandles.figure.WindowButtonMotionFcn) | isempty(oCont.pHandles.figure.WindowButtonMotionFcn)
                    oCont.pHandles.figure.WindowButtonMotionFcn = {@oComp.mImgAxisButtonMotion, oCont};
                    oCont.pHandles.figure.WindowButtonUpFcn = {@oComp.mImgAxisButtonUp, oCont};
                else
                    msgbox('WindowButtonMotionFcn and WindowButtonUpFcn are already set');
                end
            end
        end
        
        function mImgAxisButtonMotion(oComp, a, b, oCont)
            newC = [oCont.pHandles.imgAxis.CurrentPoint(1,2) oCont.pHandles.imgAxis.CurrentPoint(1,1)];
            newC = round(newC);
            if oCont.pAcfg.contour.contourTracking.enable
                img = oComp.mGetImgOfType('OutPhase').data;
                mask = zeros(size(img));
                mask(newC(1), newC(2)) = 1;
                %trackImg = double(max(max(img))-img);
                trackImg = imgradient(img);
                trackSize = oCont.pAcfg.contour.contourTracking.size;
                mask = trackImg.*imgaussfilt(mask, trackSize, 'FilterSize', trackSize*8+1);
                [v ind] = max(mask(:));
                [x y]=ind2sub(size(img), ind);
                newC = [x y];
            end
            oCont.pLineCoord = [oCont.pLineCoord; newC];
            oCont.pHandles.draw.XData = oCont.pLineCoord(:,2);
            oCont.pHandles.draw.YData = oCont.pLineCoord(:,1);
            drawnow;
        end
        
        function mImgAxisButtonUp(oComp, a, b, oCont)
            oCont.pHandles.figure.WindowButtonMotionFcn = '';
            oCont.pHandles.figure.WindowButtonUpFcn = '';
            oCont = oCont.oComp.mMergeContours(oCont);
        end
        
        function mGraphAxisButtonDown(oComp, oCont, hit)
            % calc VAT
            oComp.pFatThresh = oCont.pHandles.graphAxis.CurrentPoint(1,1);
            
            oCont.oComp(oCont.pTable.row) = oComp;
            
            oComp.mUpdateHistLine(oCont);
            
            oComp.mVatOverlay(oCont);
            
            if isempty(oCont.pHandles.figure.WindowButtonMotionFcn) | isempty(oCont.pHandles.figure.WindowButtonMotionFcn)
                oCont.pHandles.figure.WindowButtonMotionFcn = {@oComp.mGraphAxisButtonMotion, oCont};
                oCont.pHandles.figure.WindowButtonUpFcn = {@oComp.mGraphAxisButtonUp, oCont};
            else
                msgbox('WindowButtonMotionFcn and WindowButtonUpFcn are already set');
            end
        end
        
        function mGraphAxisButtonMotion(oComp, a, b, oCont)
            oComp.pFatThresh = oCont.pHandles.graphAxis.CurrentPoint(1,1);
            oComp.mVatOverlay(oCont);
            %-------------------------------%
            oCont.oComp(oCont.pTable.row) = oComp;
            oComp.mUpdateHistLine(oCont);
        end
        
        function mGraphAxisButtonUp(oComp, a, b, oCont)
            oCont.pHandles.figure.WindowButtonMotionFcn = '';
            oCont.pHandles.figure.WindowButtonUpFcn = '';
            oCont.pHandles.imgDisplay.AlphaData = 1;
            oCont.mTableCellSelect('Caller', 'mGraphAxisButtonUp');
        end
        
        % % % Experimental % % %
        % % Hemi Fat % %
        function oComp = mShowHemiLine(oComp, oCont)
            oCont.pAcfg.contour.showHemiLine = ~oCont.pAcfg.contour.showHemiLine;
            if oCont.pAcfg.contour.showHemiLine == 0
                uiwait(msgbox('Hemi Line not visible'));
            elseif oCont.pAcfg.contour.showFemurBox == 1
                uiwait(msgbox('Hemi Line now visible'));
            end
        end
        
        function oComp = mSetHemiLine(oComp, oCont)
            hLine = imline(oCont.pHandles.imgAxis);
            imgSize = size(oComp(1).mGetStandardImg.data);
            %-------------------------------%
            Message = text(0, 0, '', 'Parent', oCont.pHandles.imgAxis, 'Color', 'white');
            Message.String = 'Double click on line to confirm!';
            Message.Position = [oCont.pHandles.imgAxis.XLim(2)/2-50 oCont.pHandles.imgAxis.YLim(2)/2];
            Message.FontSize = 14;
            Message.FontWeight = 'bold';
            
            wait(hLine);
            
            % make mask
            % get points at upper and lower edge of image according to y=mx+t
            pos = hLine.getPosition;
            x1 = pos(1,1);
            x2 = pos(2,1);
            y1 = pos(1,2);
            y2 = pos(2,2);
            m = (y2-y1)/(x2-x1);
            t = y1-m*x1;
            y11 = 1;
            x11 = (y11-t)/m;
            y22 = imgSize(1);
            x22 = (y22-t)/m;
            
            [y x] = makeLinePoints(x11, y11, x22, y22);
            
            lineMask = zeros(imgSize);
            lineMask = oComp.mGetBoundMask(lineMask, [x' y']);
            
            leftMask = lineMask;
            leftMask(1, x11:end) = 1;
            leftMask(end, x22:end) = 1;
            leftMask(:, end) = 1;
            leftMask = imfill(leftMask,'holes');
            
            rightMask = lineMask;
            rightMask(1, 1:x11) = 1;
            rightMask(end, 1:x22) = 1;
            rightMask(:, 1) = 1;
            rightMask = imfill(rightMask,'holes');
            
            % store data
            % first delete old hemi line
            oComp = oComp.mDelHemiLine(oCont);
            
            datInd = oCont.pTable.row;
            oCompTmp = oCont.oComp(datInd);
            oCompTmp.pVarious.hemiLine.handle = hLine;
            oCompTmp.pVarious.hemiLine.ind = datInd;
            oCompTmp.pVarious.hemiLine.pos = pos;
            oCompTmp.pVarious.hemiLine.leftMask = leftMask;
            oCompTmp.pVarious.hemiLine.rightMask = rightMask;
            oComp(datInd) = oCompTmp;
            
            oCont.pAcfg.contour.showHemiLine = 1;
        end
        
        function oComp = mDelHemiLine(oComp, oCont)
            hemiLineInd = find(arrayfun(@(x) isfield(x.pVarious, 'hemiLine'), oCont.oComp));
            if ~isempty(hemiLineInd)
                for i =1:numel(hemiLineInd)
                    oComp(hemiLineInd(i)).pVarious = rmfield(oComp(hemiLineInd(i)).pVarious, 'hemiLine');
                end
            end
            oCont.pAcfg.contour.showHemiLine = 0;
        end
        
        function oComp = mSaveHemiResults(oComp, oCont)
            answer = questdlg('What should be saved?', 'Choose data saving:', 'only XLS', 'XLS and SessionFile', 'Abord', 'Abord');
            if strcmp(answer, 'Abord')
                return
            end
            if strcmp(answer, 'XLS and SessionFile')
                % save session file
                path = fullfile(oCont.pAcfg.lastLoadPath, [oCont.mGetSaveFilePrefix '_data.mat']);
                oCont.oComp.mSave_oComp(path, oCont);
            end
            if strcmp(answer, 'only XLS') | strcmp(answer, 'XLS and SessionFile')
                oComp.mSaveRightHemiResults(oCont);
                oComp.mSaveLeftHemiResults(oCont);
            end
        end
        
        function oComp = mSaveRightHemiResults(oComp, oCont)
            Name = 'Right';
            wb = waitbar(0, ['saving XLS ' Name ' hemiFat results']);
            wb.Name = 'saving....';
            hemiLineInd = find(arrayfun(@(x) isfield(x.pVarious, 'hemiLine'), oCont.oComp));
            hemiLinePos = oCont.oComp(hemiLineInd).pVarious.hemiLine.pos;
            mask = oCont.oComp(hemiLineInd).pVarious.hemiLine.rightMask;
            %% collect infos
                dicomInfo = oComp(1).mGetStandardImg.dicomInfo;
                try info.comment = dicomInfo.StudyComments; end
                try info.description = dicomInfo.RequestedProcedureDescription; end
                try info.physicianName = dicomInfo.ReferringPhysicianName.FamilyName; end
                try info.institution = dicomInfo.InstitutionName; end
                try info.stationName = dicomInfo.StationName; end
                try info.manufacturer = dicomInfo.Manufacturer; end
                try info.manufacturerModelName = dicomInfo.ManufacturerModelName; end
                
                try info.patientName = [dicomInfo.PatientName.FamilyName '_' dicomInfo.PatientName.GivenName];
                catch
                    try info.patientName = dicomInfo.PatientName.FamilyName;
                    catch
                        info.patientName = 'NoName';
                    end
                end
                try info.patientWeight = num2str(dicomInfo.PatientWeight); end
                try info.patientAge = dicomInfo.PatientAge; end
                try info.patientSex = dicomInfo.PatientSex; end
                try info.patientBirthDat = dicomInfo.PatientBirthDate; end
                try info.patientIoCont = dicomInfo.PatientID; end
                
                try info.creationDate = datestr(datenum(dicomInfo.InstanceCreationDate, 'yyyymmdd'), 'dd.mm.yyyy'); end
                
                % remove empty entries
                emptyInd = structfun(@isempty, info);
                infoFields = fieldnames(info);
                for i = 1:numel(emptyInd)
                    if emptyInd(i)
                        info = rmfield(info, infoFields(i));
                    end
                end
                
                waitbar(0.3, wb);
            %% create struct with oComp!
                img = oComp(1).mGetStandardImg;
                imgSize = size(img.data);
                    %% preparation
                xlsPath = fullfile(oCont.pAcfg.lastLoadPath, [oCont.mGetSaveFilePrefix 'HemiFat_' Name '_data.xlsx']);
                    %% read values from oComp to s struct
                for i = 1:numel(oComp)
                    oCompTmp = oComp(i);
                    %% modify Boundaries
                    if oCompTmp.pSliceDone
                    try
                        oBoundInd = oCompTmp.oBoundaries.getBoundInd('outerBound');
                        oBound = oCompTmp.oBoundaries(oBoundInd);
                        oBoundMask = oCompTmp.mGetBoundMask(img.data, oBound.coord);
                        oBoundMask = oBoundMask&mask;
                        tmp = bwboundaries(oBoundMask);
                        oBound.coord = tmp{1};
                    catch
                        oBound = [];
                    end
                    
                    try
                        iBoundInd = oCompTmp.oBoundaries.getBoundInd('innerBound');
                        iBound = oCompTmp.oBoundaries(iBoundInd);
                        iBoundMask = oCompTmp.mGetBoundMask(img.data, iBound.coord);
                        iBoundMask = iBoundMask&mask;
                        tmp = bwboundaries(iBoundMask);
                        iBound.coord = tmp{1};
                    catch
                        iBound = [];
                    end
                    
                    try
                        vBoundInd = oCompTmp.oBoundaries.getBoundInd('visceralBound');
                        vBound = oCompTmp.oBoundaries(vBoundInd);
                        vBoundMask = oCompTmp.mGetBoundMask(img.data, vBound.coord);
                        vBoundMask = vBoundMask&mask;
                        tmp = bwboundaries(vBoundMask);
                        vBound.coord = tmp{1};
                    catch
                        vBound = [];
                    end
                    
                    oCompTmp.oBoundaries = [oBound iBound vBound];
                    
                    end
                    
                    s.voxVol(i) = oCompTmp.oImgs(1).getVoxelVolume;
                    s.fatThresh(i) = oCompTmp.pFatThresh;
                    try 
                        oBoundCoord = oCompTmp.oBoundaries(oCompTmp.oBoundaries.getBoundInd('outerBound')); 
                    catch
                        oBoundCoord = []; 
                    end
                    try 
                        oBoundMask = oCompTmp.mGetBoundMask(oCompTmp.mWaterImg.data, oBoundcoord.coord); 
                    catch
                        oBoundMask = [];
                    end
                    tmp = regionprops(oBoundMask , 'Area', 'Perimeter');
                    if isempty(tmp)
                        tmp(1).Area = 0;
                        tmp(1).Perimeter = 0;
                    end
                    s.bodyArea(i) = tmp.Area;
                    s.bodyPerimeter(i) = tmp.Perimeter;
                    
                    s.SatArea(i) = oCompTmp.mVolumeSAT/oCompTmp.oImgs(1).dicomInfo.SpacingBetweenSlices*1000;
                    s.VatArea(i) = oCompTmp.mVolumeVAT/oCompTmp.oImgs(1).dicomInfo.SpacingBetweenSlices*1000;
                    s.SatVol(i) = oCompTmp.mVolumeSAT;
                    s.VatVol(i) = oCompTmp.mVolumeVAT;
                    s.Loc1(i) = {oCompTmp.pLoc1};
                    s.Loc2(i) = {oCompTmp.pLoc2};
                    s.sliceLoc(i) = str2num(oCompTmp.pSliceLocation);
                    s.sliceNr(i) = i;
                    if isempty(oCompTmp.pUseSlice)
                        oCompTmp.pUseSlice = false;
                    end
                    s.useSlice (i) = oCompTmp.pUseSlice;
                end
                [keepInd status message] = LMranging(oComp, oCont.pAcfg.xlsSaveRange);
                if status==1
                    answer = questdlg(message, 'Slice indexing error', 'Save marked slices', 'Abord', 'Abord')
                    switch answer
                        case 'Abord'
                            return
                        case 'Save marked slices'
                            keepInd = s.useSlice;
                    end
                end
                
                tmp = s.useSlice;
                s.useSlice = zeros(size(tmp));
                s.useSlice(keepInd) = tmp(keepInd);
                % postprocessing
                s.Loc1(ismember(s.Loc1, 'none')) = {''};
                s.Loc2(ismember(s.Loc2, 'none')) = {''};
                
                sFlip = structfun(@(x) x', s, 'Uniformoutput', 0);
                    %% write to xls sheet (use all available data)
                writetable(table(hemiLinePos), xlsPath, 'Sheet', 'Pos');
                writetable(struct2table(info), xlsPath, 'Sheet', 'infos');
                writetable(struct2table(sFlip), xlsPath, 'Sheet', 'allFatData');
                waitbar(0.6, wb);
                    %% write to xls sheet (xls file like <nikita 2016)
                % filter accordint to pUseSlice
                s.SatVol(~s.useSlice) = 0;
                s.VatVol(~s.useSlice) = 0;
            if oCont.pAcfg.doSliceSpacingInterpolation
                if ~isnan(oCont.pAcfg.sliceSpacingInterpolationDistance) && any(abs(diff(s.sliceLoc) - oCont.pAcfg.sliceSpacingInterpolationDistance) > 0.5)
                    % prepare data: interpolate to equidistant slice loc
                    x = s.sliceLoc;
                    sliceLocOld = x;
                    window = oCont.pAcfg.sliceSpacingInterpolationDistance;
                    xn = [x(1):window:x(end)];
                    sliceLocNew = xn;
                    % find correct slice width (oContescription available (oContrawIO)) for volume calculation
                    sliceWidth = diff(sliceLocNew);
                    sW1 = sliceWidth./2; sW1(end+1) = 0;
                    sW2 = sliceWidth./2; sW2 = [0 sW2];
                    sliceWidth = sW1+sW2;
                    % now special treatment for boundary images (1 and end)
                    sliceWidth(1) = sliceWidth(1)+str2num(oCont.oComp(1).oImgs(1).sliceThickness)/2;
                    sliceWidth(end) = sliceWidth(end)+str2num(oCont.oComp(end).oImgs(end).sliceThickness)/2;
                    % slice Width done!
                    
                    %_Calc VatVol
                    y = s.VatArea;  % mm^2
                    [x2 y2] = DconvV2(x,y,window,'xmeanymean',[]); x2 = x2{1}; y2 = y2{1};
                    ynVatArea = interp1(x2, y2, sliceLocNew);
                    VatVol = ynVatArea.*sliceWidth*0.001;
                    
                    %_Calc SatVol
                    y = s.SatArea;  % mm^2
                    [x2 y2] = DconvV2(x,y,window,'xmeanymean',[]); x2 = x2{1}; y2 = y2{1};
                    ynSatArea = interp1(x2, y2, sliceLocNew);
                    SatVol = ynSatArea.*sliceWidth*0.001;
                    
                    %_Set WKs
                    LocsInd = find(~ismember(s.Loc1,''));
                    LocsVal = s.Loc1(LocsInd);
                    LocsPosOrig = s.sliceLoc(LocsInd);
                    Loc1 = cell(1, numel(sliceLocNew)); Loc1(:) = {''};
                    for i = 1:numel(LocsPosOrig)
                        [val ind] = min(abs(sliceLocNew-LocsPosOrig(i)));
                        Loc1(ind) = LocsVal(i);
                    end
                    
                    %_Set Landmarks
                    LocsInd = find(~ismember(s.Loc2,''));
                    LocsVal = s.Loc2(LocsInd);
                    LocsPosOrig = s.sliceLoc(LocsInd);
                    Loc2 = cell(1, numel(sliceLocNew)); Loc2(:) = {''};
                    for i = 1:numel(LocsPosOrig)
                        [val ind] = min(abs(sliceLocNew-LocsPosOrig(i)));
                        Loc2(ind) = LocsVal(i);
                    end
                    
                    %_Set sliceNr
                    sliceNr = 1:numel(sliceLocNew);
                    
                    %_Plot Interpolation Reults
                    figure();
                    plot(diff(x),'LineStyle', 'none', 'Marker', 'o', 'Color', 'black', 'DisplayName', 'raw');
                    hold on
                    plot(diff(x2),'LineStyle', 'none', 'Marker', 'x', 'Color', 'black', 'DisplayName', 'smooth');
                    plot(diff(xn),'LineStyle', 'none', 'Marker', '.', 'Color', 'red', 'DisplayName', 'interpolated');
                    drawnow;
                    a = gca;
                    a.XLabel.String = 'sliceLocation';
                    a.YLabel.String = 'VatArea';
                    legend('show');
                    
                    figure();
                    plot(x, y, 'LineStyle', 'none', 'Marker', 'o', 'Color', 'black', 'DisplayName', 'raw');
                    hold on
                    plot(x2,y2,'LineStyle', 'none', 'Marker', 'x', 'Color', 'black', 'DisplayName', 'smooth');
                    plot(xn,ynSatArea,'LineStyle', 'none', 'Marker', '.', 'Color', 'red', 'DisplayName', 'interpolated');
                    drawnow;
                    a = gca;
                    a.XLabel.String = 'sliceLocation';
                    a.YLabel.String = 'VatArea';
                    legend('show');
                    
                    
                    
                    
                    waitfor(msgbox('Slice locations are not equally distributed -> interpolation was done as shown in figure!'));
                    
                    
                else
                    sliceLocNew = s.sliceLoc;
                    SatVol = s.SatVol;
                    VatVol = s.VatVol;
                    sliceNr = s.sliceNr;
                    Loc1 = s.Loc1;
                    Loc2 = s.Loc2;
                end
            else
                sliceLocNew = s.sliceLoc;
                SatVol = s.SatVol;
                VatVol = s.VatVol;
                sliceNr = s.sliceNr;
                Loc1 = s.Loc1;
                Loc2 = s.Loc2;    
            end
                
                Pos = 1;
                xlswrite(xlsPath, {info.patientName} , 'FAT', 'A1');
                xlswrite(xlsPath, {info.creationDate} , 'FAT', 'A2');
                xlswrite(xlsPath, {'Slice #' 'Slice Pos (S)' 'Slice Gap' 'SAT [cm^3]' 'VAT [cm^3]' 'Loc1' 'Loc2'} , 'FAT', 'A3');
                xlswrite(xlsPath, sliceNr', 'FAT', 'A4');
                xlswrite(xlsPath, sliceLocNew', 'FAT', 'B4');
                xlswrite(xlsPath, [0 diff(sliceLocNew)]', 'FAT', 'C4');
                xlswrite(xlsPath, round(SatVol)', 'FAT', 'D4');
                xlswrite(xlsPath, round(VatVol)', 'FAT', 'E4');
                xlswrite(xlsPath, Loc1', 'FAT', 'F4');
                xlswrite(xlsPath, Loc2', 'FAT', 'G4');
                Pos = 3+numel(VatVol)+1;
                xlswrite(xlsPath, {'Slice #' 'Slice Pos (S)' 'Slice Gap' 'SAT [cm^3]' 'VAT [cm^3]' 'Loc1' 'Loc2'} , 'FAT', ['A' num2str(Pos)]);
                Pos = Pos+2;
                xlswrite(xlsPath, {'Summe'}, 'FAT', ['A' num2str(Pos)]);
                xlswrite(xlsPath, round(sum(SatVol)), 'FAT', ['D' num2str(Pos)]);
                xlswrite(xlsPath, round(sum(VatVol)), 'FAT', ['E' num2str(Pos)]);
                waitbar(0.9, wb);
                pause(0.3);
                close(wb);
        end
        
        function oComp = mSaveLeftHemiResults(oComp, oCont)
            Name = 'Left';
            wb = waitbar(0, ['saving XLS ' Name ' hemiFat results']);
            wb.Name = 'saving....';
            hemiLineInd = find(arrayfun(@(x) isfield(x.pVarious, 'hemiLine'), oCont.oComp));
            hemiLinePos = oCont.oComp(hemiLineInd).pVarious.hemiLine.pos;
            mask = oCont.oComp(hemiLineInd).pVarious.hemiLine.leftMask;
            %% collect infos
                dicomInfo = oComp(1).mGetStandardImg.dicomInfo;
                try info.comment = dicomInfo.StudyComments; end
                try info.description = dicomInfo.RequestedProcedureDescription; end
                try info.physicianName = dicomInfo.ReferringPhysicianName.FamilyName; end
                try info.institution = dicomInfo.InstitutionName; end
                try info.stationName = dicomInfo.StationName; end
                try info.manufacturer = dicomInfo.Manufacturer; end
                try info.manufacturerModelName = dicomInfo.ManufacturerModelName; end
                
                try info.patientName = [dicomInfo.PatientName.FamilyName '_' dicomInfo.PatientName.GivenName];
                catch
                    try info.patientName = dicomInfo.PatientName.FamilyName;
                    catch
                        info.patientName = 'NoName';
                    end
                end
                try info.patientWeight = num2str(dicomInfo.PatientWeight); end
                try info.patientAge = dicomInfo.PatientAge; end
                try info.patientSex = dicomInfo.PatientSex; end
                try info.patientBirthDat = dicomInfo.PatientBirthDate; end
                try info.patientIoCont = dicomInfo.PatientID; end
                
                try info.creationDate = datestr(datenum(dicomInfo.InstanceCreationDate, 'yyyymmdd'), 'dd.mm.yyyy'); end
                
                % remove empty entries
                emptyInd = structfun(@isempty, info);
                infoFields = fieldnames(info);
                for i = 1:numel(emptyInd)
                    if emptyInd(i)
                        info = rmfield(info, infoFields(i));
                    end
                end
                
                waitbar(0.3, wb);
            %% create struct with oComp!
                img = oComp(1).mGetStandardImg;
                imgSize = size(img.data);
                    %% preparation
                xlsPath = fullfile(oCont.pAcfg.lastLoadPath, [oCont.mGetSaveFilePrefix 'HemiFat_' Name '_data.xlsx']);
                    %% read values from oComp to s struct
                for i = 1:numel(oComp)
                    oCompTmp = oComp(i);
                    %% modify Boundaries
                    if oCompTmp.pSliceDone
                    try
                        oBoundInd = oCompTmp.oBoundaries.getBoundInd('outerBound');
                        oBound = oCompTmp.oBoundaries(oBoundInd);
                        oBoundMask = oCompTmp.mGetBoundMask(img.data, oBound.coord);
                        oBoundMask = oBoundMask&mask;
                        tmp = bwboundaries(oBoundMask);
                        oBound.coord = tmp{1};
                    catch
                        oBound = [];
                    end
                    
                    try
                        iBoundInd = oCompTmp.oBoundaries.getBoundInd('innerBound');
                        iBound = oCompTmp.oBoundaries(iBoundInd);
                        iBoundMask = oCompTmp.mGetBoundMask(img.data, iBound.coord);
                        iBoundMask = iBoundMask&mask;
                        tmp = bwboundaries(iBoundMask);
                        iBound.coord = tmp{1};
                    catch
                        iBound = [];
                    end
                    
                    try
                        vBoundInd = oCompTmp.oBoundaries.getBoundInd('visceralBound');
                        vBound = oCompTmp.oBoundaries(vBoundInd);
                        vBoundMask = oCompTmp.mGetBoundMask(img.data, vBound.coord);
                        vBoundMask = vBoundMask&mask;
                        tmp = bwboundaries(vBoundMask);
                        vBound.coord = tmp{1};
                    catch
                        vBound = [];
                    end
                    
                    oCompTmp.oBoundaries = [oBound iBound vBound];
                    
                    end
                    
                    s.voxVol(i) = oCompTmp.oImgs(1).getVoxelVolume;
                    s.fatThresh(i) = oCompTmp.pFatThresh;
                    try 
                        oBoundCoord = oCompTmp.oBoundaries(oCompTmp.oBoundaries.getBoundInd('outerBound')); 
                    catch
                        oBoundCoord = []; 
                    end
                    try 
                        oBoundMask = oCompTmp.mGetBoundMask(oCompTmp.mWaterImg.data, oBoundcoord.coord); 
                    catch
                        oBoundMask = [];
                    end
                    tmp = regionprops(oBoundMask , 'Area', 'Perimeter');
                    if isempty(tmp)
                        tmp(1).Area = 0;
                        tmp(1).Perimeter = 0;
                    end
                    s.bodyArea(i) = tmp.Area;
                    s.bodyPerimeter(i) = tmp.Perimeter;
                    
                    s.SatArea(i) = oCompTmp.mVolumeSAT/oCompTmp.oImgs(1).dicomInfo.SpacingBetweenSlices*1000;
                    s.VatArea(i) = oCompTmp.mVolumeVAT/oCompTmp.oImgs(1).dicomInfo.SpacingBetweenSlices*1000;
                    s.SatVol(i) = oCompTmp.mVolumeSAT;
                    s.VatVol(i) = oCompTmp.mVolumeVAT;
                    s.Loc1(i) = {oCompTmp.pLoc1};
                    s.Loc2(i) = {oCompTmp.pLoc2};
                    s.sliceLoc(i) = str2num(oCompTmp.pSliceLocation);
                    s.sliceNr(i) = i;
                    if isempty(oCompTmp.pUseSlice)
                        oCompTmp.pUseSlice = false;
                    end
                    s.useSlice (i) = oCompTmp.pUseSlice;
                end
                [keepInd status message] = LMranging(oComp, oCont.pAcfg.xlsSaveRange);
                if status==1
                    answer = questdlg(message, 'Slice indexing error', 'Save marked slices', 'Abord', 'Abord')
                    switch answer
                        case 'Abord'
                            return
                        case 'Save marked slices'
                            keepInd = s.useSlice;
                    end
                end
                
                tmp = s.useSlice;
                s.useSlice = zeros(size(tmp));
                s.useSlice(keepInd) = tmp(keepInd);
                % postprocessing
                s.Loc1(ismember(s.Loc1, 'none')) = {''};
                s.Loc2(ismember(s.Loc2, 'none')) = {''};
                
                sFlip = structfun(@(x) x', s, 'Uniformoutput', 0);
                    %% write to xls sheet (use all available data)
                writetable(table(hemiLinePos), xlsPath, 'Sheet', 'Pos');
                writetable(struct2table(info), xlsPath, 'Sheet', 'infos');
                writetable(struct2table(sFlip), xlsPath, 'Sheet', 'allFatData');
                waitbar(0.6, wb);
                    %% write to xls sheet (xls file like <nikita 2016)
                % filter accordint to pUseSlice
                s.SatVol(~s.useSlice) = 0;
                s.VatVol(~s.useSlice) = 0;
            if oCont.pAcfg.doSliceSpacingInterpolation
                if ~isnan(oCont.pAcfg.sliceSpacingInterpolationDistance) && any(abs(diff(s.sliceLoc) - oCont.pAcfg.sliceSpacingInterpolationDistance) > 0.5)
                    % prepare data: interpolate to equidistant slice loc
                    x = s.sliceLoc;
                    sliceLocOld = x;
                    window = oCont.pAcfg.sliceSpacingInterpolationDistance;
                    xn = [x(1):window:x(end)];
                    sliceLocNew = xn;
                    % find correct slice width (oContescription available (oContrawIO)) for volume calculation
                    sliceWidth = diff(sliceLocNew);
                    sW1 = sliceWidth./2; sW1(end+1) = 0;
                    sW2 = sliceWidth./2; sW2 = [0 sW2];
                    sliceWidth = sW1+sW2;
                    % now special treatment for boundary images (1 and end)
                    sliceWidth(1) = sliceWidth(1)+str2num(oCont.oComp(1).oImgs(1).sliceThickness)/2;
                    sliceWidth(end) = sliceWidth(end)+str2num(oCont.oComp(end).oImgs(end).sliceThickness)/2;
                    % slice Width done!
                    
                    %_Calc VatVol
                    y = s.VatArea;  % mm^2
                    [x2 y2] = DconvV2(x,y,window,'xmeanymean',[]); x2 = x2{1}; y2 = y2{1};
                    ynVatArea = interp1(x2, y2, sliceLocNew);
                    VatVol = ynVatArea.*sliceWidth*0.001;
                    
                    %_Calc SatVol
                    y = s.SatArea;  % mm^2
                    [x2 y2] = DconvV2(x,y,window,'xmeanymean',[]); x2 = x2{1}; y2 = y2{1};
                    ynSatArea = interp1(x2, y2, sliceLocNew);
                    SatVol = ynSatArea.*sliceWidth*0.001;
                    
                    %_Set WKs
                    LocsInd = find(~ismember(s.Loc1,''));
                    LocsVal = s.Loc1(LocsInd);
                    LocsPosOrig = s.sliceLoc(LocsInd);
                    Loc1 = cell(1, numel(sliceLocNew)); Loc1(:) = {''};
                    for i = 1:numel(LocsPosOrig)
                        [val ind] = min(abs(sliceLocNew-LocsPosOrig(i)));
                        Loc1(ind) = LocsVal(i);
                    end
                    
                    %_Set Landmarks
                    LocsInd = find(~ismember(s.Loc2,''));
                    LocsVal = s.Loc2(LocsInd);
                    LocsPosOrig = s.sliceLoc(LocsInd);
                    Loc2 = cell(1, numel(sliceLocNew)); Loc2(:) = {''};
                    for i = 1:numel(LocsPosOrig)
                        [val ind] = min(abs(sliceLocNew-LocsPosOrig(i)));
                        Loc2(ind) = LocsVal(i);
                    end
                    
                    %_Set sliceNr
                    sliceNr = 1:numel(sliceLocNew);
                    
                    %_Plot Interpolation Reults
                    figure();
                    plot(diff(x),'LineStyle', 'none', 'Marker', 'o', 'Color', 'black', 'DisplayName', 'raw');
                    hold on
                    plot(diff(x2),'LineStyle', 'none', 'Marker', 'x', 'Color', 'black', 'DisplayName', 'smooth');
                    plot(diff(xn),'LineStyle', 'none', 'Marker', '.', 'Color', 'red', 'DisplayName', 'interpolated');
                    drawnow;
                    a = gca;
                    a.XLabel.String = 'sliceLocation';
                    a.YLabel.String = 'VatArea';
                    legend('show');
                    
                    figure();
                    plot(x, y, 'LineStyle', 'none', 'Marker', 'o', 'Color', 'black', 'DisplayName', 'raw');
                    hold on
                    plot(x2,y2,'LineStyle', 'none', 'Marker', 'x', 'Color', 'black', 'DisplayName', 'smooth');
                    plot(xn,ynSatArea,'LineStyle', 'none', 'Marker', '.', 'Color', 'red', 'DisplayName', 'interpolated');
                    drawnow;
                    a = gca;
                    a.XLabel.String = 'sliceLocation';
                    a.YLabel.String = 'VatArea';
                    legend('show');
                    
                    
                    
                    
                    waitfor(msgbox('Slice locations are not equally distributed -> interpolation was done as shown in figure!'));
                    
                    
                else
                    sliceLocNew = s.sliceLoc;
                    SatVol = s.SatVol;
                    VatVol = s.VatVol;
                    sliceNr = s.sliceNr;
                    Loc1 = s.Loc1;
                    Loc2 = s.Loc2;
                end
            else
                sliceLocNew = s.sliceLoc;
                SatVol = s.SatVol;
                VatVol = s.VatVol;
                sliceNr = s.sliceNr;
                Loc1 = s.Loc1;
                Loc2 = s.Loc2;    
            end
                
                Pos = 1;
                xlswrite(xlsPath, {info.patientName} , 'FAT', 'A1');
                xlswrite(xlsPath, {info.creationDate} , 'FAT', 'A2');
                xlswrite(xlsPath, {'Slice #' 'Slice Pos (S)' 'Slice Gap' 'SAT [cm^3]' 'VAT [cm^3]' 'Loc1' 'Loc2'} , 'FAT', 'A3');
                xlswrite(xlsPath, sliceNr', 'FAT', 'A4');
                xlswrite(xlsPath, sliceLocNew', 'FAT', 'B4');
                xlswrite(xlsPath, [0 diff(sliceLocNew)]', 'FAT', 'C4');
                xlswrite(xlsPath, round(SatVol)', 'FAT', 'D4');
                xlswrite(xlsPath, round(VatVol)', 'FAT', 'E4');
                xlswrite(xlsPath, Loc1', 'FAT', 'F4');
                xlswrite(xlsPath, Loc2', 'FAT', 'G4');
                Pos = 3+numel(VatVol)+1;
                xlswrite(xlsPath, {'Slice #' 'Slice Pos (S)' 'Slice Gap' 'SAT [cm^3]' 'VAT [cm^3]' 'Loc1' 'Loc2'} , 'FAT', ['A' num2str(Pos)]);
                Pos = Pos+2;
                xlswrite(xlsPath, {'Summe'}, 'FAT', ['A' num2str(Pos)]);
                xlswrite(xlsPath, round(sum(SatVol)), 'FAT', ['D' num2str(Pos)]);
                xlswrite(xlsPath, round(sum(VatVol)), 'FAT', ['E' num2str(Pos)]);
                waitbar(0.9, wb);
                pause(0.3);
                close(wb);
        end
        
        % % Femur Box % %
        function oComp = mShowFemurBox(oComp, oCont)
            oCont.pAcfg.contour.showFemurBox = ~oCont.pAcfg.contour.showFemurBox;
            if oCont.pAcfg.contour.showFemurBox == 0
                uiwait(msgbox('Femur Box not visible'));
            elseif oCont.pAcfg.contour.showFemurBox == 1
                uiwait(msgbox('Femur Box now visible'));
            end
        end
        
        function oComp = mSetFemurBox(oComp, oCont)
            hBox = imrect(oCont.pHandles.imgAxis);
            pos = hBox.getPosition;
            %-------------------------------%
            Message = text(0, 0, '', 'Parent', oCont.pHandles.imgAxis, 'Color', 'white');
            Message.String = 'Double click on box to confirm!';
            Message.Position = [pos(1)+pos(3)/2-40, pos(2)+pos(4)/3-3];
            Message.FontSize = 14;
            Message.FontWeight = 'bold';
            T1 = text(0, 0, '', 'Parent', oCont.pHandles.imgAxis, 'Color', 'white');
            T1.FontSize = 14;
            T1.FontWeight = 'bold';
            T2 = text(0, 0, '', 'Parent', oCont.pHandles.imgAxis, 'Color', 'white');
            T2.FontSize = 14;
            T2.FontWeight = 'bold';
            oComp.mShowFemurBoxSize(hBox, T1, T2);
            addNewPositionCallback(hBox, @(varargout)oComp.mShowFemurBoxSize(hBox, T1, T2));
            
            wait(hBox);
            
            oComp = oComp.mDelFemurBox(oCont);
            
            datInd = oCont.pTable.row;
            oCompTmp = oCont.oComp(datInd);
            oCompTmp.pVarious.femurBox = createMask(hBox);
            oComp(datInd) = oCompTmp;
            
            oCont.pAcfg.contour.showFemurBox = 1;
        end
        
        function oComp = mDelFemurBox(oComp, oCont)
            femurBoxInd = find(arrayfun(@(x) isfield(x.pVarious, 'femurBox'), oCont.oComp));
            if ~isempty(femurBoxInd)
                for i =1:numel(femurBoxInd)
                    oComp(femurBoxInd(i)).pVarious = rmfield(oComp(femurBoxInd(i)).pVarious, 'femurBox');
                end
            end
            oCont.pAcfg.contour.showFemurBox = 0;
        end
        
        function mShowFemurBoxSize(oComp, hBox, T1, T2)
            pos = hBox.getPosition;
            coord = bwboundaries(createMask(hBox));
            pos = [min(coord{1}(:,2)) min(coord{1}(:,1)) max(coord{1}(:,2))-min(coord{1}(:,2)) max(coord{1}(:,1))-min(coord{1}(:,1))];
            
            %-------------------------------%
            T1.String = num2str(round(pos(3)*oComp(1).oImgs(1).dicomInfo.PixelSpacing(1)));
            T1.Position = [pos(1)+pos(3)/2-5, pos(2)-6];
            T2.String = num2str(round(pos(4)*oComp(1).oImgs(1).dicomInfo.PixelSpacing(2)));
            T2.Position = [pos(1)+pos(3)+2, pos(2)+pos(4)/2];
        end
        
        function oComp = mSaveBoxResults(oComp, oCont)
            answer = questdlg('What should be saved?', 'Choose data saving:', 'only XLS', 'XLS and SessionFile', 'Abord', 'Abord');
            wb = waitbar(0, 'saving femur-box results');
            wb.Name = 'saving....';
            if strcmp(answer, 'Abord')
                return
            end
            if strcmp(answer, 'XLS and SessionFile')
                % save session file
                path = fullfile(oCont.pAcfg.lastLoadPath, [oCont.mGetSaveFilePrefix '_data.mat']);
                oCont.oComp.mSave_oComp(path, oCont);
            end
            waitbar(0.2, wb);
            if strcmp(answer, 'only XLS') | strcmp(answer, 'XLS and SessionFile')
            %% collect infos
                dicomInfo = oComp(1).mGetStandardImg.dicomInfo;
                try info.comment = dicomInfo.StudyComments; end
                try info.description = dicomInfo.RequestedProcedureDescription; end
                try info.physicianName = dicomInfo.ReferringPhysicianName.FamilyName; end
                try info.institution = dicomInfo.InstitutionName; end
                try info.stationName = dicomInfo.StationName; end
                try info.manufacturer = dicomInfo.Manufacturer; end
                try info.manufacturerModelName = dicomInfo.ManufacturerModelName; end
                
                try info.patientName = [dicomInfo.PatientName.FamilyName '_' dicomInfo.PatientName.GivenName];
                catch
                    try info.patientName = dicomInfo.PatientName.FamilyName;
                    catch
                        info.patientName = 'NoName';
                    end
                end
                try info.patientWeight = num2str(dicomInfo.PatientWeight); end
                try info.patientAge = dicomInfo.PatientAge; end
                try info.patientSex = dicomInfo.PatientSex; end
                try info.patientBirthDat = dicomInfo.PatientBirthDate; end
                try info.patientIoCont = dicomInfo.PatientID; end
                
                try info.creationDate = datestr(datenum(dicomInfo.InstanceCreationDate, 'yyyymmdd'), 'dd.mm.yyyy'); end
                
                % remove empty entries
                emptyInd = structfun(@isempty, info);
                infoFields = fieldnames(info);
                for i = 1:numel(emptyInd)
                    if emptyInd(i)
                        info = rmfield(info, infoFields(i));
                    end
                end
                
                waitbar(0.3, wb);
            %% create struct with oComp!
                img = oComp(1).mGetStandardImg;
                imgSize = size(img.data);
                    %% preparation
                xlsPath = fullfile(oCont.pAcfg.lastLoadPath, [oCont.mGetSaveFilePrefix 'FemurBox_data.xlsx']);
                    %% get masks
                femurBoxInd = find(arrayfun(@(x) isfield(x.pVarious, 'femurBox'), oCont.oComp));
                mask = oCont.oComp(femurBoxInd).pVarious.femurBox;
                coord = bwboundaries(mask);
                femurBoxPos = [min(coord{1}(:,2)) min(coord{1}(:,1)) max(coord{1}(:,2))-min(coord{1}(:,2)) max(coord{1}(:,1))-min(coord{1}(:,1))];
                    %% read values from oComp to s struct
                    for i = 1:numel(oComp)
                        oCompTmp = oComp(i);
                        %% modify Boundaries for left side only
                        if oCompTmp.pSliceDone
                            try
                                oBoundInd = oCompTmp.oBoundaries.getBoundInd('outerBound');
                                oBound = oCompTmp.oBoundaries(oBoundInd);
                                oBoundMask = oCompTmp.mGetBoundMask(img.data, oBound.coord);
                                oBoundMask = oBoundMask&mask;
                                tmp = bwboundaries(oBoundMask);
                                oBound.coord = tmp{1};
                            catch
                                oBound = [];
                            end
                            
                            try
                                iBoundInd = oCompTmp.oBoundaries.getBoundInd('innerBound');
                                iBound = oCompTmp.oBoundaries(iBoundInd);
                                iBoundMask = oCompTmp.mGetBoundMask(img.data, iBound.coord);
                                iBoundMask = iBoundMask&mask;
                                tmp = bwboundaries(iBoundMask);
                                iBound.coord = tmp{1};
                            catch
                                iBound = [];
                            end
                            
                            try
                                vBoundInd = oCompTmp.oBoundaries.getBoundInd('visceralBound');
                                vBound = oCompTmp.oBoundaries(vBoundInd);
                                vBoundMask = oCompTmp.mGetBoundMask(img.data, vBound.coord);
                                vBoundMask = vBoundMask&mask;
                                tmp = bwboundaries(vBoundMask);
                                vBound.coord = tmp{1};
                            catch
                                vBound = [];
                            end
                            
                            oCompTmp.oBoundaries = [oBound iBound vBound];
                            
                        end
                        
                        s.voxVol(i) = oCompTmp.oImgs(1).getVoxelVolume;
                        s.fatThresh(i) = oCompTmp.pFatThresh;
                        try
                            oBoundCoord = oCompTmp.oBoundaries(oCompTmp.oBoundaries.getBoundInd('outerBound'));
                        catch
                            oBoundCoord = [];
                        end
                        try
                            oBoundMask = oCompTmp.mGetBoundMask(oCompTmp.mWaterImg.data, oBoundcoord.coord);
                        catch
                            oBoundMask = [];
                        end
                        tmp = regionprops(oBoundMask , 'Area', 'Perimeter');
                        if isempty(tmp)
                            tmp(1).Area = 0;
                            tmp(1).Perimeter = 0;
                        end
                        s.bodyArea(i) = tmp.Area;
                        s.bodyPerimeter(i) = tmp.Perimeter;
                        
                        s.SatArea(i) = oCompTmp.mVolumeSAT/oCompTmp.oImgs(1).dicomInfo.SpacingBetweenSlices*1000;
                        s.VatArea(i) = oCompTmp.mVolumeVAT/oCompTmp.oImgs(1).dicomInfo.SpacingBetweenSlices*1000;
                        s.SatVol(i) = oCompTmp.mVolumeSAT;
                        s.VatVol(i) = oCompTmp.mVolumeVAT;
                        s.Loc1(i) = {oCompTmp.pLoc1};
                        s.Loc2(i) = {oCompTmp.pLoc2};
                        s.sliceLoc(i) = str2num(oCompTmp.pSliceLocation);
                        s.sliceNr(i) = i;
                        if isempty(oCompTmp.pUseSlice)
                            oCompTmp.pUseSlice = false;
                        end
                        s.useSlice (i) = oCompTmp.pUseSlice;
                    end
                    [keepInd status message] = LMranging(oComp, oCont.pAcfg.xlsSaveRange);
                    if status==1
                        answer = questdlg(message, 'Slice indexing error', 'Save marked slices', 'Abord', 'Abord')
                        switch answer
                            case 'Abord'
                                return
                            case 'Save marked slices'
                                keepInd = s.useSlice;
                        end
                    end
                    
                    tmp = s.useSlice;
                    s.useSlice = zeros(size(tmp));
                    s.useSlice(keepInd) = tmp(keepInd);
                % postprocessing
                s.Loc1(ismember(s.Loc1, 'none')) = {''};
                s.Loc2(ismember(s.Loc2, 'none')) = {''};
                
                sFlip = structfun(@(x) x', s, 'Uniformoutput', 0);
                    %% write to xls sheet (use all available data)
                writetable(table(femurBoxPos), xlsPath, 'Sheet', 'Pos');
                writetable(struct2table(info), xlsPath, 'Sheet', 'infos');
                writetable(struct2table(sFlip), xlsPath, 'Sheet', 'allFatData');
                waitbar(0.6, wb);
                    %% write to xls sheet (xls file like <nikita 2016)
                % filter accordint to pUseSlice
                s.SatVol(~s.useSlice) = 0;
                s.VatVol(~s.useSlice) = 0;
            if oCont.pAcfg.doSliceSpacingInterpolation
                if ~isnan(oCont.pAcfg.sliceSpacingInterpolationDistance) && any(abs(diff(s.sliceLoc) - oCont.pAcfg.sliceSpacingInterpolationDistance) > 0.5)
                    % prepare data: interpolate to equidistant slice loc
                    x = s.sliceLoc;
                    sliceLocOld = x;
                    window = oCont.pAcfg.sliceSpacingInterpolationDistance;
                    xn = [x(1):window:x(end)];
                    sliceLocNew = xn;
                    % find correct slice width (oContescription available (oContrawIO)) for volume calculation
                    sliceWidth = diff(sliceLocNew);
                    sW1 = sliceWidth./2; sW1(end+1) = 0;
                    sW2 = sliceWidth./2; sW2 = [0 sW2];
                    sliceWidth = sW1+sW2;
                    % now special treatment for boundary images (1 and end)
                    sliceWidth(1) = sliceWidth(1)+str2num(oCont.oComp(1).oImgs(1).sliceThickness)/2;
                    sliceWidth(end) = sliceWidth(end)+str2num(oCont.oComp(end).oImgs(end).sliceThickness)/2;
                    % slice Width done!
                    
                    %_Calc VatVol
                    y = s.VatArea;  % mm^2
                    [x2 y2] = DconvV2(x,y,window,'xmeanymean',[]); x2 = x2{1}; y2 = y2{1};
                    ynVatArea = interp1(x2, y2, sliceLocNew);
                    VatVol = ynVatArea.*sliceWidth*0.001;
                    
                    %_Calc SatVol
                    y = s.SatArea;  % mm^2
                    [x2 y2] = DconvV2(x,y,window,'xmeanymean',[]); x2 = x2{1}; y2 = y2{1};
                    ynSatArea = interp1(x2, y2, sliceLocNew);
                    SatVol = ynSatArea.*sliceWidth*0.001;
                    
                    %_Set WKs
                    LocsInd = find(~ismember(s.Loc1,''));
                    LocsVal = s.Loc1(LocsInd);
                    LocsPosOrig = s.sliceLoc(LocsInd);
                    Loc1 = cell(1, numel(sliceLocNew)); Loc1(:) = {''};
                    for i = 1:numel(LocsPosOrig)
                        [val ind] = min(abs(sliceLocNew-LocsPosOrig(i)));
                        Loc1(ind) = LocsVal(i);
                    end
                    
                    %_Set Landmarks
                    LocsInd = find(~ismember(s.Loc2,''));
                    LocsVal = s.Loc2(LocsInd);
                    LocsPosOrig = s.sliceLoc(LocsInd);
                    Loc2 = cell(1, numel(sliceLocNew)); Loc2(:) = {''};
                    for i = 1:numel(LocsPosOrig)
                        [val ind] = min(abs(sliceLocNew-LocsPosOrig(i)));
                        Loc2(ind) = LocsVal(i);
                    end
                    
                    %_Set sliceNr
                    sliceNr = 1:numel(sliceLocNew);
                    
                    %_Plot Interpolation Reults
                    figure();
                    plot(diff(x),'LineStyle', 'none', 'Marker', 'o', 'Color', 'black', 'DisplayName', 'raw');
                    hold on
                    plot(diff(x2),'LineStyle', 'none', 'Marker', 'x', 'Color', 'black', 'DisplayName', 'smooth');
                    plot(diff(xn),'LineStyle', 'none', 'Marker', '.', 'Color', 'red', 'DisplayName', 'interpolated');
                    drawnow;
                    a = gca;
                    a.XLabel.String = 'sliceLocation';
                    a.YLabel.String = 'VatArea';
                    legend('show');
                    
                    figure();
                    plot(x, y, 'LineStyle', 'none', 'Marker', 'o', 'Color', 'black', 'DisplayName', 'raw');
                    hold on
                    plot(x2,y2,'LineStyle', 'none', 'Marker', 'x', 'Color', 'black', 'DisplayName', 'smooth');
                    plot(xn,ynSatArea,'LineStyle', 'none', 'Marker', '.', 'Color', 'red', 'DisplayName', 'interpolated');
                    drawnow;
                    a = gca;
                    a.XLabel.String = 'sliceLocation';
                    a.YLabel.String = 'VatArea';
                    legend('show');
                    
                    
                    
                    
                    waitfor(msgbox('Slice locations are not equally distributed -> interpolation was done as shown in figure!'));
                    
                    
                else
                    sliceLocNew = s.sliceLoc;
                    SatVol = s.SatVol;
                    VatVol = s.VatVol;
                    sliceNr = s.sliceNr;
                    Loc1 = s.Loc1;
                    Loc2 = s.Loc2;
                end
            else
                sliceLocNew = s.sliceLoc;
                SatVol = s.SatVol;
                VatVol = s.VatVol;
                sliceNr = s.sliceNr;
                Loc1 = s.Loc1;
                Loc2 = s.Loc2;    
            end
                
                Pos = 1;
                xlswrite(xlsPath, {info.patientName} , 'FAT', 'A1');
                xlswrite(xlsPath, {info.creationDate} , 'FAT', 'A2');
                xlswrite(xlsPath, {'Slice #' 'Slice Pos (S)' 'Slice Gap' 'SAT [cm^3]' 'VAT [cm^3]' 'Loc1' 'Loc2'} , 'FAT', 'A3');
                xlswrite(xlsPath, sliceNr', 'FAT', 'A4');
                xlswrite(xlsPath, sliceLocNew', 'FAT', 'B4');
                xlswrite(xlsPath, [0 diff(sliceLocNew)]', 'FAT', 'C4');
                xlswrite(xlsPath, round(SatVol)', 'FAT', 'D4');
                xlswrite(xlsPath, round(VatVol)', 'FAT', 'E4');
                xlswrite(xlsPath, Loc1', 'FAT', 'F4');
                xlswrite(xlsPath, Loc2', 'FAT', 'G4');
                Pos = 3+numel(VatVol)+1;
                xlswrite(xlsPath, {'Slice #' 'Slice Pos (S)' 'Slice Gap' 'SAT [cm^3]' 'VAT [cm^3]' 'Loc1' 'Loc2'} , 'FAT', ['A' num2str(Pos)]);
                Pos = Pos+2;
                xlswrite(xlsPath, {'Summe'}, 'FAT', ['A' num2str(Pos)]);
                xlswrite(xlsPath, round(sum(SatVol)), 'FAT', ['D' num2str(Pos)]);
                xlswrite(xlsPath, round(sum(VatVol)), 'FAT', ['E' num2str(Pos)]);
                waitbar(0.9, wb);
                pause(0.3);
                close(wb);
            end
        end
        
        % % other % %
        function oCont = mImportNikitaBound(oComp, oCont)
            % import only roi data from old Nikita tool < 2016
            answer = questdlg('klick OK to load an old Nikita mat file and overwrite the current segmentation!', 'load segmentation data', 'OK', 'Abord', 'Abord');
            switch answer
                case 'OK'
                    [file folder] = uigetfile(fullfile(oCont.pAcfg.lastLoadPath, '*.mat'));
                    load(fullfile(folder, file), 'threshold');
                    load(fullfile(folder, file), 'X_seg');
                    names = {'outerBound', 'innerBound', 'visceralBound'};
                    for i = 1:numel(names) %run through names
                        for j = 1:numel(X_seg(1,1,1,:,i))   %run through slices
                            cBoundCoords = bwboundaries(X_seg(:,:,1,j,i));
                            if isempty(cBoundCoords)
                            else
                                cBoundIn = oComp(j).oBoundaries.empty;
                                cBoundIn(1).name = names{i};
                                cBoundIn.coord = cBoundCoords{1};
                                
                                cBounds = oComp(j).oBoundaries;
                                cBounds = cBounds.setBound(cBoundIn);
                                oComp(j).oBoundaries = cBounds;
                                %% translate threshold
                                th = threshold(j)/255;
                                img = oComp(j).mWaterImg;
                                bInd = oComp(j).oBoundaries.getBoundInd('visceralBound');
                                if ~bInd==0
                                    bMask = oComp(j).oBoundaries(bInd).getBoundMask(img);
                                    oComp(j).pFatThresh = max(img.data(find(bMask)))*th;
                                end
                                %% check if slice is done
                                oComp(j) = oComp(j).mSet_pSliceDone;
                            end
                            
                        end
                        
                    end
                    
                    oCont.oComp = oComp;
                    oCont.mTableCellSelect('Caller', 'mImportNikitaBound');
                case 'Abord'
                    return
            end
        end
        
        % % % Image Organisation % % %
        function img = mFatImg(oComp)
            img = cImageDcm;
            inPhaseImg = oComp.mGetImgOfType('InPhase');
            outPhaseImg = oComp.mGetImgOfType('OutPhase');
            
            img.data = (inPhaseImg.data-outPhaseImg.data)./2;
            img.imgType = 'Fat';
        end
        
        function img = mWaterImg(oComp)
            img = cImageDcm;
            inPhaseImg = oComp.mGetImgOfType('InPhase');
            outPhaseImg = oComp.mGetImgOfType('OutPhase');
            
            img.data = (inPhaseImg.data+outPhaseImg.data)./2;
            img.imgType = 'Water';
        end
        
        function mVatOverlay(oComp, oCont)
            vbInd = oComp.oBoundaries.getBoundInd('visceralBound');
            VatImg = oComp.mWaterImg;
            VatMask = ~logical(oComp.mGetBoundMask(VatImg.data, oComp.oBoundaries(vbInd).coord));
            %-------------------------------%
            VatImg.data(VatMask) = 0;   % only Vat pxls are not 0
            VatFatPxl = VatImg.data>oComp.pFatThresh;  % Vat is only pxlValues bigger than hist cursor
            % indicate the VatFat pixels with overlay
            maskInd = find(VatFatPxl);
            color = 1-oCont.pAcfg.color3;
            %oComp.pHistory.plotImgRGB = oComp.pHistory.plotImg.conv2RGB;
            img = oComp.mGetImg2Display(oCont);
            img = img.data;
%             if isfield(oComp.pHistory, 'vatOverlay') && ~isempty(oComp.pHistory.vatOverlay)
%                 img = oCont.pHandles.imgDisplay.CData;
%                 % remove overlay
%                 oldInd = oComp.pHistory.vatOverlay;
%                 Rimg = img(:,:,1);
%                 Rimg(oldInd) = Rimg(oldInd)./(1-color(1));
%                 Gimg = img(:,:,2);
%                 Gimg(oldInd) = Gimg(oldInd)./(1-color(1));
%                 Bimg = img(:,:,3);
%                 Bimg(oldInd) = Bimg(oldInd)./(1-color(1));
%             else
%                 img = oCont.pHandles.imgDisplay.CData;
%             end
            
            Rimg = img(:,:,1);
            Rimg(maskInd) = Rimg(maskInd) - Rimg(maskInd).*color(1);
            Gimg = img(:,:,1);
            Gimg(maskInd) = Gimg(maskInd) - Gimg(maskInd).*color(2);
            Bimg = img(:,:,1);
            Bimg(maskInd) = Bimg(maskInd) - Bimg(maskInd).*color(3);
            
            RGBimg = [Rimg Gimg Bimg];
            RGBimg = reshape(RGBimg, size(Rimg,1), size(Rimg,2), 3);
            
            %-------------------------------%
            oCont.pHandles.imgDisplay.CData = uint8(RGBimg);
        end
        
        % % % Segmentaion Organisation % % %
        function oComp = mVisceralFromInnerBound(oComp, oCont)
            if numel(oComp)==1
                Ind = 1;
            else
                Ind = oCont.pTable.row;
            end
            iB = oComp(Ind).oBoundaries.getBoundOfType('innerBound');
            vB = oComp(Ind).oBoundaries.getBoundOfType('visceralBound');
            
            imgTmp = oComp(Ind).mGetImgOfType('OutPhase');
            iBMask = oComp(Ind).mGetBoundMask(imgTmp.data, iB.coord);
            
            %% % find visceralBound
            visceralBound.Mask = imerode(iBMask, strel('disk', 2, 0));
            cc = bwconncomp(visceralBound.Mask,4);
            [a maxInd] = max(cellfun(@(x) size(x,1), cc.PixelIdxList));    % find maximum area
            visceralBound.Mask = zeros(size(iBMask));
            visceralBound.Mask(cc.PixelIdxList{maxInd}) = 1;
            
            visceralBound.coord = cCompute.mGetBoundaryImageCoord(visceralBound.Mask);
            
            vB.coord = visceralBound.coord;
            
            oComp(Ind).oBoundaries = oComp(Ind).oBoundaries.setBound(vB);
            
            if numel(oComp(Ind).oBoundaries)==3
                if oComp(Ind).pFatThresh==0
                    oComp(Ind) = oComp(Ind).mFindThreshLvl(oCont);
                end
            end
            oComp(Ind) = oComp(Ind).mSet_pSliceDone;
        end
        
        function oComp = mAutoSegment(oComp, segProps)
            iBound = oComp.oBoundaries.getBoundOfType('innerBound');
            oBound = oComp.oBoundaries.getBoundOfType('outerBound');
            vBound = oComp.oBoundaries.getBoundOfType('visceralBound');
            %-------------------------------%
            switch segProps.name %segProps.name
                case 'NikitaFat160322Segmenting.m'
                    imgBase = oComp.mGetImgOfType('OutPhase');
                    imgInfo = imgBase.dicomInfo;
                    imgData = imgBase.data;
                    
                    [outerBound innerBound visceralBound] = NikitaFat160322Segmenting(imgData, imgInfo);
                    outerBound = logical(outerBound);
                    innerBound = logical(innerBound);
                    visceralBound = logical(visceralBound);
                    
                    oBound.coord = bwboundaries(outerBound);
                    iBound.coord = bwboundaries(innerBound);
                    vBound.coord = bwboundaries(visceralBound);
                case 'RS_BodyBounds'
                    imgOut = oComp.mGetImgOfType('OutPhase');
                    iDat = imgOut.data;
                    iZero = zeros(size(iDat));
                    
                    oB.ThreshMultip = 1;
                    oB.ThreshWindowSizeMultip = 1;
                    
                    debug = false;
                    if debug
                        f = figure;
                        im = imshow(iDat,'DisplayRange', []);
                        ax = f.Children;
                    end
                    tic
                    code = 'start';
                    while ~strcmp(code, 'Done')
                        disp(code)
                        switch code
                            % RS_BodyBounds(img, mode, outerBoundThreshMultip, innerBoundThreshMultip, innerBoundThreshWindowSizeMultip, visceralBoundErodeSize)
                            case 'start'
                                oB.status = 'to be done';
                                code = 'AutoSegment_oB';
                            case 'AutoSegment_oB'
                                iB = [];
                                vB = [];
                                [oB iB vB code] = RS_BodyBounds(imgOut, 'MRT_OutPhase', 'outerBound', oB, iB, vB);
                            case 'Precheck: outerBound in dark area'
                                oB.ThreshWindowSizeMultip = oB.ThreshWindowSizeMultip+0.05;
                                code = 'AutoSegment_oB';
                            case 'Precheck: outerBound not round'
                                oB.ThreshMultip = oB.ThreshMultip-0.05;
                                code = 'AutoSegment_oB';
                            case 'outerBound OK'
                                code = 'AutoSegment_iB_1';
                                %code = 'checkBoundaryAppropriateness';
                                
                            case 'AutoSegment_iB_1'
                                iB.ThreshMultip = 0.7;
                                iB.ThreshWindowSizeMultip = 0.9;
                                iB.ThreshSensitivity = 0.1;
                                iB.maxCircularity = 2.3;
                                vB = [];
                                [oB iB vB code] = RS_BodyBounds(imgOut, 'MRT_OutPhase', 'innerBound', oB, iB, vB);
                                iBs(1) = iB;
                                code = 'AutoSegment_iB_2';
                            case 'AutoSegment_iB_2'
                                iB.ThreshMultip = 0.7;
                                iB.ThreshWindowSizeMultip = 4;
                                iB.ThreshSensitivity = 0.1;
                                iB.maxCircularity = 2.3;
                                vB = [];
                                [oB iB vB code] = RS_BodyBounds(imgOut, 'MRT_OutPhase', 'innerBound', oB, iB, vB);
                                iBs(2) = iB;
                                code = 'AutoSegment_iB_3';
                            case 'AutoSegment_iB_3'
                                iB.ThreshMultip = 0.8;
                                iB.ThreshWindowSizeMultip = 0.7;
                                iB.ThreshSensitivity = 0.6;
                                iB.maxCircularity = 2.3;
                                vB = [];
                                [oB iB vB code] = RS_BodyBounds(imgOut, 'MRT_OutPhase', 'innerBound', oB, iB, vB);
                                iBs(3) = iB;
                                code = 'AutoSegment_iB_4';
                            case 'AutoSegment_iB_4'
                                iB.ThreshMultip = 0.72;
                                iB.ThreshWindowSizeMultip = 0.4;
                                iB.ThreshSensitivity = 0.1;
                                iB.maxCircularity = 2.3;
                                vB = [];
                                [oB iB vB code] = RS_BodyBounds(imgOut, 'MRT_OutPhase', 'innerBound', oB, iB, vB);
                                iBs(4) = iB;
                                code = 'AutoSegment_iB_5';
                            case 'AutoSegment_iB_5'
                                iB.ThreshMultip = 1.1;
                                iB.ThreshWindowSizeMultip = 0.7;
                                iB.ThreshSensitivity = 0.1;
                                iB.maxCircularity = 2.3;
                                vB = [];
                                [oB iB vB code] = RS_BodyBounds(imgOut, 'MRT_OutPhase', 'innerBound', oB, iB, vB);
                                iBs(5) = iB;
                                code = 'AutoSegment_iB_6';
                            case 'AutoSegment_iB_6'
                                iB.ThreshMultip = 0.8;
                                iB.ThreshWindowSizeMultip = 0.65;
                                iB.ThreshSensitivity = 0.4;
                                iB.maxCircularity = 2.3;
                                vB = [];
                                [oB iB vB code] = RS_BodyBounds(imgOut, 'MRT_OutPhase', 'innerBound', oB, iB, vB);
                                iBs(6) = iB;
                                code = 'selectBestBound_iB';
                                
                            case 'selectBestBound_iB'
                                [ind quality] = checkBoundaryAppropriateness('innerBound', oB, iBs , vB, iDat, iZero);
                                %ind=1
                                oComp.pVarious.AutoSegment_ind = ind;
                                oComp.pVarious.AutoSegment_quality = quality;
                                iB = iBs(ind);
                                code = 'AutoSegment_vB_1'
                                
                            case 'AutoSegment_vB_1'
                                vB.ErodeSize = 2;
                                [oB iB vB code] = RS_BodyBounds(imgOut, 'MRT_OutPhase', 'visceralBound', oB, iB, vB);
                                vBs(1) = vB;
                                code = 'selectBestBound_vB';
                            case 'selectBestBound_vB'
                                vB = vB;
                                code = 'Done'
                                
                            
                        end
                        
                        if debug
                            axChilds = ax.Children;
                            for i = 1:numel(axChilds)
                                if isa(axChilds(i), 'matlab.graphics.primitive.Image')
                                else
                                    delete(axChilds(i));
                                end
                            end
                            try
                                cControl.mDrawContour(ax, cCompute.mGetBoundaryImageCoord(oB.Mask), {[1 1 1]});
                                cControl.mDrawContour(ax, cCompute.mGetBoundaryImageCoord(iB.Mask), {[1 0 0]}, 'drawMode', {'dot'});
                                cControl.mDrawContour(ax, cCompute.mGetBoundaryImageCoord(vB.Mask), {[0 1 0]}, 'drawMode', {'dot'});
                            end
                            pause
                        end
                    end
                    
                    toc
                    try
                        oBound.coord = {oB.coord};
                    end
                    try
                        iBound = iBound(1);
                        iBound.coord = {iB.coord};
                    end
                    try
                        vBound.coord = {vB.coord};
                    end
                case 'OuterBound_RS'
                    imgIn = oComp.mGetImgOfType('InPhase');
                    imgOut = oComp.mGetImgOfType('OutPhase');
                    oBound.coord = {RSouterBound(imgIn, 'MRT_OutPhase', segProps.magThreshold)};
                case 'OuterBound_RS+InnerBound_RS'
                    imgIn = oComp.mGetImgOfType('InPhase');
                    imgOut = oComp.mGetImgOfType('OutPhase');
                    segProps.bwThreshold = 2.3;
                    segProps.erodeLvl = 2;
                    oBound.coord = {RSouterBound(imgIn, 'MRT_OutPhase', segProps.magThreshold)};
                    iBound.coord = {RSinnerBound(imgOut, 'MRT_OutPhase', oBound.coord, segProps.bwThreshold)};
                    %vBound.coord = {RSvisceralBound(imgOut, 'MRT_OutPhase', iBound.coord, segProps.erodeLvl)};
                case 'OuterBound_RS+NikitaFat'
                    imgIn = oComp.mGetImgOfType('InPhase');
                    imgOut = oComp.mGetImgOfType('OutPhase');
                    outerBoundCoord = {RSouterBound(imgIn, 'MRT_OutPhase', 8)};
                    
                    imgBody = imgOut.data;
                    indBody = find(~oComp.mGetBoundMask(imgOut.data, outerBoundCoord));
                    imgBody(indBody) = 0;
                    
                    
                    [outerBound innerBound visceralBound] = NikitaFat161213Segmenting(imgBody, imgOut.dicomInfo);
                    %visceralBound = imerode(innerBound, strel('diamond', 1));
                    
                    outerBound = logical(outerBound);
                    innerBound = logical(innerBound);
                    visceralBound = logical(visceralBound);
                    
                    oBound.coord = outerBoundCoord;
                    iBound.coord = bwboundaries(innerBound);
                    vBound.coord = bwboundaries(visceralBound);
                case 'NikitaFat161213Segmenting.m'
                    imgBase = oComp.mGetImgOfType('OutPhase');
                    imgInfo = imgBase.dicomInfo;
                    imgData = imgBase.data;
                    
                    [outerBound innerBound visceralBound] = NikitaFat161213Segmenting(imgData, imgInfo);
                    outerBound = logical(outerBound);
                    innerBound = logical(innerBound);
                    visceralBound = logical(visceralBound);
                    
                    oBound.coord = bwboundaries(outerBound);
                    iBound.coord = bwboundaries(innerBound);
                    vBound.coord = bwboundaries(visceralBound);
            end
            %-------------------------------%
            try
                oBound.coord = oBound.coord{1};
                oComp.oBoundaries = oComp.oBoundaries.setBound(oBound);
            end
            try
                iBound.coord = iBound.coord{1};
                oComp.oBoundaries = oComp.oBoundaries.setBound(iBound);
            end
            try
                vBound.coord = vBound.coord{1};
                oComp.oBoundaries = oComp.oBoundaries.setBound(vBound);
            end
            oComp = oComp.mSet_pSliceDone;
        end
        
        function oComp = mAutoSegmentSingle(oComp, oCont)
            if numel(oComp)==1
                Ind = 1;
            else
                Ind = oCont.pTable.row;
            end
            oComp(Ind) = oComp(Ind).mAutoSegment(oCont.pAcfg.segProps);
            
            if numel(oComp(Ind).oBoundaries)==3
                if oComp(Ind).pFatThresh<1
                    oComp(Ind) = oComp(Ind).mFindThreshLvl(oCont);
                end
            end
            
%             if oComp(Ind).pVarious.AutoSegment_quality>0.75
%                 oComp(Ind) = oComp(Ind).mSet_pSliceDone;
%             end
            oComp(Ind) = oComp(Ind).mSet_pSliceDone;
            
            %-------------------------------%
        end
        
        function oComp = mAutoSegmentAll(oComp, oCont)
            imgInd = LMranging(oComp, oCont.pAcfg.xlsSaveRange);
            
            wb = waitbar(0, ['Image 1 of ' num2str(numel(imgInd)) '. Time left: inf']);
            wb.Name = 'segmenting image data';
            n = 0;
            for i=imgInd
                n=n+1
                if ~oComp(i).pSliceDone
                    t1 = tic;
                    oComp(i) = oComp(i).mAutoSegment(oCont.pAcfg.segProps);
                    
                    if numel(oComp(i).oBoundaries)==3
                        if oComp(i).pFatThresh<1
                            oComp(i) = oComp(i).mFindThreshLvl(oCont);
                        end
                    end
                    oComp(i) = oComp(i).mSet_pSliceDone;
                    t2 = toc(t1);
                    time(n) = t2;
                    timeLeft = median(time)*(numel(imgInd)-i);
                    if ishandle(wb)
                        waitbar(n/numel(imgInd), wb, ['Image ' num2str(n+1) ' of ' num2str(numel(imgInd)) '. Time left: ' num2str(timeLeft, '%.0f') ' sec']);
                    else
                        disp('aha');
                        return % user abord
                    end
                end
            end
            close(wb);
        end
        
        function oComp = mSet_pSliceDone(oComp)
            if ~oComp.oBoundaries.getBoundInd('outerBound')==0 & ~oComp.oBoundaries.getBoundInd('innerBound')==0 & ~oComp.oBoundaries.getBoundInd('visceralBound')==0 & oComp.pFatThresh>1
                % use slice if it switches to Done from notDone
                if oComp.pSliceDone == false
                    oComp.pUseSlice = true;
                end
                oComp.pSliceDone = true;
            else
                oComp.pSliceDone = false;
            end
            
        end
        
            function oComp = mClearSegment(oComp)
            oComp
                % to be done
            
        end
        
        function oComp = mFindThreshLvl(oComp, oCont)
            if numel(oComp)==1
                Ind = 1;
            else
                Ind = oCont.pTable.row;
            end
            oCompTmp = oComp(Ind);
            %-------------------------------%
            if ~oCompTmp.oBoundaries.getBoundInd('outerBound')==0 & ~oCompTmp.oBoundaries.getBoundInd('innerBound')==0
                %AUTOMATISCH THRESHOLDERMITTLUNG
                outerCoord = oCompTmp.oBoundaries(oCompTmp.oBoundaries.getBoundInd('outerBound')).coord;
                innerCoord = oCompTmp.oBoundaries(oCompTmp.oBoundaries.getBoundInd('innerBound')).coord;
                visceralCoord = oCompTmp.oBoundaries(oCompTmp.oBoundaries.getBoundInd('visceralBound')).coord;
                
                outerMask = oCompTmp.mGetBoundMask(oCompTmp.mWaterImg.data, outerCoord);
                innerMask = oCompTmp.mGetBoundMask(oCompTmp.mWaterImg.data, innerCoord);
                satMask = outerMask&~innerMask;
                
                satMask(:,1:min(visceralCoord(:,2))) = 0;
                satMask(:,max(visceralCoord(:,2)):end) = 0;
                satMask = imerode(satMask, strel('diamond', 6));
                
                
                satIntensities = oCompTmp.mWaterImg.data(find(satMask));
                
                %                 figure();
                %                 imshow(satMask);
                %                 imshow(imerode(satMask, strel('diamond',5)));
                %                 histogram(satIntensities, 100);
                
                %oCompTmp.pFatThresh = prctile(satIntensities, 0.5);
                oCompTmp.pFatThresh = median(satIntensities)*0.7;
                
            else
                uiwait(msgbox('SAT must be segmented for FatThreshold determination!'));
            end
            
            oComp(Ind) = oCompTmp;
            %oCont.mTableCellSelect('Caller', 'mFindThreshLvl');
            
            %
            %             oCont.pAcfg.autoThresh = ~oCont.pAcfg.autoThresh;
            %             if oCont.pAcfg.autoThresh == 0
            %                 uiwait(msgbox('AutoThresholding is Off now! Please set threshold level manually.'));
            %             elseif oCont.pAcfg.autoThresh == 1
            %                 uiwait(msgbox('AutoThresholding is On now!'));
%                 oComp = oComp.determineFatThresh(oCont);
%             end
        end
        
        function volumeSAT = mVolumeSAT(oComp) % in cm^3
            try
                obInd = oComp.oBoundaries.getBoundInd('outerBound');
            catch
                uiwait(msgbox(['OuterBoundary could not be found for slice location ' num2str(oComp.pSliceLocation)]));
                volumeSAT = 0;
                return
            end
            try
                ibInd = oComp.oBoundaries.getBoundInd('innerBound');
            catch
                uiwait(msgbox(['InnerBoundary could not be found for slice location ' num2str(oComp.pSliceLocation)]));
                volumeSAT = 0;
                return
            end
            %-------------------------------%
            if obInd==0 | ibInd==0
                volumeSAT = 0;
            else
                pxlcount = (oComp.oBoundaries(obInd).areaBound - oComp.oBoundaries(ibInd).areaBound);
                volumeSAT = pxlcount*oComp.oImgs(1).getVoxelVolume;
            end
        end
        
        function volumeVAT = mVolumeVAT(oComp) % in cm^3
            try
                vbInd = oComp.oBoundaries.getBoundInd('visceralBound');
            catch
                uiwait(msgbox(['VisceralBoundary could not be found for slice location ' num2str(oComp.pSliceLocation)]));
                volumeVAT = 0;
                return
            end
            %-------------------------------%
            if vbInd==0
                volumeVAT = 0;
            else
                VatFatPxl = oComp.mGetVatFatPxl;  % Vat is only pxlValues bigger than hist cursor
                
                volumeVAT = sum(sum(VatFatPxl))*oComp.oImgs(1).getVoxelVolume;
            end
        end
        
        function VatFatPxl = mGetVatFatPxl(oComp)
            vbInd = oComp.oBoundaries.getBoundInd('visceralBound');
            VatImg = oComp.mWaterImg;
            VatMask = ~logical(oComp.mGetBoundMask(VatImg.data, oComp.oBoundaries(vbInd).coord, 'fillHoles', true));
            %-------------------------------%
            VatImg.data(VatMask) = 0;   % only Vat pxls are not 0
            VatFatPxl = VatImg.data>oComp.pFatThresh;  % Vat is only pxlValues bigger than hist cursor
        end
        
        function DiffSlice = mGetDiffSlices(oComp,sliceLoc,i)
            
%             maxSizeoComp = size(oCont.oComp,2);
%             
%             SlicLocTmp = cell2mat(cellfun(@(x) str2double(x), {oCont.oComp(:).pSliceLocation}, 'un', 0));
%             DiffSliceTmp = diff(SlicLocTmp);
            
            switch i 
                case 1
                    DiffSlice = 0;
                
                otherwise
                    DiffSlice = sliceLoc(end-1)-sliceLoc(end);
            end
            
            
        end
                
        % % % Application Management % % %
        function mSaveXls(oComp, oCont)
            %% preparation
            xlsPath = fullfile(oCont.pAcfg.lastLoadPath, [oCont.mGetSaveFilePrefix '_data.xlsx']);
            [file, path] = uiputfile(xlsPath);
            if path==0
                return
            end
            xlsPath = fullfile(path,file);
            %% collect infos and store in variables
            dicomInfo = oComp(1).mGetStandardImg.dicomInfo;
            try info.comment = dicomInfo.StudyComments; end
            try info.description = dicomInfo.RequestedProcedureDescription; end
            try info.physicianName = dicomInfo.ReferringPhysicianName.FamilyName; end
            try info.institution = dicomInfo.InstitutionName; end
            try info.stationName = dicomInfo.StationName; end
            try info.manufacturer = dicomInfo.Manufacturer; end
            try info.manufacturerModelName = dicomInfo.ManufacturerModelName; end
            
            try info.patientName = [dicomInfo.PatientName.FamilyName '_' dicomInfo.PatientName.GivenName]; 
            catch
                try info.patientName = dicomInfo.PatientName.FamilyName;
                catch
                    info.patientName = 'NoName';
                end
            end
            try info.patientWeight = num2str(dicomInfo.PatientWeight); end
            try info.patientAge = dicomInfo.PatientAge; end
            try info.patientSex = dicomInfo.PatientSex; end
            try info.patientBirthDat = dicomInfo.PatientBirthDate; end
            try info.patientIoCont = dicomInfo.PatientID; end
            
            try info.creationDate = datestr(datenum(dicomInfo.InstanceCreationDate, 'yyyymmdd'), 'dd.mm.yyyy'); end
            
            % remove empty entries
            emptyInd = structfun(@isempty, info);
            infoFields = fieldnames(info);
            for i = 1:numel(emptyInd)
                if emptyInd(i)
                    info = rmfield(info, infoFields(i));
                end
            end
            
            % create struct with oComp data
            for i = 1:numel(oComp)
                oCompTmp = oComp(i);
                s.voxVol(i) = oCompTmp.oImgs(1).getVoxelVolume;
                s.fatThresh(i) = oCompTmp.pFatThresh;
                try 
                    oBoundCoord = oCompTmp.oBoundaries(oComp.oBoundaries.getBoundInd('outerBound')); 
                catch
                    oBoundCoord = []; 
                end
                try 
                    oBoundMask = oCompTmp.mGetBoundMask(oComp.mWaterImg.data, oBoundcoord.coord); 
                catch
                    oBoundMask = []; 
                end
                tmp = regionprops(oBoundMask , 'Area', 'Perimeter');
                if isempty(tmp)
                    tmp(1).Area = 0;
                    tmp(1).Perimeter = 0;
                end
                s.bodyArea(i) = tmp.Area;
                s.bodyPerimeter(i) = tmp.Perimeter;
                
                s.SatArea(i) = oCompTmp.mVolumeSAT/oCompTmp.oImgs(1).dicomInfo.SpacingBetweenSlices*1000;
                s.VatArea(i) = oCompTmp.mVolumeVAT/oCompTmp.oImgs(1).dicomInfo.SpacingBetweenSlices*1000;
                s.SatVol(i) = oCompTmp.mVolumeSAT;
                s.VatVol(i) = oCompTmp.mVolumeVAT;
                s.Loc1(i) = {oCompTmp.pLoc1};
                s.Loc2(i) = {oCompTmp.pLoc2};
                s.sliceLoc(i) = str2num(oCompTmp.pSliceLocation);
                s.sliceDiff(i) = oCompTmp.mGetDiffSlices(s.sliceLoc,i);
                s.sliceNr(i) = i;
                
                if isempty(oCompTmp.pUseSlice)
                        oCompTmp.pUseSlice = false;
                end
                s.useSlice (i) = oCompTmp.pUseSlice;
                
                
            end
            [keepInd status message] = LMranging(oComp, oCont.pAcfg.xlsSaveRange);
            if status==1
                answer = questdlg(message, 'Slice indexing error', 'Save marked slices', 'Abord', 'Abord')
                switch answer
                    case 'Abord'
                        return
                    case 'Save marked slices'
                        keepInd = s.useSlice;
                end
            end
                    
            tmp = s.useSlice;
            s.useSlice = zeros(size(tmp));
            s.useSlice(keepInd) = tmp(keepInd);
            
            % postprocessing
            s.Loc1(ismember(s.Loc1, 'none')) = {''};
            s.Loc2(ismember(s.Loc2, 'none')) = {''};
            
            sFlip = structfun(@(x) x', s, 'Uniformoutput', 0);
            
            %% write to xls sheet (use all available data)
            writetable(struct2table(info), xlsPath, 'Sheet', 'infos');
            writetable(struct2table(sFlip), xlsPath, 'Sheet', 'allFatData');
            
            %% write to xls sheet (xls file like <nikita 2016)
            % filter accordint to pUseSlice
            s.SatVol(~s.useSlice) = 0;
            s.VatVol(~s.useSlice) = 0;
            if oCont.pAcfg.doSliceSpacingInterpolation
            inconsistentSliceSpacing = any(abs(diff(s.sliceLoc) - oCont.pAcfg.sliceSpacingInterpolationDistance) > 0.1);
            if inconsistentSliceSpacing
                if isnan(oCont.pAcfg.sliceSpacingInterpolationDistance)
                    sliceSpacingInterpolationDistance = str2num(inputdlg('Set a interpolation value for slice spacing in mm:', 'SliceSpacing inconsistent'));
                    oCont.pAcfg.sliceSpacingInterpolationDistance = sliceSpacingInterpolationDistance;
                end
                % prepare data: interpolate to equidistant slice loc
                x = s.sliceLoc;
                sliceLocOld = x;
                window = oCont.pAcfg.sliceSpacingInterpolationDistance;
                xn = [x(1):window:x(end)];
                sliceLocNew = xn;
                % find correct slice width (oContescription available (oContrawIO)) for volume calculation
                sliceWidth = diff(sliceLocNew);
                sW1 = sliceWidth./2; sW1(end+1) = 0;
                sW2 = sliceWidth./2; sW2 = [0 sW2];
                sliceWidth = sW1+sW2;
                % now special treatment for boundary images (1 and end)
                sliceWidth(1) = sliceWidth(1)+str2num(oCont.oComp(1).oImgs(1).sliceThickness)/2;
                sliceWidth(end) = sliceWidth(end)+str2num(oCont.oComp(end).oImgs(end).sliceThickness)/2;
                % slice Width done!
                
                %_Calc VatVol
                y = s.VatArea;  % mm^2
                [x2 y2] = DconvV2(x,y,window,'xmeanymean',[]); x2 = x2{1}; y2 = y2{1};
                ynVatArea = interp1(x2, y2, sliceLocNew);
                VatVol = ynVatArea.*sliceWidth*0.001;
                
                %_Calc SatVol
                y = s.SatArea;  % mm^2
                [x2 y2] = DconvV2(x,y,window,'xmeanymean',[]); x2 = x2{1}; y2 = y2{1};
                ynSatArea = interp1(x2, y2, sliceLocNew);
                SatVol = ynSatArea.*sliceWidth*0.001;
                
                %_Set WKs
                LocsInd = find(~ismember(s.Loc1,''));
                LocsVal = s.Loc1(LocsInd);
                LocsPosOrig = s.sliceLoc(LocsInd);
                Loc1 = cell(1, numel(sliceLocNew)); Loc1(:) = {''};
                for i = 1:numel(LocsPosOrig)
                    [val ind] = min(abs(sliceLocNew-LocsPosOrig(i)));
                    Loc1(ind) = LocsVal(i);
                end
                
                %_Set Landmarks
                LocsInd = find(~ismember(s.Loc2,''));
                LocsVal = s.Loc2(LocsInd);
                LocsPosOrig = s.sliceLoc(LocsInd);
                Loc2 = cell(1, numel(sliceLocNew)); Loc2(:) = {''};
                for i = 1:numel(LocsPosOrig)
                    [val ind] = min(abs(sliceLocNew-LocsPosOrig(i)));
                    Loc2(ind) = LocsVal(i);
                end
                
                %_Set sliceNr
                sliceNr = 1:numel(sliceLocNew);
                
                %_Plot Interpolation Reults
                figure();
                plot(diff(x),'LineStyle', 'none', 'Marker', 'o', 'Color', 'black', 'DisplayName', 'raw');
                hold on
                plot(diff(x2),'LineStyle', 'none', 'Marker', 'x', 'Color', 'black', 'DisplayName', 'smooth');
                plot(diff(xn),'LineStyle', 'none', 'Marker', '.', 'Color', 'red', 'DisplayName', 'interpolated');
                drawnow;
                a = gca;
                a.XLabel.String = 'sliceLocation';
                a.YLabel.String = 'VatArea';
                legend('show');
                
                figure();
                plot(x, y, 'LineStyle', 'none', 'Marker', 'o', 'Color', 'black', 'DisplayName', 'raw');
                hold on
                plot(x2,y2,'LineStyle', 'none', 'Marker', 'x', 'Color', 'black', 'DisplayName', 'smooth');
                plot(xn,ynSatArea,'LineStyle', 'none', 'Marker', '.', 'Color', 'red', 'DisplayName', 'interpolated');
                drawnow;
                a = gca;
                a.XLabel.String = 'sliceLocation';
                a.YLabel.String = 'VatArea';
                legend('show');
                
                
                
                
                waitfor(msgbox('Slice locations are not equally distributed -> interpolation was done as shown in figure!'));
                
                else
                sliceLocNew = s.sliceLoc;
                SatVol = s.SatVol;
                VatVol = s.VatVol;
                sliceNr = s.sliceNr;
                Loc1 = s.Loc1;
                Loc2 = s.Loc2;
            
            end
            else
                sliceLocNew = s.sliceLoc;
                SatVol = s.SatVol;
                VatVol = s.VatVol;
                sliceNr = s.sliceNr;
                Loc1 = s.Loc1;
                Loc2 = s.Loc2;
            end

            Pos = 1;
            xlswrite(xlsPath, {info.patientName} , 'FAT', 'A1');
            xlswrite(xlsPath, {info.creationDate} , 'FAT', 'A2');
            xlswrite(xlsPath, {'Slice #' 'Slice Pos (S)' 'Slice Gap' 'SAT [cm^3]' 'VAT [cm^3]' 'Loc1' 'Loc2'} , 'FAT', 'A3');
            xlswrite(xlsPath, sliceNr', 'FAT', 'A4');
            xlswrite(xlsPath, sliceLocNew', 'FAT', 'B4');
            xlswrite(xlsPath, [0 diff(sliceLocNew)]', 'FAT', 'C4');
            xlswrite(xlsPath, round(SatVol)', 'FAT', 'D4');
            xlswrite(xlsPath, round(VatVol)', 'FAT', 'E4');
            xlswrite(xlsPath, Loc1', 'FAT', 'F4');
            xlswrite(xlsPath, Loc2', 'FAT', 'G4');
            Pos = 3+numel(VatVol)+1;
            xlswrite(xlsPath, {'Slice #' 'Slice Pos (S)' 'Slice Gap' 'SAT [cm^3]' 'VAT [cm^3]' 'Loc1' 'Loc2'} , 'FAT', ['A' num2str(Pos)]);
            Pos = Pos+2;
            xlswrite(xlsPath, {'Summe'}, 'FAT', ['A' num2str(Pos)]);
            xlswrite(xlsPath, round(sum(SatVol)), 'FAT', ['D' num2str(Pos)]);
            xlswrite(xlsPath, round(sum(VatVol)), 'FAT', ['E' num2str(Pos)]);
            
            
        end
        
        function [ind status message] = LMranging(oComp, range)
            % ranging according to landmark entries
            message = 'OK';
            status = 0;
            ind1 = find(ismember({oComp.pLoc2}, range{1}));
            if strcmp(range{1}, 'none')
                ind2 = numel(oComp);
            elseif isempty(ind1)
                ind1 = 1;       
                message = [range{1} ' not found'];
                status = 1;
            elseif numel(ind1)>1
                ind1 = ind1(1);
                message = [range{1}  'found multiple times'];
                status = 1;
            end
            
            ind2 = find(ismember({oComp.pLoc2}, range{2}));
            if strcmp(range{2}, 'none')
                ind2 = numel(oComp);
            elseif isempty(ind2)
                ind2 = numel(oComp);       
                message = [range{1} ' not found'];
                status = 1;
            elseif numel(ind2)>1
                ind2 = ind2(1);
                message = [range{1}  'found multiple times'];
                status = 1;
            end
            ind = ind1:ind2;
        end
        
        function mCloseReq(oComp, oCont)
            b = questdlg('save data before exit?', 'saveData yes/no', 'yes', 'no', 'yes');
            switch b
                case 'yes'
                    oCont.saveData;
                case 'no'
            end
        end
        
        % % % Object Management % % %
        function dataNew = mTabDat2oCompApp(oComp, dataOld, saveDate)
            %update from ISAT to Dicomflex (only for old ISAT datasets <20180401)
            dataNew = struct(oComp);
            dataNew = repmat(dataNew,1,numel(dataOld));
            % first update data to newest ISAT-tabDat version:
            tabDat = tabDatFatSegmentIOphase;
            tabDat = tabDat.updatetabDatFatSegmentIOphase(dataOld, saveDate);
            tabDat = arrayfun(@struct, tabDat);
            % convert ISAT to Dicomflex
            fields = fieldnames(dataNew);
            for i = 1:numel(fields)
                field = fields{i};
                switch field
                    case {'oBoundaries'}
                        dataNew = setStructArrayField(dataNew, field, {tabDat.Boundaries});
                    case {'pSliceDone'}
                        dataNew = setStructArrayField(dataNew, field, {tabDat.segmentDone});
                        dataNew = setStructArrayField(dataNew, 'pUseSlice', {tabDat.segmentDone});
                    case {'pLoc1'}
                        dataNew = setStructArrayField(dataNew, field, {tabDat.Loc1});
                    case {'pLoc2'}
                        dataNew = setStructArrayField(dataNew, field, {tabDat.Loc2});
                    case {'pFatThresh'}
                        dataNew = setStructArrayField(dataNew, field, {tabDat.fatThresh});
                    case {'oImgs'}
                        oImgTmp = cImageDcm;
                        for j = 1:numel(tabDat)
                            tabDat(j).imgs = oImgTmp.imgDcm2oImageDcm(tabDat(j).imgs);
                        end
                        dataNew = setStructArrayField(dataNew, field, {tabDat.imgs});
                    case {'pDataName'}
                        dataNew = setStructArrayField(dataNew, field, {tabDat.dataName});
                    case {'pHistory'}
                        dataNew = setStructArrayField(dataNew, field, {tabDat.history});
                    case {'pStandardImgType'}
                        dataNew = setStructArrayField(dataNew, field, {tabDat.standardImgType});
                    case {'pPatientName'}
                        dataNew = setStructArrayField(dataNew, field, {tabDat.patientName});
                    case {'pSliceLocation'}
                        dataNew = setStructArrayField(dataNew, field, {tabDat.sliceLocation});
                    case {'pVarious'}
                        dataNew = setStructArrayField(dataNew, field, {tabDat.Various});
                    otherwise
                end
            end
        end
        
        function oComp = mUpdate_cCompute_segfatMR(oComp, data, saveDate)
            % here each slice gets implemented in the current oComp and
            % oComp_segfatMR structure
            for i = 1:numel(data)
                oComp(i) = cCompute_segfatMR;  % object
                oCompTmp = data(i);  % simple variable (struct)
                switch oCompTmp.pVersion_cCompute_segfatMR
                    case {'1.2'}
                        for f = fieldnames(oComp(i))'
                            f = f{1};
                            oComp(i).(f) = oCompTmp.(f);
                        end
                        % take care about cfgversion!!?!!?!!?!!
                        sc = superclasses(class(oComp(i)));
                        oComp(i) = oComp(i).(['mUpdate_' sc{1}]);
                    otherwise
                        msgbox('oCompFatSegment version problem in oCompFatSegment_updateFcn!');
                end
            end
        end
        
        function oComp = cCompute_segfatMR(cCompArray)
            
        end
    end
end


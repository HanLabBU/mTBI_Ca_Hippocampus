function varargout = imsc(im)
im = squeeze(im);
nFrames = size(im,3);
if nFrames > 1
   frameMean = mean2(getDataSample(im));
   frameStd = std2(getDataSample(im));
   frameNum = ceil(rand*nFrames);
   im = squeeze(im(:,:,frameNum));
   txtString = sprintf('Frame: %i\nClass: %s',frameNum,class(im));
else
   frameMean = mean2(im);
   frameStd = mean2(im);
   txtString = sprintf('Class: %s',class(im));
end
warning('off','MATLAB:Figure:SetPosition')

hIm = handle(imagesc(im));
hAx = handle(hIm.Parent);
hFig = handle(hAx.Parent);

hFig.Units = 'normalized';
hFig.Position = [0.01 0.39 0.97 0.56];
hAx.Position = [0 .02 .95 .95];
hAx.PlotBoxAspectRatioMode = 'manual';
hAx.XTick = [];
hAx.YTick = [];
% hAx.CLim = frameMean + [-frameStd frameStd];
hAx.CLim = [frameMean/3 frameMean+5*frameStd];
hFig.Renderer = 'opengl';

n = 256;
redtrans = 50;
greentrans = round(n/20);
bluetrans = 10;
chan.red = [fliplr(linspace(0, .8, redtrans))' ; zeros(n-redtrans,1) ];
chan.green = [zeros(greentrans,1) ; linspace(0, 1, n-greentrans)'];
chan.blue = [fliplr( logspace(1, 2, n-bluetrans)./250)'-log(2)/500 ; linspace(log(2)./500, .5, bluetrans)'];
cmap = [chan.red(:) chan.green(:) chan.blue(:)];
colormap(cmap)
colorbar
text(20,50,txtString);

if nargout
   varargout{1} = hIm;
end
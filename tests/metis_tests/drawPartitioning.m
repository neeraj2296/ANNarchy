function drawPartitioning(csvfile, layersToDraw, layerDistance, azimut, elevation)
% drawPartitioning(csvfile, layersToDraw=allLayers, layerDistance=0,
%                  azimut=105, elevation=15)
%
% draws a three-dimensional partitioning of the given csv-file with the
% following colums layer, x, y, z, part
% the view is set to the angle azimut and elevation

if nargin<5 elevation = 15; end
if nargin<4 azimut = 105; end
if nargin<3 layerDistance = 0; end
if nargin<1 'no csv file specified'; end

% read the csv file
graph = csvread(csvfile);
graphSize = size(graph);
numberNeurons = graphSize(1);

% separate the values
layer = graph(:, 1)';
x = graph(:, 2)';
y = graph(:, 3)';
z = graph(:, 4)';
part = graph(:, 5)';


% define the layers and neurons to be drawn
if nargin<2 layersToDraw=unique(layer); end
neuronNumbers=[1:numberNeurons];
neuronsToDraw = neuronNumbers(ismember(layer, layersToDraw));

layerZSize = zeros(1,length(layersToDraw));
for i = 1:length(layersToDraw),
    tmpZ = z(layer == layersToDraw(i));
    layerZSize(i) = max(tmpZ) - min(tmpZ) + 1 + layerDistance;
end

% the z offset caused by stacking the layers for each neuron
neuronsZOffset = zeros(1,length(neuronsToDraw));
for i = 1:length(neuronsToDraw),
    tmpLayer = layer(neuronsToDraw(i));
    neuronsZOffset(i) = sum(layerZSize(layersToDraw < tmpLayer));
end

% the coordinates of the markers
X = x(neuronsToDraw);
Y = y(neuronsToDraw);
Z = z(neuronsToDraw) + neuronsZOffset;

% size of the markers
S = ones(1, length(neuronsToDraw))*50;

% color of the markers
C = part(neuronsToDraw);

% draw it
scatter3(X(:),Y(:),Z(:),S(:),C(:),'filled'), view(azimut, elevation);

% set the color axis to an absolute value
caxis([min(part) max(part)]);

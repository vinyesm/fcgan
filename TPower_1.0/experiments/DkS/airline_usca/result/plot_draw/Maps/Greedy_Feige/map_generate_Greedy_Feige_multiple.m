clear all;
load('..\..\..\..\data\data_city');
load('..\..\..\Greedy_Feige_Multiple_Result.mat');

dense_subgraph_idx_mul = Greedy_Feige_Multiple_Result.dense_subgraph_idx;
center_city_idx_mul = Greedy_Feige_Multiple_Result.center_city_idx;

axesm('mercator','grid','off','MapLatLimit',[-75 75]); tightmap;
% Map the geostruct with the continent outlines
geoshow('landareas.shp');

track_colors = ['bgrcmykw']; 
for (ds_id = 1:length(dense_subgraph_idx_mul))
    dense_subgraph_idx = dense_subgraph_idx_mul{ds_id};
    center_city_idx = center_city_idx_mul(ds_id);
    % The first field by convention is Geometry (dimensionality).
    % As Geometry is the same for all elements, assign it with deal:
    nCity = length(dense_subgraph_idx);
    Cities = [];
    [Cities(1:nCity).Geometry] = deal('Point');
    
    % Add the latitudes and longitudes to the geostruct:
    

    for (i=1:nCity)
        city_idx = dense_subgraph_idx(i);
        Cities(i).Lat = city_lat_lon(city_idx,1);
        Cities(i).Lon = city_lat_lon(city_idx,2);
        Cities(i).Name = city_name{city_idx};
        if (city_idx == center_city_idx_mul(ds_id))
            center_city_idx = i;
        end
    end
    
    % Map the City locations with filled circular markers
    geoshow(Cities,'Marker','o','MarkerSize', 2,...
        'MarkerFaceColor','c','MarkerEdgeColor','k');
    
    % Call the new geostruct Tracks and give it a line geometry:
    city_track_sub = [];
    for (i=1:nCity)
        city_idx = dense_subgraph_idx(i);
        idx_1 = find(city_track(:,1)==city_idx);
        idx_2 = find(ismember(city_track(idx_1,2) , dense_subgraph_idx));
        city_track_sub = [city_track_sub; city_track(idx_1(idx_2),:)];
    end
    nTracks = size(city_track_sub, 1);
    Tracks = [];
    [Tracks(1:nTracks).Geometry] = deal('Line');
    
    % Create a text field identifying kind of track each entry is.
    % Here they all will be great circles, identified as 'gc'
    % (string signifying great circle arc to certain functions)
    trackType = 'gc';
    [Tracks.Type] = deal(trackType);
    
    % Give each track an identifying name
    for (i=1:nTracks)
        Tracks(i).Name = [city_name{city_track_sub(i,1)}, ' - ', city_name{city_track_sub(i,2)}]; %'Paris-Santiago';
        [Tracks(i).Lat Tracks(i).Lon] = ...
            track2(trackType, city_lat_lon(city_track_sub(i,1),1),city_lat_lon(city_track_sub(i,1),2),city_lat_lon(city_track_sub(i,2),1),city_lat_lon(city_track_sub(i,2),2));
        
    end
    
    
    % The distance function computes distance and azimuth between
    % given points, in degrees. Store both in the geostruct.
    for j = 1:nTracks
        [dist az] = ...
            distance(trackType,Tracks(j).Lat(1),...
            Tracks(j).Lon(1),...
            Tracks(j).Lat(end),...
            Tracks(j).Lon(end));
        [Tracks(j).Length] = dist;
        [Tracks(j).Azimuth] = az;
    end
    
    % On cylindrical projections like Mercator, great circle tracks
    % are curved except those that follow the Equator or a meridian.
    
    % Graphically differentiate the tracks by creating a symbolspec;
    % key line color to track length, using the 'summer' colormap.
    % Symbolspecs make it easy to vary color and linetype by
    % attribute values. You can also specify default symbologies.
    
    colorRange = makesymbolspec('Line',...
        {'Length',[min([Tracks.Length]) ...
        max([Tracks.Length])],...
        'Color',track_colors(ds_id)});
    geoshow(Tracks,'SymbolSpec',colorRange);
    
    % Display the city names using data in the geostruct field Name.
    % Note that you must treat the Name field as a cell array.
    geoshow(Cities(center_city_idx),'Marker','o','MarkerSize', 6,...
        'MarkerFaceColor','r','MarkerEdgeColor','k');
    
    textm([Cities(center_city_idx).Lat],[Cities(center_city_idx).Lon],...
        {Cities(center_city_idx).Name},'FontWeight','bold','FontSize',6);
end

set(gcf, 'PaperPositionMode', 'auto');

print('-depsc2', 'Gready_Feige_multiple.eps');

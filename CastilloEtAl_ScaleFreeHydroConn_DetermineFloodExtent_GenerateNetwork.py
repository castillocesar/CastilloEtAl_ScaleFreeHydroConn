'''
This script is provided as part of the methods section for the manuscript entitled "Scale-free
structure of surface-water connectivity within a lowland river floodplain" by Cesar R Castillo,
İnci Güneralp, Billy Hales, and Burak Güneralp. The corresponding author is Cesar R Castillo
and can be contacted at castillocesar@tamu.edu.  

This script uses tif rasters of inundation depth exported from hec-ras to determine the
surface-water connections between landscape patches in the river-floodplain landscape. The
landscape patches in this particular case are composed of an irregular tessalation that are
represented by a polygon shapefile/feature class with a number of specific fields in the
attribute table. The connectivity analysis is only really focused on portions of the landscape
that are connected to the main channel by surface water and this script also requires a polyline
shapefile/feature class the centerline or thalwag for the main channel to filter out 
inundation that is not connected to the main channel.

It is assumed that all inputs and outputs will be within a source directory that is specified
below (variable is "prj_path").

It is assumed that all the depth rasters are found within a single folder and that they all have
the same spatial reference as the vector data used in the analysis.

All sets of outputs (intermediate and final) will be sent to their own respective folders that
will hold their pertinent outputs.

This code was designed to be used with Python 3.6 and ArcGIS Pro 2.x
'''

# loading libraries/extensions and configuring the analysis environment
import os, arcpy, glob, datetime
start = datetime.datetime.now()
print('script started at', str(datetime.datetime.now()))
arcpy.CheckOutExtension('spatial')
arcpy.env.overwriteOutput = True

###############################################################################################
###############################################################################################
# specifying the inputs to the script

### this is the source directory that will contain all the inputs and outputs
prj_path = 'DIRECTORY' 
# the minimum depth of inundation that will allowed for that inundation pixel to be used in the
# analysis. our model had a calibration deviation with previous estimates of 20 cms and we used
# this value as the minimum amount of inundation to be included in the analysis
min_inun = 'VALUE'
# the subdirectory within the source path that contains the inundation depth rasters
ras_path = os.path.join(os.path.split(prj_path)[0],
                        'SUBDIRECTORY')
# filename of the polygon shapefile/feature class that contains the landscape patches (also
# known as facets). it is assumed that this file is within the source directory.
facets = os.path.join(prj_path,
                      'SHAPEFILE OR FEATURE CLASS')
# filename of the polyline shapefile/feature class that represents the centerline or thalwag for
# the main channel in the landscape. it is assumed that this file is within the source directory.
river_line_shp = os.path.join(prj_path,
                              'SHAPEFILE OR FEATURE CLASS')
# filename of the geoTIFF raster of the digital terrain model (dtm) for the landscape of interest.
# it is assumed that this file is within the source directory.
dtm = os.path.join(prj_path,
                   'TIFF')
# table the depicts that unique identifier that pertains to each landscape classifications used 
# in the analysis
landform_codes = os.path.join(prj_path,
                              'landform_codes.dbf')
# filename of the polygon shapefile/feature class that depicts the study are that will be
# included in the analysis
study_area = os.path.join(prj_path,
                          'StudyArea.shp')
# filename for the polygon shapefile/feature class of the landscape patches/facets that have been 
# clipped to the study area for the anlysis
facets_clipped = os.path.join(prj_path,
                              'EMST_landformclass_withchan_chanagg_soil_studyarea.shp')
# subdirectory depicting where a cleaned version of the geoTIFF raster of the inundation depth is
# contained
depth_ras_path = os.path.join(prj_path,
                              'depth_Ras')

###############################################################################################
###############################################################################################
# specifying the output subdirectories where final and intermediate outputs will be sent
print('started setting up the project environment at', str(datetime.datetime.now()))
# subdirectory where the filtered inundation raster into a polygon shapefile or feature class
outfloodextent_path = os.path.join(prj_path,
                                   'floodextent_poly')
# subdirectory where the polygon shapefiles/feature classes that represent the landscape
# patches/facets that are connected to main channel are sent
facets_connected_path = os.path.join(prj_path,
                                     'facets_connected')
# subdirectory where the point shapefiles/feature classes that represent the connected landscape
# patches/facets are sent
facets_connected_pt_path = os.path.join(prj_path,
                                        'facets_connected_pt')
# subdirectory where the polyline shapefiles/feature classes that represent connections between
# landscape patches/facets are sent
facets_connected_ln_path = os.path.join(prj_path,
                                        'facets_connected_ln')
# subdirectory where the table that depicts which landscape patches/facets are adjacent to each
# other with regard to their spatial location
facets_connected_nbrtbl_path = os.path.join(prj_path,
                                            'facets_connected_nbrtbl')
# sudirectory where the attribute table (in dbf format) of the polyline shapefiles/feature
# classes that represent connections between landscape patches/facets are sent. some of the
# fields in these tables are used to construct an edge-list that is used in the network analysis
facets_connected_ln_tbl_path = os.path.join(prj_path,
                                            'facets_connected_ln_tbl')
# subdirectory where the binary/boolean geoTIFF raster of the inundation extent is sent
bool_ras_path = os.path.join(prj_path,
                             'bool_floodextent')
# subdirectory where the portions of the bindary/boolean geoTIFF raster of the inundation extent
# that are connected to the main channel are sent
bool_ras_conn_path = os.path.join(prj_path,
                                 'bool_floodextent_conn')
# subdirectory where the polygon shapefiles/feature classes of the landscape patches/feature
# classes that intersect the inundation extent polygon shapefile/feature class
facets_inun_intersect_path = os.path.join(prj_path,
                                          'facets_inun_intersect')
# subdirectory where the summary the dbf file of the summary statistics performed on the 
# intersected polygon shapefile/feature class between landscape patches/facets and inundation
# extent
facets_inun_intersect_summary_path = os.path.join(prj_path,
                                                  'facets_inun_intersect_summary')
# subdirectory where the cleaned landscape patches/facets after filtering are sent
facets_clean_path = os.path.join(prj_path,
                                 'landforms_clean.shp')
# subdirectory where geoTIFFs rasters of flow accumulation are sent
flow_acc_ras = os.path.join(prj_path,
                            'flow_acc.tif')

#####################################################################################################
# specifying some of the "environment" parameters for the geoprocessing
arcpy.env.outputCoordinateSystem = arcpy.Describe(facets).spatialReference
arcpy.env.extent = arcpy.Describe(facets).extent
arcpy.env.snapRaster = dtm
arcpy.env.cellSize = arcpy.Describe(dtm).meanCellHeight

#####################################################################################################
# formating some of the inputs in order to make sure they have the appropriate fields and fieldnames
print('add necessary fields to the facets feature class at', str(datetime.datetime.now()))
# clipping the landforms to the study area
arcpy.analysis.Clip(in_features=facets,
                    clip_features=study_area,
                    out_feature_class=facets_clipped)
facets = facets_clipped
del facets_clipped

# assigning a lf_code to each landform from ecognition
if len(arcpy.ListFields(facets, 'lf_code')) == 0:
    landformcodes_list = []
    with arcpy.da.SearchCursor(in_table=landform_codes, field_names=['landform', 'lf_code']) as cursor:
        for row in cursor:
            landformcodes_list.append(row)

    arcpy.management.AddField(in_table=facets,
                              field_name='lf_code',
                              field_type='LONG')
    with arcpy.da.UpdateCursor(in_table=facets, field_names=['landform', 'lf_code']) as cursor:
        for row in cursor:
            for landformcodes in landformcodes_list:
                if row[0] == landformcodes[0]:
                    row[1] = landformcodes[1]
            cursor.updateRow(row)

# end of if statement block
# adding an ID field to the facet attribute table
if len(arcpy.ListFields(facets, 'ID')) == 0:
    arcpy.management.AddField(in_table=facets,
                              field_name='ID',
                              field_type='LONG')
    arcpy.management.CalculateField(in_table=facets,
                                    field='ID',
                                    expression='!FID! + 1',
                                    expression_type='PYTHON_9.3')

# end of if statement block

if len(arcpy.ListFields(facets, 'conn_h_m')) == 0:
    arcpy.management.AddField(in_table=facets,
                              field_name='conn_h_m',
                              field_type='FLOAT')
else:
    arcpy.management.CalculateField(in_table=facets,
                                    field='conn_h_m',
                                    expression=0.0,
                                    expression_type='PYTHON_9.3')
# end of if else statement block

# calculating some geometry attributes for the facets
arcpy.management.AddGeometryAttributes(Input_Features=facets,
                                       Geometry_Properties=['AREA',
                                                            'CENTROID_INSIDE',
                                                            'PERIMETER_LENGTH'],
                                       Length_Unit='METERS',
                                       Area_Unit='SQUARE_METERS')

####################################################################################

# this block of code converts the depth raster outputed from hec-ras to a polygon
# shapefile that will be used in the connectivity analysis

print('started the flood inundation conversion to polygons at', str(datetime.datetime.now()))
wrr_list = [folders for folders in os.listdir(ras_path) if folders.startswith('WRR')]
for i in range(0, len(wrr_list)):
    ras = [file for file in os.listdir(os.path.join(ras_path, wrr_list[i])) if (file.startswith('Depth') and file.endswith('.tif'))]
    shp_name = 'floodextent_' + wrr_list[i].split('_')[1] + '.shp'
    ras_name = 'floodextent_' + wrr_list[i].split('_')[1] + '.tif'
    depth_ras = os.path.join(os.path.join(ras_path,
                                          wrr_list[i]),
                             ras[0])
    shp_path = os.path.join(outfloodextent_path,
                            shp_name)
    ras_bool_name = os.path.join(bool_ras_path,
                                 ras_name)
    ras_bool_conn_name = os.path.join(bool_ras_conn_path,
                                      shp_name)
    ras_arcpy = arcpy.Raster(depth_ras)
    ras_bool = arcpy.sa.Con(in_conditional_raster=ras_arcpy > min_inun,
                            in_true_raster_or_constant=1)
    ras_bool.save(ras_bool_name)
    
    # converting the boolean inundation raster into a polygon shapefile
    floodextent_shp = arcpy.conversion.RasterToPolygon(in_raster=ras_bool,
                                                       out_polygon_features=shp_path,
                                                       simplify='NO_SIMPLIFY',
                                                       raster_field='VALUE')

    floodextent_lyr = arcpy.management.MakeFeatureLayer(in_features=floodextent_shp,
                                                        out_layer='floodextent_lyr')
    arcpy.management.SelectLayerByLocation(in_layer=floodextent_lyr,
                                           overlap_type='INTERSECT',
                                           select_features=river_line_shp,
                                           selection_type='NEW_SELECTION')

    floodextent_chan_shp = arcpy.management.CopyFeatures(in_features=floodextent_lyr,
                                                         out_feature_class=ras_bool_conn_name)

    # delete unnecessary fields from the polygon shapefile
    arcpy.management.DeleteField(in_table=floodextent_chan_shp,
                                 drop_field=['Join_Count',
                                             'TARGET_FID',
                                             'gridcode',
                                             'Id_1'])

    print('finished the flood inundation conversion for', ras_name, 'at', str(datetime.datetime.now()))
    del ras, shp_name, ras_name, depth_ras, shp_path, floodextent_shp, floodextent_chan_shp, floodextent_lyr
        

# end of for loop

######################################################################################################

'''
defining the graphscript function
this function draws lines that connects facet nodes (centroids in this case) that are identified as
being connected (i.e. neighbors) in order to construct graphs/networks for connected and not
connected portions

inputs to this particular function include:
shp is a point shapefile that identifies the centroid of a respective facet
dbf is a dbf table that contains records with a source facet identifier code and the identifier for
one of its neighbors
out is file path with file type (.shp) where the output line shapefile will be saved
'''
def graphscript(pt_in, nbr_tbl, ln_out):
        # constructing a search cursor that returns the facet code and the associated geometry
	sCurs1 = arcpy.da.SearchCursor(in_table=pt_in,
                                       field_names=['PID',
                                                    'shape@'])
	# building an empty dictionary that will hold the geometric coordinates for facets centroids
	xydict = {}
	# iterating over the elements in the search cursor and populating the previously constructed dictionary
	for item in sCurs1:
            xydict[item[0]] = (item[1].getPart(0).X,
                               item[1].getPart(0).Y)

	del sCurs1
	# end of loop for cursor
	
	# building an empty polyline shapefile that wil make-up the graph for the particular stage
	# this shapefile will also contain two that identify the source and neighbour facets
	arcpy.management.CreateFeatureclass(out_path=os.path.split(ln_out)[0],
                                            out_name=os.path.split(ln_out)[1],
                                            geometry_type='POLYLINE')
	newfields_list = [['src_PID', 'LONG'],
                          ['nbr_PID', 'LONG']]
	for newfields in newfields_list:
            arcpy.management.AddField(in_table=ln_out,
                                      field_name=newfields[0],
                                      field_type=newfields[1])

        del newfields, newfields_list
        # end of loop for cursor

	# constructing an insert cursor that inserts facet centroid attributes to out shapefile
	iCurs = arcpy.da.InsertCursor(in_table=ln_out,
                                      field_names=['shape@',
                                                   'src_PID',
                                                   'nbr_PID'])

	# constructing a search cursor that compiles the source and neighbour facets codes from dbf table
	sCurs2 = arcpy.da.SearchCursor(in_table=nbr_tbl,
                                       field_names=['src_PID',
                                                    'nbr_PID',
                                                    'src_c_hm',
                                                    'nbr_c_hm',
                                                    'src_dmean',
                                                    'nbr_dmean',
                                                    'src_drng',
                                                    'nbr_drng',
                                                    'src_fa',
                                                    'nbr_fa'])
	
	# constructing an empty list and populating it with data from dbf table search cursor
	# this operation also eliminates any duplicates in the list
	nbrList = []
	rejectList = []
        for item in sCurs2:
            if rejectList.count((item[0], item[1])) == 0:
                if item[2] < item[3]:
                    rejectList.append((item[1],
                                       item[0]))
                    nbrList.append((item[0],
                                    item[1]))
                elif item[2] == item[3]:
                    if item[4] > item[5]:
                        rejectList.append((item[1],
                                           item[0]))
                        nbrList.append((item[0],
                                        item[1]))
                    elif item[6] > item[7]:
                        rejectList.append((item[1],
                                           item[0]))
                        nbrList.append((item[0],
                                        item[1]))
                    elif item[8] < item[9]:
                        rejectList.append((item[1],
                                           item[0]))
                        nbrList.append((item[0],
                                        item[1]))
        
        # end of for loop

	pnt = arcpy.Point()
	ary = arcpy.Array()
	for item in nbrList:
            pnt.X = xydict[int(item[0])][0]
            pnt.Y = xydict[int(item[0])][1]
            ary.add(pnt)
            pnt.X = xydict[int(item[1])][0]
            pnt.Y = xydict[int(item[1])][1]
            ary.add(pnt)
	    ply = arcpy.Polyline(ary)
	    iCurs.insertRow([ply,
                             item[0],
                             item[1]])
	    ary.removeAll()
	del sCurs2, iCurs

# end of function
print('finished the creation of graphscript at', str(datetime.datetime.now()))

####################################################################################################################
# the code in this box loops through the set of geospatial files and performs the overlay and 
# connectivity analysis

# creating list of shp files that will be processed using the for loop below
print('started the connection and neighbor analysis at', str(datetime.datetime.now()))
shps = glob.glob(os.path.join(bool_ras_conn_path,
                              '*.shp'))
for shp in shps:
    # selecting the facets and floodextent polygons that will be considering connected to the river
    selfloodextentname = os.path.join(outfloodextent_path,
                                      os.path.splitext(os.path.split(shp)[1])[0])
    inun_lyr = arcpy.management.MakeFeatureLayer(in_features=shp,
                                                 out_layer='floodextentlayer')
    arcpy.management.SelectLayerByLocation(in_layer=inun_lyr,
                                           overlap_type='INTERSECT',
                                           select_features=river_line_shp,
                                           selection_type='NEW_SELECTION')
    inun_shp = arcpy.management.CopyFeatures(in_features=inun_lyr,
                                             out_feature_class=os.path.join(outfloodextent_path,
                                                                            os.path.split(shp)[1]))
    arcpy.management.AddGeometryAttributes(Input_Features=inun_shp,
                                           Geometry_Properties=['AREA',
                                                                'CENTROID_INSIDE',
                                                                'PERIMETER_LENGTH'],
                                           Length_Unit='METERS',
                                           Area_Unit='SQUARE_METERS')
    facets_lyr = arcpy.management.MakeFeatureLayer(in_features=facets,
                                                   out_layer='facetslayer')
                                                   
    arcpy.management.SelectLayerByLocation(in_layer=facets_lyr,
                                           overlap_type='INTERSECT',
                                           select_features=inun_shp,
                                           selection_type='SUBSET_SELECTION')
    # use a temporary facets layer to do the filtering minimum inundated area
    facets_temp = os.path.join(prj_path,
                               'facets_temp.shp')
    arcpy.management.CopyFeatures(in_features=facets_lyr,
                                  out_feature_class=facets_temp)
    # intersecting the connected facets with the river inundation shapefile to determine what fraction of the
    # facet area is inundated for the inundation scenario
    facets_inun_intersect_name = os.path.splitext(os.path.split(facets)[1].split('_')[1])[0] + '_inun_intersect_' + os.path.splitext(os.path.split(shp)[1].split('_')[1])[0] + '.shp'
    facets_inun_intersect = os.path.join(facets_inun_intersect_path,
                                         facets_inun_intersect_name)
    arcpy.analysis.Intersect(in_features=[facets_temp,
                                          inun_shp],
                             out_feature_class=facets_inun_intersect,
                             join_attributes='ALL',
                             output_type='INPUT')
    facets_inun_intersect_singlepart = os.path.join(prj_path,
                                                    'facets_inun_intersect.shp')
    arcpy.management.MultipartToSinglepart(in_features=facets_inun_intersect,
                                           out_feature_class=facets_inun_intersect_singlepart)
    facets_inun_intersect = facets_inun_intersect_singlepart
    del facets_inun_intersect_singlepart
    # add flow accumulation to the facets attribute table
    flow_acc_tbl = os.path.join(prj_path,
                                'flow_acc_tbl.dbf')
    arcpy.sa.ZonalStatisticsAsTable(in_zone_data=facets_inun_intersect,
                                    zone_field='FID',
                                    in_value_raster=flow_acc_ras,
                                    out_table=flow_acc_tbl,
                                    ignore_nodata='DATA',
                                    statistics_type='ALL')

    if len(arcpy.ListFields(facets_inun_intersect, 'flow_acc')) == 0:
        flowacc_list = []
        with arcpy.da.SearchCursor(in_table=flow_acc_tbl, field_names=['OID', 'SUM']) as cursor:
            for row in cursor:
                flowacc_list.append(row)
                
        arcpy.management.AddField(in_table=facets_inun_intersect,
                                  field_name='flow_acc',
                                  field_type='DOUBLE')
        with arcpy.da.UpdateCursor(in_table=facets_inun_intersect, field_names=['FID', 'flow_acc']) as cursor:
            for row in cursor:
                for flowacc in flowacc_list:
                    if row[0] == flowacc[0]:
                        row[1] = flowacc[1]
                cursor.updateRow(row)
    # end of code block

    # removing unnecessary fields
    fields = ['FID_facets',
              'POLY_AREA',
              'INSIDE_X',
              'INSIDE_Y',
              'PERIMETER',
              'FID_floode',
              'ID_1',
              'POLY_ARE_1',
              'INSIDE_X_1',
              'INSIDE_Y_1',
              'PERIMETE_1']
    for i in range(0, len(fields)):
        arcpy.management.DeleteField(in_table=facets_inun_intersect,
                                     drop_field=fields[i])
    # end of for loop

    # re-calculate the feature geometry information after the intersect
    arcpy.management.AddGeometryAttributes(Input_Features=facets_inun_intersect,
                                           Geometry_Properties=['AREA',
                                                                'CENTROID_INSIDE',
                                                                'PERIMETER_LENGTH'],
                                           Length_Unit='METERS',
                                           Area_Unit='SQUARE_METERS')
    # removing the patches that are less than 100m2 in size because the
    # EMST pixels is 10mX10m in size
    inun_facets_lyr = arcpy.management.MakeFeatureLayer(in_features=facets_inun_intersect,
                                                        out_layer='inunfacetslyr')
    arcpy.management.SelectLayerByAttribute(in_layer_or_view=inun_facets_lyr,
                                            selection_type='NEW_SELECTION',
                                            where_clause='"POLY_AREA" > 100')
    
    # final facets layer
    facets_connected_name = os.path.splitext(os.path.split(facets)[1].split('_')[1])[0] + '_' + os.path.split(shp)[1].split('_')[1]
    facets_connected = os.path.join(facets_connected_path,
                                    facets_connected_name)
    # use inundated facets layer to create the final facets layer
    arcpy.management.CopyFeatures(in_features=inun_facets_lyr,
                                  out_feature_class=facets_connected)
    # creating a patch ID for the connected facets
    arcpy.management.AddField(in_table=facets_connected,
                                  field_name='PID',
                                  field_type='LONG')
    arcpy.management.CalculateField(in_table=facets_connected,
                                    field='PID',
                                    expression='!FID! + 1',
                                    expression_type='PYTHON')
    
    # need to instill some memory into the facets attribute table
    # should add a new field to the attribute table that identifies
    # the flow in which this particular facet first becomes inundated
    stage = float(os.path.split(facets_connected)[1].split('_')[1][:3])/10
    with arcpy.da.UpdateCursor(in_table=facets_connected, field_names=['conn_h_m']) as cursor:
        for row in cursor:
            if row[0] == 0:
                row[0] = stage
            cursor.updateRow(row)

    # add inundation depth to the facets attribute table
    depth_ras_list = [files for files in os.listdir(depth_ras_path) if files.endswith('.tif')]
    for depth_ras in depth_ras_list:
        depth_ras_stage = float(depth_ras.split('_')[1].split('cm')[0])/100
        if depth_ras_stage == stage:
            depth_zstats_tbl = os.path.join(prj_path,
                                            'depth_zstats_tbl.dbf')
            arcpy.sa.ZonalStatisticsAsTable(in_zone_data=facets_connected,
                                            zone_field='PID',
                                            in_value_raster=os.path.join(depth_ras_path,
                                                                         depth_ras),
                                            out_table=depth_zstats_tbl,
                                            statistics_type='ALL')
    # end of for loop

    # adding the depth stats 
    keeper_fields = ['PID',
                     'MEAN',
                     'RANGE']
    attributes_list2 = []
    with arcpy.da.SearchCursor(in_table=depth_zstats_tbl, field_names=keeper_fields) as cursor:
        for row in cursor:
            attributes_list2.append(row)
    # end of search cursor
    newfields_list = [['depth_mean', 'DOUBLE'],
                      ['depth_rng', 'DOUBLE']]
    for newfields in newfields_list:
        arcpy.management.AddField(in_table=facets_connected,
                                  field_name=newfields[0],
                                  field_type=newfields[1])
    # end of for loop adding new fields to the nbr table
    fieldnames = ['PID',
                  'depth_mean',
                  'depth_rng']
    with arcpy.da.UpdateCursor(in_table=facets_connected, field_names=fieldnames) as cursor:
        for row in cursor:
            for attributes in attributes_list2:
                if row[0] == attributes[0]:
                    row[1] = attributes[1]
                    row[2] = attributes[2]
            cursor.updateRow(row)
    # end of code block
    # constructing the nbr table that is used by graphscript and other code to determine the connections
    facets_connected_nbrtbl_name = os.path.splitext(os.path.split(facets)[1].split('_')[1])[0] + '_nbrtbl_' + os.path.splitext(os.path.split(shp)[1])[0].split('_')[1] + '.dbf'
    facets_connected_nbrtbl = arcpy.management.CreateTable(out_path=facets_connected_nbrtbl_path,
                                                           out_name=facets_connected_nbrtbl_name)
    arcpy.analysis.PolygonNeighbors(in_features=facets_connected,
                                    out_table=facets_connected_nbrtbl,
                                    in_fields='PID',
                                    area_overlap='NO_AREA_OVERLAP',
                                    both_sides='BOTH_SIDES',
                                    out_linear_units='METERS',
                                    out_area_units='SQUARE_METERS')

    # addin and calculating the s_n_code attribute to facets_connected_nbrtbl that
    # is used in graphscript to identify where to begin and end lines
    arcpy.management.AddField(in_table=facets_connected_nbrtbl,
                              field_name='s_n_code',
                              field_type='TEXT',
                              field_is_nullable='NULLABLE',
                              field_is_required='NON_REQUIRED')
    arcpy.management.CalculateField(in_table=facets_connected_nbrtbl,
                                    field='s_n_code',
                                    expression='\'s\' + str(int( !src_PID!)) + \'n\' + str(int( !nbr_PID!))',
                                    expression_type='PYTHON_9.3')

    # adding more attributes to the nbr table that will be used to determine the connections using graph script
    keeper_fields = ['PID',
                     'ID',
                     'flow_acc',
                     'conn_h_m',
                     'POLY_AREA',
                     'PERIMETER',
                     'landform',
                     'lf_code',
                     'depth_mean',
                     'depth_rng']
    attributes_list2 = []
    with arcpy.da.SearchCursor(in_table=facets_connected, field_names=keeper_fields) as cursor:
        for row in cursor:
            attributes_list2.append(row)
    # end of search cursor
    newfields_list = [['src_ID', 'LONG'],
                      ['nbr_ID', 'LONG'],
                      ['src_c_hm', 'FLOAT'],
                      ['nbr_c_hm', 'FLOAT'],
                      ['src_fa', 'DOUBLE'],
                      ['nbr_fa', 'DOUBLE'],
                      ['src_area', 'DOUBLE'],
                      ['nbr_area', 'DOUBLE'],
                      ['src_per', 'DOUBLE'],
                      ['nbr_per', 'DOUBLE'],
                      ['src_lf', 'TEXT'],
                      ['nbr_lf', 'TEXT'],
                      ['src_lfc', 'LONG'],
                      ['nbr_lfc', 'LONG'],
                      ['src_dmean', 'DOUBLE'],
                      ['nbr_dmean', 'DOUBLE'],
                      ['src_drng', 'DOUBLE'],
                      ['nbr_drng', 'DOUBLE']]
    for newfields in newfields_list:
        arcpy.management.AddField(in_table=facets_connected_nbrtbl,
                                  field_name=newfields[0],
                                  field_type=newfields[1])
    # end of for loop adding new fields to the nbr table
    fieldnames = ['src_PID',
                  'nbr_PID',
                  'src_ID',
                  'nbr_ID',
                  'src_fa',
                  'nbr_fa',
                  'src_c_hm',
                  'nbr_c_hm',
                  'src_area',
                  'nbr_area',
                  'src_per',
                  'nbr_per',
                  'src_lf',
                  'nbr_lf',
                  'src_lfc',
                  'nbr_lfc',
                  'src_dmean',
                  'nbr_dmean',
                  'src_drng',
                  'nbr_drng']
    with arcpy.da.UpdateCursor(in_table=facets_connected_nbrtbl, field_names=fieldnames) as cursor:
        for row in cursor:
            for attributes in attributes_list2:
                if row[0] == attributes[0]:
                    row[2] = attributes[1]
                    row[4] = attributes[2]
                    row[6] = attributes[3]
                    row[8] = attributes[4]
                    row[10] = attributes[5]
                    row[12] = attributes[6]
                    row[14] = attributes[7]
                    row[16] = attributes[8]
                    row[18] = attributes[9]
                if row[1] == attributes[0]:
                    row[3] = attributes[1]
                    row[5] = attributes[2]
                    row[7] = attributes[3]
                    row[9] = attributes[4]
                    row[11] = attributes[5]
                    row[13] = attributes[6]
                    row[15] = attributes[7]
                    row[17] = attributes[8]
                    row[19] = attributes[9]
            cursor.updateRow(row)
    # end of code block

    facets_connected_pt_name = os.path.splitext(os.path.split(facets)[1].split('_')[1])[0] + '_pt_' + os.path.split(shp)[1].split('_')[1]
    facets_connected_pt = os.path.join(facets_connected_pt_path,
                                       facets_connected_pt_name)
    arcpy.management.CreateFeatureclass(out_path=os.path.split(facets_connected_pt)[0],
                                        out_name=os.path.split(facets_connected_pt)[1],
                                        geometry_type='POINT')
    # adding the necesarry attributes to the point feature class
    fields = arcpy.ListFields(facets_connected)
    fields_list = []
    for i in range(0, len(fields)):
        fields_list.append([fields[i].name, fields[i].type])
    # end of the for loop
    fields_list = fields_list[2:]
    for item in fields_list:
        arcpy.management.AddField(in_table=facets_connected_pt,
                                  field_name=item[0],
                                  field_type=item[1])
    # end of for loop
    
    # gathering list of fields to use in cursors
    fields = arcpy.ListFields(facets_connected)
    fields_list = []
    for i in range(0, len(fields)):
        fields_list.append(fields[i].name)
    # end of the for loop
    # removing uneccessary fields and adding the geometry access keyword
    fields_list = ['shape@'] + fields_list[2:]
    scurs = arcpy.da.SearchCursor(in_table=facets_connected,
                                  field_names=fields_list)
    icurs = arcpy.da.InsertCursor(in_table=facets_connected_pt,
                                  field_names=fields_list)
    for item in scurs:
        icurs.insertRow((item[0].labelPoint,
                         item[1],
                         item[2],
                         item[3],
                         item[4],
                         item[5],
                         item[6],
                         item[7],
                         item[8],
                         item[9],
                         item[10],
                         item[11],
                         item[12],
                         item[13])) 
    del scurs, icurs
    
    # construct the line shapefile for the facets classified as being connected to the main channel
    facets_connected_ln_name = os.path.splitext(os.path.split(facets)[1].split('_')[1])[0] + '_ln_' + os.path.split(shp)[1].split('_')[1]
    facets_connected_ln = os.path.join(facets_connected_ln_path,
                                       facets_connected_ln_name)
    graphscript(pt_in=facets_connected_pt,
                nbr_tbl=facets_connected_nbrtbl,
                ln_out=facets_connected_ln)

    # adding and calculating the s_n_code field that is used to understand the connections
    arcpy.management.AddField(in_table=facets_connected_ln,
                              field_name='s_n_code',
                              field_type='TEXT',
                              field_is_nullable='NULLABLE',
                              field_is_required='NON_REQUIRED')
    arcpy.management.CalculateField(in_table=facets_connected_ln,
                                    field='s_n_code',
                                    expression='\'s\' + str(int( !src_PID!)) + \'n\' + str(int( !nbr_PID!))',
                                    expression_type='PYTHON_9.3')

    # adding attributes to facets_connected_ln that will allow 
    newfields_list = [['src_ID', 'LONG'],
                      ['nbr_ID', 'LONG'],
                      ['src_fa', 'DOUBLE'],
                      ['nbr_fa', 'DOUBLE'],
                      ['src_c_hm', 'FLOAT'],
                      ['nbr_c_hm', 'FLOAT'],
                      ['src_area', 'DOUBLE'],
                      ['nbr_area', 'DOUBLE'],
                      ['src_per', 'DOUBLE'],
                      ['nbr_per', 'DOUBLE'],
                      ['src_lf', 'TEXT'],
                      ['nbr_lf', 'TEXT'],
                      ['src_lfc', 'LONG'],
                      ['nbr_lfc', 'LONG'],
                      ['src_dmean', 'DOUBLE'],
                      ['nbr_dmean', 'DOUBLE'],
                      ['src_drng', 'DOUBLE'],
                      ['nbr_drng', 'DOUBLE']]
    for newfields in newfields_list:
        arcpy.management.AddField(in_table=facets_connected_ln,
                                  field_name=newfields[0],
                                  field_type=newfields[1])
    fieldnames = ['src_PID',
                  'nbr_PID',
                  'src_ID',
                  'nbr_ID',
                  'src_fa',
                  'nbr_fa',
                  'src_c_hm',
                  'nbr_c_hm',
                  'src_area',
                  'nbr_area',
                  'src_per',
                  'nbr_per',
                  'src_lf',
                  'nbr_lf',
                  'src_lfc',
                  'nbr_lfc',
                  'src_dmean',
                  'nbr_dmean',
                  'src_drng',
                  'nbr_drng']
    with arcpy.da.UpdateCursor(in_table=facets_connected_ln, field_names=fieldnames) as cursor:
        for row in cursor:
            for attributes in attributes_list2:
                if row[0] == attributes[0]:
                    row[2] = attributes[1]
                    row[4] = attributes[2]
                    row[6] = attributes[3]
                    row[8] = attributes[4]
                    row[10] = attributes[5]
                    row[12] = attributes[6]
                    row[14] = attributes[7]
                    row[16] = attributes[8]
                    row[18] = attributes[9]
                if row[1] == attributes[0]:
                    row[3] = attributes[1]
                    row[5] = attributes[2]
                    row[7] = attributes[3]
                    row[9] = attributes[4]
                    row[11] = attributes[5]
                    row[13] = attributes[6]
                    row[15] = attributes[7]
                    row[17] = attributes[8]
                    row[19] = attributes[9]
            cursor.updateRow(row)
    # end of code block

    # calculating the length of each of the lines constructed in graph script that can be used
    # to infer things about the network and its topology
    arcpy.management.AddGeometryAttributes(Input_Features=facets_connected_ln,
                                           Geometry_Properties='LENGTH',
                                           Length_Unit='METERS')

    attributes_list3 = []
    with arcpy.da.SearchCursor(in_table=facets_connected, field_names=['ID', 'conn_h_m']) as cursor:
        for row in cursor:
            attributes_list3.append(row)
    # end of code block

    # using an update cursor to change the values of conn_h_m and connmj_h_m in the facets shapefile
    stage = float(os.path.split(facets_connected)[1].split('_')[1][:3])/10
    with arcpy.da.UpdateCursor(in_table=facets, field_names=['ID', 'conn_h_m']) as cursor:
        for row in cursor:
            if row[1] == 0:
                for attributes in attributes_list3:
                    if row[0] == attributes[0]:
                        row[1] = stage
            cursor.updateRow(row)
    # end of code block
    
    # exporting the attribute table of the facets_connected_ln shapefile so that
    # it can be used to process the connections in R
    fieldnames = ['src_PID',
                  'nbr_PID',
                  'src_ID',
                  'nbr_ID',
                  'src_fa',
                  'nbr_fa',
                  'src_c_hm',
                  'nbr_c_hm',
                  'src_area',
                  'nbr_area',
                  'src_per',
                  'nbr_per',
                  'src_lf',
                  'nbr_lf',
                  'src_lfc',
                  'nbr_lfc',
                  'src_dmean',
                  'nbr_dmean',
                  'src_drng',
                  'nbr_drng']
    facets_connected_ln_tbl_name = os.path.splitext(os.path.split(facets)[1].split('_')[1])[0] + '_ln_tbl_' + os.path.splitext(os.path.split(shp)[1].split('_')[1])[0] + '.txt'
    arcpy.stats.ExportXYv(Input_Feature_Class=facets_connected_ln,
                          Value_Field=fieldnames,
                          Delimiter='COMMA',
                          Output_ASCII_File=os.path.join(facets_connected_ln_tbl_path,
                                                         facets_connected_ln_tbl_name),
                          Add_Field_Names_to_Output='ADD_FIELD_NAMES')
    
    print('finished the connectivity and neighbor analysis for', os.path.split(shp)[1], 'at', str(datetime.datetime.now()))

# end of this crazy for loop

##################################################################################################################

# end of the script and printing some messages depicting how long the processing time was
end = datetime.datetime.now()
ex_time = end - start
print('finished the entire script at', str(datetime.datetime.now()))
print('execution time for script:', str(ex_time))



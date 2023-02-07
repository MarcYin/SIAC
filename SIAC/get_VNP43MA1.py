import os
import json
import getpass
import requests
import numpy as np
from retry import retry
from functools import partial
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timedelta
from osgeo import gdal, ogr, osr

home = os.path.expanduser("~")
file_path = os.path.dirname(os.path.realpath(__file__))
test_url = 'https://e4ftl01.cr.usgs.gov/VIIRS/VNP43MA1.001/2022.11.08/VNP43MA1.A2022312.h19v08.001.2022321114550.h5'

def get_auth():
    try:
        auth = tuple([os.environ['Earthdata_user'], os.environ['Earthdata_pass']])
        with requests.Session() as s:                                            
            s.auth = auth                                                        
            r1     = s.get(test_url)                                             
            r      = s.get(r1.url, stream=True, headers={'user-agent': 'My app'}) 
        if r.status_code == 401:
            print('Wrong username and password are set for Earthdata_user and Earthdata_pass Environment variables.') 
        else:
            with open(file_path + '/data/.earthdata_auth', 'wb') as f:
                for i in auth:               
                    f.write((i + '\n').encode())
    except:
        #print('Environment variables Earthdata_user and Earthdata_pass for accessing NASA Earthdata are not set, and trying to read from file or input.')
        pass
    if os.path.exists(file_path + '/data/.earthdata_auth'):
        try:
            username, password = np.loadtxt(file_path + '/data/.earthdata_auth', dtype=str)
            auth = tuple([username, password])
            with requests.Session() as s:
                s.auth = auth
                r1     = s.get(test_url)
                r      = s.get(r1.url, stream=True, headers={'user-agent': 'My app'}) 
            while r.status_code == 401:                 
                print('Wrong username and password stored, please enter again')
                username = input('Username for NASA Earthdata: ')      
                password = getpass.getpass('Password for NASA Earthdata: ')
                auth = tuple([username, password])
                with requests.Session() as s:
                    s.auth = auth
                    r1     = s.get(test_url)
                    r      = s.get(r1.url, stream=True, headers={'user-agent': 'My app'}) 
                
            os.remove(file_path + '/data/.earthdata_auth')    
            with open(file_path + '/data/.earthdata_auth', 'wb') as f:
                for i in auth:                                 
                    f.write((i + '\n').encode())              
            auth = tuple([username, password])

        except:
            print('Please provide NASA Earthdata username and password for downloading MCD43 data, which can be applied here: https://urs.earthdata.nasa.gov.')
            username = input('Username for NASA Earthdata: ')
            password = getpass.getpass('Password for NASA Earthdata: ')
            auth = username, password                              
            with requests.Session() as s:
                s.auth = auth
                r1     = s.get(test_url)
                r      = s.get(r1.url, stream=True, headers={'user-agent': 'My app'}) 
            while r.status_code == 401:                 
                print('Wrong username and password typed, please enter again')
                username = input('Username for NASA Earthdata: ')      
                password = getpass.getpass('Password for NASA Earthdata: ')
                auth = tuple([username, password])
                with requests.Session() as s:
                    s.auth = auth
                    r1     = s.get(test_url)
                    r      = s.get(r1.url, stream=True, headers={'user-agent': 'My app'}) 
            os.remove(file_path + '/data/.earthdata_auth')
            with open(file_path + '/data/.earthdata_auth', 'wb') as f:     
                for i in auth:                                     
                    f.write((i + '\n').encode())
            auth = tuple([username, password])
    else:
        username = input('Username for NASA Earthdata: ')
        password = getpass.getpass('Password for NASA Earthdata: ')
        auth = username, password
        with requests.Session() as s:
            s.auth = auth
            r1     = s.get(test_url)
            r      = s.get(r1.url, stream=True, headers={'user-agent': 'My app'}) 
        while r.status_code == 401:                 
            print('Wrong username and password typed, please enter again')
            username = input('Username for NASA Earthdata: ')      
            password = getpass.getpass('Password for NASA Earthdata: ')
            auth = tuple([username, password])
            with requests.Session() as s:
                s.auth = auth
                r1     = s.get(test_url)
                r      = s.get(r1.url, stream=True, headers={'user-agent': 'My app'}) 
        password_file = file_path + '/data/.earthdata_auth'
        os.makedirs(os.path.dirname(password_file), exist_ok=True)
        with open(password_file, 'wb') as f:
            for i in auth: 
                f.write((i + '\n').encode())
        auth = tuple([username, password])
    return auth

def get_boundary(raster, to_wgs84 = True):
    '''
    Take a raster and find boundary of the raster 
    if to_wgs84 = True, transform them into wgs84 coordinates
    with at least meters resoluton.
    Return geojson(s)
    '''
    g    = gdal.Open(raster)
    geom = g.GetGeoTransform()
    x_coords = geom[0] + np.arange(g.RasterXSize + 1) * geom[1]
    y_coords = geom[3] + np.arange(g.RasterYSize + 1) * geom[5]
    north_boundary = np.dstack([x_coords,       np.ones(g.RasterXSize + 1) * y_coords[0]])
    south_boundary = np.dstack([x_coords[::-1], np.ones(g.RasterXSize + 1) * y_coords[-1]])
    east_boundary  = np.dstack([x_coords[-1] * np.ones(g.RasterYSize + 1), y_coords])
    west_boundary  = np.dstack([x_coords[0]  * np.ones(g.RasterYSize + 1), y_coords[::-1]])
    boundary = np.hstack([north_boundary, east_boundary, south_boundary, west_boundary]).tolist()   
    polygon  = {'type': 'Polygon', 'coordinates': boundary}
    geometry = ogr.CreateGeometryFromJson(json.dumps(polygon))
    sgeometry = geometry.SimplifyPreserveTopology(abs(geom[1])/1000)
    spolygon = json.loads(sgeometry.ExportToJson())
    sourceprj = osr.SpatialReference()
    sourceprj.ImportFromWkt(g.GetProjection())
    try:
        urn = sourceprj.ExportToXML().split('srsID')[1].split('=')[1].split('</gml:name>')[0].replace('"', '').replace('>', '')
        geojson = {                       
          "type": "FeatureCollection",                               
          "name": "AOI",          
          "crs": { "type": "name", "properties": { "name": urn} },
          "features": [           
          { "type": "Feature", "properties": {}, "geometry": spolygon}
                      ]                       
         }
    except:
        proj4 = sourceprj.ExportToProj4()
        geojson = {                       
          "type": "FeatureCollection",                               
          "name": "AOI",          
          "crs": { "type": "name", "properties": { "name": proj4, "type": "proj4"} },
          "features": [           
          { "type": "Feature", "properties": {}, "geometry": spolygon}
                      ]                       
         }

    geojson = json.dumps(geojson)

    if to_wgs84:
        sourceprj = osr.SpatialReference()
        sourceprj.ImportFromWkt(g.GetProjection())

        targetprj = osr.SpatialReference()
        targetprj.ImportFromProj4('+proj=longlat +datum=WGS84 +no_defs ')
        if int(gdal.VersionInfo()) > 2500000:
          targetprj.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
        transform = osr.CoordinateTransformation(sourceprj, targetprj)
        sgeometry.Transform(transform)

        # stgeometry = geometry.SimplifyPreserveTopology(abs(geom[1]) / 10. * 0.000001) 
        stpolygon = json.loads(sgeometry.ExportToJson())
        
        tgeojson = {                       
          "type": "FeatureCollection",                               
          "name": "AOI",          
          "crs": { "type": "name", "properties": { "name": "urn:ogc:def:crs:EPSG::4326"} },
          "features": [           
          { "type": "Feature", "properties": {}, "geometry": stpolygon}
                      ]                       
         }
        tgeojson = json.dumps(tgeojson)

        return tgeojson, geojson
    else:
        return geojson

def get_file_type(input_file):
    '''
    Check a file is vector or raster
    '''
    # from https://gis.stackexchange.com/questions/330215/gdal-2-x-method-for-determing-if-an-object-is-a-raster-or-a-vector
    ds = gdal.OpenEx(input_file)
    
    file_type = None
    if ds.RasterCount > 0 and ds.GetLayerCount() > 0:
        #  Raster and vector
        raise IOError('File contains both raster and vector')
    elif ds.RasterCount > 0:
        file_type = 'raster'
    elif ds.GetLayerCount() > 0:
        file_type = 'vector'
    else:
        raise IOError('Can only read raster/vector file included in gdal')
    
    return file_type
    
def sort_coordinates(list_of_xy_coords):
    # from https://pavcreations.com/clockwise-and-counterclockwise-sorting-of-coordinates/
    cx, cy = list_of_xy_coords.mean(0)
    x, y = list_of_xy_coords.T
    angles = np.arctan2(x-cx, y-cy)
    indices = np.argsort(-angles)
    return list_of_xy_coords[indices]

def cal_modis_tile(input_file):
    '''
    Check MODIS tile names based on the geometries from a raster or vector file.
    It will reproject the geometry to WGS84 coordinate system to filter
    the MODIS tiles to the ones intersecting with each geometry.
    Input:
        input_file: path or /viscurl/ link to the rater/vector file
    return:
        tiles: numpy array of string tiles in the format: h%02dv%02d, like h27v05
    '''
    
    file_type = get_file_type(input_file)
    if file_type == 'raster':
        aoi = get_boundary(input_file, to_wgs84=True)[0]
    else:
        aoi = input_file
    ds = ogr.Open(aoi)
    layer = ds.GetLayer(0)
    
    inSpatialRef = layer.GetSpatialRef()
    input_epsg_code = inSpatialRef.GetAttrValue('AUTHORITY',1)
    
    
    # g = ogr.Open('/vsicurl/vsizip/MODIS_tiles/tile_wgs.shp')
    # shapefile from https://github.com/yosukefk/modis_tile/blob/master/tile_wgs.zip
    g = ogr.Open('/vsizip/vsicurl/https://github.com/yosukefk/modis_tile/raw/master/tile_wgs.zip')
    modis_tile_layer = g.GetLayer(0)
    # output SpatialReference
    # outSpatialRef = osr.SpatialReference()
    # outSpatialRef.ImportFromEPSG(4326)
    outSpatialRef = modis_tile_layer.GetSpatialRef()
    output_epsg_code = outSpatialRef.GetAttrValue('AUTHORITY',1)
    
    tiles = []
    if input_epsg_code != output_epsg_code:
        for feat in layer:
            feat_geom = feat.geometry()
    
            # transform to the same coordinate system
            coordTrans = osr.CoordinateTransformation(inSpatialRef, outSpatialRef)
            feat_geom.Transform(coordTrans)
            
            modis_tile_layer.SetSpatialFilter(feat_geom)    
            tiles += [feature['tilename'] for feature in modis_tile_layer]
    else:
        for feat in layer:
            feat_geom = feat.geometry()
            modis_tile_layer.SetSpatialFilter(feat_geom)
            tiles += [feature['tilename'] for feature in modis_tile_layer]

def get_polygon_strs(input_file):  
    '''
        Turn polygon coordinates in the input file into counterclockwise convexhull strings
        Input file can be vector or raster file, as long as gdal can open it
    '''
    file_type = get_file_type(input_file)
    if file_type == 'raster':
        aoi = get_boundary(input_file, to_wgs84=True)[0]
    else:
        aoi = input_file
    ds = ogr.Open(aoi)
    layer = ds.GetLayer(0)
    
    inSpatialRef = layer.GetSpatialRef()
    input_epsg_code = inSpatialRef.GetAttrValue('AUTHORITY',1)
    
    
    # output SpatialReference
    outSpatialRef = osr.SpatialReference()
    outSpatialRef.ImportFromEPSG(4326)
    output_epsg_code = outSpatialRef.GetAttrValue('AUTHORITY',1)
    
    # tiles = []
    counterclockwise_polygon_strs = []
    if input_epsg_code != output_epsg_code:
        for feat in layer:
            feat_geom = feat.geometry()
    
            # transform to wgs84 coordinate system
            coordTrans = osr.CoordinateTransformation(inSpatialRef, outSpatialRef)
            feat_geom.Transform(coordTrans)
            
            coords = np.array(json.loads(feat_geom.ConvexHull().ExportToJson())['coordinates'][0])[:-1]
            coords = sort_coordinates(coords)
            coords = np.append(coords, [coords[0]], axis=0)
            polygon_str_counterclockwise = ','.join(coords.ravel().astype(str).tolist())
            counterclockwise_polygon_strs.append(polygon_str_counterclockwise)

    else:
        for feat in layer:
            feat_geom = feat.geometry()

            # feature geometry could be multipolygon, 
            # but only convexhull outline is needed
            # to find modis files
            coords = np.array(json.loads(feat_geom.ConvexHull().ExportToJson())['coordinates'][0])[:-1]
            coords = sort_coordinates(coords)
            coords = np.append(coords, [coords[0]], axis=0)
            polygon_str_counterclockwise = ','.join(coords.ravel().astype(str).tolist())
            counterclockwise_polygon_strs.append(polygon_str_counterclockwise)
    return counterclockwise_polygon_strs

@retry(tries=4, delay=1, backoff=2)
def query_VNP43MA1_temporal_spatial(polygon_str_counterclockwise, temporal_start, temporal_end):
    
    # temporal_filter = '%s,%s'%(temporal_start, temporal_end)
    
    temporal_filter = '%s,%s'%(temporal_start.strftime('%Y-%m-%dT%H:%M:%SZ'), temporal_end.strftime('%Y-%m-%dT%H:%M:%SZ'))
    
    url = 'https://cmr.earthdata.nasa.gov/search/granules.json'
    page_size = 2000
    response = requests.get(url,
                            params={
                                'short_name	': 'VNP43MA1',
                                'temporal': temporal_filter,
                                'page_size': page_size,
                                'provider': 'LPDAAC_ECS',
                                'version': '001',
                                'polygon': polygon_str_counterclockwise
                                    },
                                headers={ 'Accept': 'application/json'}
                    )

    filename_urls = []
    if response.ok:
        response_json = response.json()
        for i in response_json['feed']['entry']:
            
            fname = i['producer_granule_id']
            sensing_date = datetime.strptime(fname.split('.')[1], 'A%Y%j')
            if (sensing_date >= temporal_start) & (sensing_date <= temporal_end):
                filename_urls.append([i['producer_granule_id'], i['links'][0]['href']])
    

        # check for if there is more pages
        total_file_num = int(response.headers['CMR-Hits'])    
        if total_file_num > page_size:
            # print('There are multiple pages')
            num_pages = int(np.ceil(total_file_num / page_size))
            for i in range(1, num_pages):
                # print('Searching page %d of %d'%(i+1, num_pages))
                response = requests.get(url,
                            params={
                                'short_name	': 'VNP43MA1',
                                'equator_crossing_date': temporal_filter,
                                'page_size': page_size,
                                'page_num' : i+1,
                                'provider': 'LPDAAC_ECS',
                                'version': '001',
                                'polygon': polygon_str_counterclockwise
                                    },
                                headers={ 'Accept': 'application/json'}
                    )
                
                if response.ok:
                    response_json = response.json()
                    for i in response_json['feed']['entry']:
                        # filename_urls.append([i['producer_granule_id'], i['links'][0]['href']])
                        fname = i['producer_granule_id']
                        sensing_date = datetime.strptime(fname.split('.')[1], 'A%Y%j')
                        if (sensing_date >= temporal_start) & (sensing_date <= temporal_end):
                            filename_urls.append([i['producer_granule_id'], i['links'][0]['href']])
                
                        # print(i['id'], i['producer_granule_id'], i['links'][0]['href'])
                else:
                    raise ValueError(response.reason, response_json)
    else:
        raise ValueError(response.reason, response_json)

    return filename_urls

def find_files(aoi, obs_time, temporal_window = 16):
    counterclockwise_polygon_strs = get_polygon_strs(aoi)

    temporal_start = (obs_time - timedelta(days = int(temporal_window)))#.strftime('%Y-%m-%dT%H:%M:%SZ')
    temporal_end   = (obs_time + timedelta(days = int(temporal_window)))#.strftime('%Y-%m-%dT%H:%M:%SZ')

    # filename_urls = []
    # for polygon_str_counterclockwise in counterclockwise_polygon_strs[:2]:
    #     filename_urls += query_VNP43MA1_temporal_spatial(polygon_str_counterclockwise, temporal_filter)
    
    par = partial(query_VNP43MA1_temporal_spatial, temporal_start=temporal_start, temporal_end = temporal_end)
    with ThreadPoolExecutor(10) as executor:
        filename_urls = executor.map(par, counterclockwise_polygon_strs)
    
    filename_urls = list(filename_urls)
    filename_urls = np.concatenate(filename_urls)
    filename_urls = np.unique(filename_urls, axis=0)
    return filename_urls


@retry(tries=4, delay=1, backoff=2)
def downloader(url_fname, auth):
    url, fname = url_fname
    fname = os.path.abspath(fname)
    
    with requests.Session() as s:
        s.max_redirects = 100000
        s.auth = auth
        r1     = s.get(url)
        r      = s.get(r1.url, stream=True, headers={'user-agent': 'My app'})
        if r.ok:
            remote_size = int(r.headers['Content-Length'])
            if os.path.exists(fname):
                local_size = os.path.getsize(fname)
                if local_size != remote_size:
                    os.remove(fname)
                    data = r.content
                    if len(data) == remote_size:
                        with open(fname, 'wb') as f:
                            f.write(data)
                    else:
                        raise IOError('Failed to download the whole file.')
            else:
                data = r.content
                if len(data) == remote_size:
                    with open(fname, 'wb') as f:
                        f.write(data)
                else:           
                    raise IOError('Failed to download the whole file.')
        else:
            print(r.content)

def download_VNP43MA1(aoi, obs_time, VNP43_dir, temporal_window = 16):
    
    filename_urls = find_files(aoi, obs_time,  temporal_window = temporal_window)
    url_fnames_to_get = [[i[1], os.path.abspath(os.path.join(VNP43_dir, i[0]))] for i in filename_urls if not os.path.exists(os.path.join(VNP43_dir, i[0]))]
    if len(url_fnames_to_get) > 0:
        auth = get_auth()
        par = partial(downloader, auth = auth)
        with ThreadPoolExecutor(2) as executor:
            executor.map(par, url_fnames_to_get)
    filenames = [os.path.abspath(os.path.join(VNP43_dir, i[0])) for i in filename_urls]
    return filenames

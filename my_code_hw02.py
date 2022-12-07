# -- my_code_hw02.py
# -- hw02 GEO1015.2022
# -- [YOUR NAME] CHEN Mengying
# -- [YOUR STUDENT NUMBER] 5837324

from datetime import datetime
from pytz import timezone
import suncalc  # -- https://pypi.org/project/suncalc/
import pyproj
import math
import numpy
import rasterio
from rasterio import features


def is_sunny(dataset, px, py, dt):
    """
    !!! TO BE COMPLETED !!!
     
    Does the sun shine there at that moment?
     
    Input:
        dataset: the rasterio input dataset
        px:  x-coordinate of p
        py:  y-coordinate of p
        dt:  ISO-formatted datetime ('YYYY-MM-DD HH:MM'), eg '2022-08-12 13:32'
    Output:
        True/False: True = ground is illuminated; False = ground is not illuminated
           (raise Exception if p is outside the extent of the dataset)
           (raise Exception if p is no_data value)
    """
    # raise Exception("Point given is outside the extent of the dataset")
    # raise Exception("Point given has no_data value")
    '''
    第一步：用some_code_to_help_with_rasterio得到地面栅格每一个点的坐标和高度（nodata的赋值为0）
    第二步：根据输入的dt，用some_code_to_help_with_suncalc得到太阳的高度角，方位角,在转成地理坐标
    第三步：对点p（有坐标有高度，第一步得到），和太阳点，用bresenham_with_rasterio得到一个栅格线段
    第四步：比较栅格线段上的点对应的vp值和地面栅格高度值，若地面>vp，返回false，计算下一个地面点；若地面值<vp，继续，直到
    '''
    # step 1 get every points' (x,y,h)



def some_code_to_help_with_rasterio(dataset, px, py):
    """
    !!! USE THIS CODE !!!
     
    Example code that can be useful: use it for your function, 
    copy part of it, it's allowed.
    """
    # -- numpy of input
    n1 = dataset.read(1)
    # -- index of p in the numpy raster
    row, col = dataset.index(px, py)
    # -- an empty array with all values=0
    n2 = numpy.zeros(dataset.shape, dtype=numpy.int8)
    for i in range(n1.shape[0]):
        for j in range(n1.shape[1]):
            z = n1[i][j]
            if z == dataset.nodatavals:
                n2[i][j] = 50
    # -- put p with value=99
    n2[row, col] = 99
    # -- write this to disk
    output_file = 'testing.tiff'
    with rasterio.open(output_file, 'w',
                       driver='GTiff',
                       height=n2.shape[0],
                       width=n2.shape[1],
                       count=1,
                       dtype=rasterio.uint8,
                       crs=dataset.crs,
                       transform=dataset.transform) as dst:
        dst.write(n2.astype(rasterio.uint8), 1)
    print("File written to '%s'" % output_file)


def some_code_to_help_with_suncalc(dt, px, py):
    """
    !!! USE THIS CODE !!!
     
    Example code that can be useful: use it for your function, 
    copy part of it, it's allowed.
    """
    # -- define the timezone for where Delft is located
    ams_tz = timezone('Europe/Amsterdam')
    # -- get the datetime object for a given local time at a given date
    # dt = '2022-11-12 13:37'
    dto = ams_tz.localize(datetime.fromisoformat(dt))
    # print(dto)
    # -- now let's get the time UTC time, which must be used with suncalc
    # -- notice that this gives us either +01:00 in winter
    # -- and +02:00 because of EU summer time the local time
    time_utc = dto.astimezone(timezone('UTC'))
    print(time_utc)

    # -- get the position of the sun according to the input point p
    # -- to interpret the results https://github.com/mourner/suncalc#sun-position
    # -- suncalc.get_position(datetime_in_utc, latitude, longitude)
    pr = pyproj.Proj(init=4326)
    lon_p, lat_p = pr(px, py, inverse=True)  # get lon and lat of point p

    possun = suncalc.get_position(time_utc, lon_p, lat_p)
    breakpoint()
    # print(possun)
    if possun

    # -- convert from EPSG:28992 to EPSG:4326
    transfo = pyproj.Transformer.from_crs("EPSG:28992", "EPSG:4326")
    l = transfo.transform(81687.0, 447480.0)
    print(l)
    '''
    控制台得到四个值：
    dto（阿姆斯特丹时区时间） 
    time_utc（utc时间，用来给suncalc）
    possun（utc时间下太阳的altitude（高度，弧度制），azimuth（方位角弧度制）
    l(从28992转到4326之后的坐标)
    '''


def bresenham_with_rasterio(px, py, sx, sy):
    """
    !!! USE THIS CODE !!!
     
    Example code that can be useful: use it for your function, 
    copy part of it, it's allowed.
    """
    # -- https://rasterio.readthedocs.io/en/latest/topics/features.html#burning-shapes-into-a-raster
    # d = rasterio dataset as read in geo1015_hw02.py
    a = (px, py)
    b = (sx, sy)  # 这里应该输入两个坐标，一个是array里的p，一个是太阳坐标，坐标系都是当地坐标'''
    # -- create in-memory a simple GeoJSON LineString
    v = {}
    v["type"] = "LineString"
    v["coordinates"] = []
    v["coordinates"].append(d.xy(a[0], a[1]))
    v["coordinates"].append(d.xy(b[0], b[1]))  # 以线的格式，存两个端点的坐标，得到矢量线段'''
    shapes = [(v, 1)]
    re = features.rasterize(shapes,
                            out_shape=d.shape,
                            # all_touched=True,
                            transform=d.transform)
    # re is a numpy with d.shape where the line is rasterised (values != 0)
    '''最后得到了那条vp线段的栅格数组'''

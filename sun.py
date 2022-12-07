from datetime import datetime
from pytz import timezone
import suncalc  # -- https://pypi.org/project/suncalc/
import pyproj
import math
import numpy as np
import rasterio
from rasterio import features
import time


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
    第四步：比较栅格线段上的点对应的vp值和地面栅格高度值，若地面>vp，返回false，计算下一个地面点；若地面值<vp，继续，直到包含另一个端点
    '''
    # step 1 get point p's (x,y,h)
    n1 = dataset.read(1)
    row, col = dataset.index(px, py)  # get the row and col number of given point p
    try:
        height_p = n1[row, col]  # get the height value of given point p
    except Exception:
        raise Exception("Outside of the bounding box")
    if height_p == dataset.nodatavals:
        raise Exception("No data value for input point")

    # step 2 get sun's lan,lon
    # 1. change time to utc
    ams_tz = timezone('Europe/Amsterdam')
    dto = ams_tz.localize(datetime.fromisoformat(dt))
    time_utc = dto.astimezone(timezone('UTC'))

    # 2. get point p's lat and lon
    transfo = pyproj.Transformer.from_crs("EPSG:28992", "EPSG:4326")
    lat, lon = transfo.transform(px, py)
    # breakpoint()

    # 3. get the altitude and azimuth of the sun
    possun = suncalc.get_position(time_utc, lat=lat, lng=lon)

    # breakpoint()

    # 4. calculate the azimuth of p to the four corner
    bo = dataset.bounds
    a1 = math.atan2(bo.bottom - py, bo.left - px) + math.pi
    a2 = math.pi * 3 / 2 - math.atan2(bo.top - py, bo.left - px)
    a3 = 0 - (math.atan2(bo.top - py, bo.right - px) + math.pi / 2)
    a4 = 0 - (math.atan2(bo.bottom - py, bo.right - px) + math.pi / 2)
    # breakpoint()

    # 5. discuss azimuth on different edges to get sun's lat, lon
    if a4 <= possun.get('azimuth') < a1:
        sy = bo.bottom
        sx = px - math.atan(possun.get('azimuth')) * (py - sy)
    elif a1 <= possun.get('azimuth') < a2:
        sx = bo.left
        sy = py - px / math.atan(possun.get('azimuth'))
    elif a2 <= possun.get('azimuth') < math.pi or -math.pi < possun.get('azimuth') < a3:
        sy = bo.top
        sx = px + math.tan(possun.get('azimuth')) * (bo.top - py)
    elif a3 <= possun.get('azimuth') < a4:
        sx = bo.right
        sy = py - (bo.right - px) / math.tan(possun.get('azimuth'))

    # 6. turn sun's lat lon to EPSG:28992

    # step 3 get the intersected raster
    pr = (px, py)
    sr = (sx, sy)
    v = {}
    v["type"] = "LineString"
    v["coordinates"] = []
    v["coordinates"].append((pr[0], pr[1]))
    v["coordinates"].append((sr[0], sr[1]))
    shapes = [(v, 1)]
    re = features.rasterize(shapes,
                            out_shape=dataset.shape, transform=dataset.transform)
    # breakpoint()

    # step 4 compare the value on the line and the value on the ground item
    # 1. get the altitude of sun
    dist_sun_p = math.dist(pr, sr)
    sun_al = dist_sun_p * math.tan(possun.get('altitude'))

    # 2. get the height of point v ( v can be any point in the segment of p_sun)
    # and point g (g is the corresponding ground point to v)
    # breakpoint()
    line = np.multiply(re, n1)  # get the height value for raster_line
    n2 = np.where(line != dataset.nodatavals[0], line, -9999)  # no data points on the line get value -9999

    idx = np.where(re == 1)

    height_g = []
    for i in range(len(idx[0])):
        height_g.append(n2[idx[0][i]][idx[1][i]])
    # breakpoint()
    height_v = []
    for i in range(len(idx[0])):
        vh = height_p + abs(dataset.xy(idx[0][i], idx[1][i])[0] - px) * (sun_al - height_p) / (sx - px)
        height_v.append(vh)
    breakpoint()

    j = 0
    while j <= len(idx[0]):
        if height_v <= height_g:
            return False
        else:
            j = j + 1
            continue

    return True


def main():
    # -- this gives you a Rasterio dataset
    # -- https://rasterio.readthedocs.io/en/latest/quickstart.html
    d = rasterio.open('../ahn3_data/ahn3_dsm50cm_bk_small.tif')
    px, py, dt = 85278.83, 447036.12, '2022-08-12 07:32'
    start_time = time.time()
    re = is_sunny(d, px, py, dt)
    end_time = time.time()

    runtime = end_time - start_time

    print("The runtime of the project is: ", runtime, "seconds")
    print("YES it's sunny 😎") if re == True else print("NO it's not sunny 🔦")


main()

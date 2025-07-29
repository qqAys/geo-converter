import math

def gcj02_to_wgs84(gcj_lon, gcj_lat):
    """
    GCJ-02坐标转换为WGS-84坐标（近似反向转换）
    """
    # 迭代逼近法
    wgs_lon, wgs_lat = gcj_lon, gcj_lat
    for i in range(30):
        tmp_lon, tmp_lat = wgs84_to_gcj02(wgs_lon, wgs_lat)
        wgs_lon = wgs_lon - (tmp_lon - gcj_lon)
        wgs_lat = wgs_lat - (tmp_lat - gcj_lat)
    return wgs_lon, wgs_lat

def wgs84_to_gcj02(wgs_lon, wgs_lat):
    """
    WGS-84坐标转换为GCJ-02坐标（火星坐标系）
    
    Args:
        wgs_lon (float): WGS-84经度
        wgs_lat (float): WGS-84纬度
    
    Returns:
        tuple: (gcj_lon, gcj_lat) GCJ-02坐标
    """
    
    # 判断是否在中国范围内
    def out_of_china(lon, lat):
        if lon < 72.004 or lon > 137.8347:
            return True
        if lat < 0.8293 or lat > 55.8271:
            return True
        return False
    
    # 如果不在中国范围内，直接返回原坐标
    if out_of_china(wgs_lon, wgs_lat):
        return wgs_lon, wgs_lat
    
    # 转换参数
    a = 6378245.0  # 长半轴
    ee = 0.00669342162296594323  # 扁率
    
    def transform_lat(x, y):
        """纬度转换"""
        ret = -100.0 + 2.0 * x + 3.0 * y + 0.2 * y * y + 0.1 * x * y + 0.2 * math.sqrt(abs(x))
        ret += (20.0 * math.sin(6.0 * x * math.pi) + 20.0 * math.sin(2.0 * x * math.pi)) * 2.0 / 3.0
        ret += (20.0 * math.sin(y * math.pi) + 40.0 * math.sin(y / 3.0 * math.pi)) * 2.0 / 3.0
        ret += (160.0 * math.sin(y / 12.0 * math.pi) + 320 * math.sin(y * math.pi / 30.0)) * 2.0 / 3.0
        return ret
    
    def transform_lon(x, y):
        """经度转换"""
        ret = 300.0 + x + 2.0 * y + 0.1 * x * x + 0.1 * x * y + 0.1 * math.sqrt(abs(x))
        ret += (20.0 * math.sin(6.0 * x * math.pi) + 20.0 * math.sin(2.0 * x * math.pi)) * 2.0 / 3.0
        ret += (20.0 * math.sin(x * math.pi) + 40.0 * math.sin(x / 3.0 * math.pi)) * 2.0 / 3.0
        ret += (150.0 * math.sin(x / 12.0 * math.pi) + 300.0 * math.sin(x / 30.0 * math.pi)) * 2.0 / 3.0
        return ret
    
    # 计算偏移量
    d_lat = transform_lat(wgs_lon - 105.0, wgs_lat - 35.0)
    d_lon = transform_lon(wgs_lon - 105.0, wgs_lat - 35.0)
    
    # 纬度弧度
    rad_lat = wgs_lat / 180.0 * math.pi
    magic = math.sin(rad_lat)
    magic = 1 - ee * magic * magic
    sqrt_magic = math.sqrt(magic)
    
    # 调整偏移量
    d_lat = (d_lat * 180.0) / ((a * (1 - ee)) / (magic * sqrt_magic) * math.pi)
    d_lon = (d_lon * 180.0) / (a / sqrt_magic * math.cos(rad_lat) * math.pi)
    
    # 计算最终坐标
    gcj_lat = wgs_lat + d_lat
    gcj_lon = wgs_lon + d_lon
    
    return gcj_lon, gcj_lat


if __name__ == "__main__":
    wgs_lon = 116.397428  # 北京天安门经度
    wgs_lat = 39.90923   # 北京天安门纬度

    print(f"WGS-84坐标: ({wgs_lon}, {wgs_lat})")
    
    gcj_lon, gcj_lat = wgs84_to_gcj02(wgs_lon, wgs_lat)
    
    print(f"GCJ-02坐标: ({gcj_lon}, {gcj_lat})")
    print(f"偏移量: ({gcj_lon - wgs_lon}, {gcj_lat - wgs_lat})")

    wgs_lon, wgs_lat = gcj02_to_wgs84(gcj_lon, gcj_lat)

    print(f"WGS-84坐标(转换): ({wgs_lon}, {wgs_lat})")

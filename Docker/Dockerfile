FROM osgeo/gdal:ubuntu-small-3.6.0 as base


WORKDIR /app
COPY run_siac.py /app/run_siac.py
RUN set -xe \
   && apt-get update \
    && apt-get install -y --no-install-recommends python3-pip libgomp1 \
    && pip install --no-cache-dir wheel lightgbm scikit-learn \
    numba scikit-image pyproj requests scipy psutil six retry \
    earthengine-api pyproj https://github.com/MarcYin/SIAC/archive/master.zip \
    && rm -rf /var/lib/apt/lists/*



CMD [ "python", "run_siac.py"]


# docker run -it --name SIAC -v ~/Downloads/MCD43:/app/MCD43 \
#                           -v ~/Downloads/CAMS:/app/CAMS \
#                           -v ~/Downloads/DEM:/app/DEM \
#                           -v ~/Downloads/water_mask:/app/water_mask \
#                           -v "$HOME"/.config/earthengine:/root/.config/earthengine \
#                           -v ~/Downloads/S2_L1C:/app/S2_L1C \
#                           -v ~/Downloads/S2_L2A:/app/S2_L2A \
#                           --rm marcyin/siac_app 

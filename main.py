import gpxpy.gpx
import numpy as np
import matplotlib.pyplot as plt
import geopandas as gpd
from shapely.geometry import Point, LineString
from geopy.distance import geodesic

# Function to calculate bearing between two points
def calculate_initial_compass_bearing(pointA, pointB):
    lat1 = np.deg2rad(pointA.latitude)
    lat2 = np.deg2rad(pointB.latitude)
    diffLong = np.deg2rad(pointB.longitude - pointA.longitude)

    x = np.sin(diffLong) * np.cos(lat2)
    y = np.cos(lat1) * np.sin(lat2) - (np.sin(lat1) * np.cos(lat2) * np.cos(diffLong))

    initial_bearing = np.arctan2(x, y)

    # Now we have the initial bearing but math.atan2 return values
    # from -180° to +180° which is not what we want for a compass bearing
    initial_bearing = np.rad2deg(initial_bearing)
    compass_bearing = (initial_bearing + 360) % 360

    return compass_bearing

# List of GPX files for each boat and corresponding boat names
gpx_files = [
    'tracks/thora-riga-salacgriva.gpx',
    'tracks/nolax-riga-salacgriva.gpx',
    'tracks/shadow-riga-salacgriva.gpx',
]

boat_names = ['Thora', 'Nola X', 'Shadow']

# Wind and course information
wind_data = [
    {'leg': 'Leg 1', 'length_nm': 2.67, 'wind_angle': 230, 'course_angle': 113},
    {'leg': 'Leg 2', 'length_nm': 2.39, 'wind_angle': 230, 'course_angle': 261},
    {'leg': 'Leg 3', 'length_nm': 2.39, 'wind_angle': 245, 'course_angle': 81},
    {'leg': 'Leg 4', 'length_nm': 2.76, 'wind_angle': 245, 'course_angle': 310},
    {'leg': 'Leg 5', 'length_nm': 2.76, 'wind_angle': 245, 'course_angle': 130},
    {'leg': 'Leg 6', 'length_nm': 2.67, 'wind_angle': 245, 'course_angle': 293},
]

# Initialize storage for all boats' data
all_boats_segments = []

# Function to parse and segment GPX data
def process_gpx(gpx_file):
    with open(gpx_file, 'r') as file:
        gpx = gpxpy.parse(file)
    
    track_points = []
    for track in gpx.tracks:
        for segment in track.segments:
            for point in segment.points:
                track_points.append(point)
    
    bearings = []
    distances = []
    for i in range(1, len(track_points)):
        bearing = calculate_initial_compass_bearing(track_points[i-1], track_points[i])
        distance = geodesic((track_points[i-1].latitude, track_points[i-1].longitude),
                            (track_points[i].latitude, track_points[i].longitude)).nautical
        bearings.append(bearing)
        distances.append(distance)

    # Segment track points based on wind data
    segments = [[] for _ in range(len(wind_data))]
    current_leg = 0
    current_leg_distance = 0

    for i in range(1, len(track_points)):
        if current_leg < len(wind_data):
            segments[current_leg].append(track_points[i])
            current_leg_distance += distances[i-1]

            if current_leg_distance >= wind_data[current_leg]['length_nm']:
                current_leg += 1
                current_leg_distance = 0
    
    return segments

# Function to convert track points to GeoDataFrame
def points_to_gdf(points):
    return gpd.GeoDataFrame(
        {'geometry': [Point(p.longitude, p.latitude) for p in points]},
        crs="EPSG:4326"
    )

# Process each GPX file and store segments
for idx, gpx_file in enumerate(gpx_files):
    segments = process_gpx(gpx_file)
    all_boats_segments.append({
        'name': boat_names[idx],
        'segments': segments
    })

# Visualize the tracks of all boats
fig, ax = plt.subplots(figsize=(10, 10))
colors = ['b', 'g', 'r']  # Color list for different boats

for boat_idx, boat_data in enumerate(all_boats_segments):
    boat_segments = boat_data['segments']
    full_track = []
    for segment in boat_segments:
        full_track.extend(segment)  # Combine all segments into a single track
    gdf = points_to_gdf(full_track)
    gdf.plot(ax=ax, color=colors[boat_idx % len(colors)], markersize=2, label=boat_data['name'])

plt.legend()
plt.title("Sailing Track Comparison")
plt.xlabel("Longitude")
plt.ylabel("Latitude")
plt.show()
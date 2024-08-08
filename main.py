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

# Load GPX file
with open('tracks/shadow-riga-salacgriva.gpx', 'r') as gpx_file:
    gpx = gpxpy.parse(gpx_file)

# Extract track points
track_points = []
for track in gpx.tracks:
    for segment in track.segments:
        for point in segment.points:
            track_points.append(point)

# Calculate bearings and distances
bearings = []
distances = []
for i in range(1, len(track_points)):
    bearing = calculate_initial_compass_bearing(track_points[i-1], track_points[i])
    distance = geodesic((track_points[i-1].latitude, track_points[i-1].longitude),
                        (track_points[i].latitude, track_points[i].longitude)).nautical
    bearings.append(bearing)
    distances.append(distance)

# Assign track points to legs
wind_data = [
    {'leg': 'Leg 1', 'length_nm': 2.67, 'wind_angle': 230, 'course_angle': 113},
    {'leg': 'Leg 2', 'length_nm': 2.39, 'wind_angle': 230, 'course_angle': 261},
    {'leg': 'Leg 3', 'length_nm': 2.39, 'wind_angle': 245, 'course_angle': 81},
    {'leg': 'Leg 4', 'length_nm': 2.76, 'wind_angle': 245, 'course_angle': 310},
    {'leg': 'Leg 5', 'length_nm': 2.76, 'wind_angle': 245, 'course_angle': 130},
    {'leg': 'Leg 6', 'length_nm': 2.67, 'wind_angle': 245, 'course_angle': 293},
]

segments = [[] for _ in range(len(wind_data))]
current_leg = 0
current_leg_distance = 0

for i in range(1, len(track_points)):
    if current_leg < len(wind_data):
        segments[current_leg].append(track_points[i])
        current_leg_distance += distances[i-1]
        
        # Check if we reached the length of the current leg
        if current_leg_distance >= wind_data[current_leg]['length_nm']:
            current_leg += 1
            current_leg_distance = 0

# Print segment information
for i, segment in enumerate(segments):
    print(f"Leg {i+1} (Length: {wind_data[i]['length_nm']} nm, Wind Angle: {wind_data[i]['wind_angle']}, Course Angle: {wind_data[i]['course_angle']}): {len(segment)} track points")

# Visualize segments
# Convert track points to GeoDataFrame
def points_to_gdf(points):
    return gpd.GeoDataFrame(
        {'geometry': [Point(p.longitude, p.latitude) for p in points]},
        crs="EPSG:4326"
    )

# Create GeoDataFrame for each segment
gdfs = [points_to_gdf(segment) for segment in segments]

# Plot segments
fig, ax = plt.subplots(figsize=(10, 10))
for i, gdf in enumerate(gdfs):
    gdf.plot(ax=ax, label=f"Leg {i+1}", markersize=2)

plt.legend()
plt.title("Sailing Track Segmentation")
plt.xlabel("Longitude")
plt.ylabel("Latitude")
plt.show()
#!/usr/bin/env python3

import geopandas as gpd
from shapely.geometry import Point
from shapely.validation import explain_validity
import re


# ---------------------------------------------------------------------
# Input / output paths
# ---------------------------------------------------------------------

src_path = "./coast1_polygons.shp"

poly_out = "./atlantic_polygons_diagnostic.shp"
points_out = "./atlantic_error_points.shp"

report_csv_out = "./atlantic_invalid_geometry_report.csv"

# Your source file is NAD83 geographic coordinates.
# If another source file uses a different CRS, update this value.
default_crs = "EPSG:4269"


# ---------------------------------------------------------------------
# Helper function: get point at actual validity error location if possible
# ---------------------------------------------------------------------

def get_validity_error_point(geom):
    """
    Try to extract the actual geometry-validity error coordinate from Shapely/GEOS.

    For invalid polygons, explain_validity() often returns text like:
        Self-intersection[-70.123456 42.123456]
        Ring Self-intersection[-71.234567 43.234567]

    If a coordinate is found, return that point.

    If no coordinate is found, fall back to the bounding-box center.
    """

    reason = explain_validity(geom)

    # Look for coordinates inside square brackets: [x y]
    match = re.search(r"\[([^\]]+)\]", reason)

    if match:
        coord_text = match.group(1)

        # Extract numeric values from the coordinate text
        nums = re.findall(
            r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?",
            coord_text
        )

        if len(nums) >= 2:
            x = float(nums[0])
            y = float(nums[1])
            return Point(x, y), reason, "validity_location"

    # Fallback if Shapely does not provide an exact error coordinate
    bounds = geom.bounds

    cx = (bounds[0] + bounds[2]) / 2.0
    cy = (bounds[1] + bounds[3]) / 2.0

    return Point(cx, cy), reason, "bbox_center_fallback"


# ---------------------------------------------------------------------
# Read source data
# ---------------------------------------------------------------------

print(f"Reading original shapefile: {src_path}")

gdf = gpd.read_file(src_path)

# If the source CRS is missing, assign the expected CRS.
# If it already has a CRS, keep it.
if gdf.crs is None:
    print(f"Input CRS is missing. Assigning default CRS: {default_crs}")
    gdf = gdf.set_crs(default_crs)
else:
    print(f"Input CRS found: {gdf.crs}")

print("Input bounds:", gdf.total_bounds)


# ---------------------------------------------------------------------
# Preserve original row IDs before exploding multipart geometries
# ---------------------------------------------------------------------

gdf["src_row"] = gdf.index

# Preserve common GIS ID fields if they exist
possible_id_fields = ["OBJECTID", "ObjectID", "FID", "fid", "ID", "id"]

existing_id_fields = [
    field for field in possible_id_fields
    if field in gdf.columns
]

if existing_id_fields:
    print("Existing ID fields found:", existing_id_fields)
else:
    print("No common OBJECTID/FID field found. Using source row index only.")


# ---------------------------------------------------------------------
# Explode multipart polygons into individual features
# ---------------------------------------------------------------------

print("\nExploding multipart geometries...")

gdf = gdf.explode(index_parts=False).reset_index(drop=True)

# Store row number after explode
gdf["expl_row"] = gdf.index

print("Total features after explode:", len(gdf))


# ---------------------------------------------------------------------
# Detect invalid geometries automatically
# ---------------------------------------------------------------------

print("\nChecking geometry validity...")

# Handle null and empty geometries separately
gdf["geom_null"] = gdf.geometry.isna()
gdf["geom_empty"] = gdf.geometry.is_empty

# For normal non-null/non-empty geometries, check validity
valid_check_mask = (~gdf["geom_null"]) & (~gdf["geom_empty"])

gdf["is_valid"] = True
gdf.loc[valid_check_mask, "is_valid"] = gdf.loc[
    valid_check_mask,
    "geometry"
].is_valid

# Invalid geometries are non-null, non-empty, and not valid
invalid_mask = valid_check_mask & (~gdf["is_valid"])

# Create invalid geometry subset
gdf_errors = gdf.loc[invalid_mask].copy()

print("Total invalid geometries found:", len(gdf_errors))

if len(gdf_errors) == 0:
    print("\nNo invalid geometries were detected by Shapely/GEOS.")
    print("No diagnostic shapefiles were created.")
    print("If ArcGIS reports errors but this script finds none, the tools may be using different validation rules.")
    raise SystemExit(0)


# ---------------------------------------------------------------------
# Add validity reasons
# ---------------------------------------------------------------------

gdf_errors["err_reason"] = gdf_errors.geometry.apply(explain_validity)
gdf_errors["is_invalid"] = 1

# Make sure CRS is preserved
gdf_errors = gdf_errors.set_crs(gdf.crs, allow_override=True)

print("\nInvalid geometry reason counts:")
print(gdf_errors["err_reason"].value_counts())

print("\nError polygon CRS:", gdf_errors.crs)
print("Error polygon bounds:", gdf_errors.total_bounds)


# ---------------------------------------------------------------------
# Save CSV report
# ---------------------------------------------------------------------

print(f"\nSaving invalid geometry report to: {report_csv_out}")

report_fields = ["src_row", "expl_row", "is_valid", "err_reason"]

# Include existing OBJECTID/FID fields if available
report_fields = existing_id_fields + report_fields

gdf_errors[report_fields].to_csv(report_csv_out, index=False)


# ---------------------------------------------------------------------
# Save invalid polygons
# ---------------------------------------------------------------------

print(f"\nSaving invalid polygon features to: {poly_out}")

# Note: Shapefile field names are limited to 10 characters.
# GeoPandas/Fiona may truncate long field names automatically.
gdf_errors.to_file(poly_out, driver="ESRI Shapefile")


# ---------------------------------------------------------------------
# Generate point markers at actual validity-error coordinates when possible
# ---------------------------------------------------------------------

print("\nGenerating error-location points...")

error_points = []

for idx, row in gdf_errors.iterrows():
    pt, reason, point_method = get_validity_error_point(row.geometry)

    error_num = len(error_points) + 1

    print(
        f"Feature {error_num:02d}, "
        f"source row {row['src_row']}, "
        f"exploded row {row['expl_row']}: "
        f"{reason} | point method: {point_method} | "
        f"point: ({pt.x}, {pt.y})"
    )

    point_record = {
        "geometry": pt,
        "err_num": error_num,
        "src_row": int(row["src_row"]),
        "expl_row": int(row["expl_row"]),
        "method": point_method[:50],
        "err_reason": reason[:250]
    }

    # Also carry over OBJECTID/FID fields if available
    for field in existing_id_fields:
        point_record[field] = row[field]

    error_points.append(point_record)


points_gdf = gpd.GeoDataFrame(error_points, crs=gdf.crs)

print("\nPoint CRS:", points_gdf.crs)
print("Point bounds:", points_gdf.total_bounds)

print(f"Saving error points to: {points_out}")
points_gdf.to_file(points_out, driver="ESRI Shapefile")


# ---------------------------------------------------------------------
# Final message
# ---------------------------------------------------------------------

print("\nSUCCESS:")
print(f"  Invalid geometries found: {len(gdf_errors)}")
print(f"  Polygon shapefile:        {poly_out}")
print(f"  Point shapefile:          {points_out}")
print(f"  CSV report:               {report_csv_out}")

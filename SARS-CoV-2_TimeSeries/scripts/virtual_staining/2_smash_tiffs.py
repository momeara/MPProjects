
import typing
import glob
import re
import pandas as pd
import tifffile
import tqdm

    
def parse_image_paths(
        image_paths : [str],
        regex : str,
        field_types : typing.Optional[typing.Dict[str,type]] = None):
    """
    Parse image paths into a pandas data frame

    image_paths : file paths to images
    regex : a regular expression to extract fields
    field_types : a dictionary giving the types of the fields
    """    
    rows = []
    for image_path in image_paths:
        row = {"image_path" : image_path}
        m = re.match(regex, image_path)
        row.update(m.groupdict())
        rows.append(row)
    image_info = pd.DataFrame(rows)
    if field_types is not None:
        image_info = image_info.astype(field_types)
    return image_info

def merge_tiffs(
        image_info,
        image_path_field,
        group_by,
        output_template):
    for group_name, image_info_group in tqdm.tqdm(image_info.groupby(group_by)):
        with tifffile.TiffWriter(
                output_template.format(
                    **{field : value for field, value in zip(group_by, group_name)}),
                bigtiff=True) as stack:
            for one_image_info in image_info_group.to_dict(orient = "records"):
                stack.save(
                    tifffile.imread(one_image_info['image_path']),
                    photometric='minisblack', 
                    contiguous=True)

image_info = parse_image_paths(
    image_paths = glob.glob("raw_data/virtual_staining/20200531T160447_48_hour_DPC/*tif"),
    regex = r".*W(?P<well_id>\d{4})F(?P<field_id>\d{4})T(?P<time_point>\d{4})Z(?P<z_plane>\d{3})C(?P<channel>\d{1}).tif$",
    field_types = {
        "well_id" : int,
        "field_id" : int,
        "time_point" : int,
        "z_plane" : int,
        "channel" : int})

merge_tiffs(
    image_info = image_info,
    image_path_field = 'image_path',
    group_by = ["well_id", "field_id", "time_point", "channel"],
    output_template = "raw_data/virtual_staining/20200531T160447_48_hour_DPC_merged/W{well_id:04d}F{field_id:04d}T{time_point:04d}C{channel:01d}.tif")


#######
image_info = parse_image_paths(
    image_paths = glob.glob("raw_data/virtual_staining/20200530T160342_48_hour/*tif"),
    regex = r".*W(?P<well_id>\d{4})F(?P<field_id>\d{4})T(?P<time_point>\d{4})Z(?P<z_plane>\d{3})C(?P<channel>\d{1}).tif$",
    field_types = {
        "well_id" : int,
        "field_id" : int,
        "time_point" : int,
        "z_plane" : int,
        "channel" : int})

merge_tiffs(
    image_info = image_info,
    image_path_field = 'image_path',
    group_by = ["well_id", "field_id", "time_point", "channel"],
    output_template = "raw_data/virtual_staining/20200530T160342_48_hour_merged/W{well_id:04d}F{field_id:04d}T{time_point:04d}C{channel:01d}.tif")





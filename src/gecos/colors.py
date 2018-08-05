import numpy as np
from colormath.color_objects import LabColor, XYZColor, sRGBColor
from colormath.color_conversions import convert_color


def convert_lab_to_rgb(lab):
    lab_color = LabColor(lab[0], lab[1], lab[2])
    rgb_color = convert_color(lab_color, sRGBColor)
    if          rgb_color.rgb_r > 0 and rgb_color.rgb_r < 1 \
            and rgb_color.rgb_g > 0 and rgb_color.rgb_g < 1 \
            and rgb_color.rgb_b > 0 and rgb_color.rgb_b < 1:
                return np.array(
                    [rgb_color.rgb_r, rgb_color.rgb_g, rgb_color.rgb_b]
                )
    else:
                return np.full(3, np.nan)

convert_lab_to_rgb = np.vectorize(convert_lab_to_rgb, signature="(3)->(3)")

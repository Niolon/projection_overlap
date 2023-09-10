
from collections import OrderedDict

import pandas as pd
import numpy as np
import fractions
import os
import re
import warnings
import shapely
from itertools import product
from configparser import ConfigParser
import logging
import argparse
from typing import Any, Dict, Union, Tuple, List

log = logging.getLogger('logger')
log.setLevel(logging.DEBUG)

formatter = logging.Formatter('%(message)s')

config_template = """[cif]
; Full path to the cif file you want to use
path = 1-SP-pyr.cif

; label of the dataset within the cif file. An integer number can be used to
; select the first (0) or second(1) dataset. Interpretation as a string name
; will always be tried first (in case there is a dataset named 1)
dataset_label = 0

; First of your reference molecules, this is the molecule that determines the
; reference plane.
; Beginning is a symmetry code in cif convention then a colon and then atom names
; For the atom names, the ordering matters! Edges of the polygon are between
; neighbouring atoms and the final and first polygon. If multiple symmetry codes
; are needed these can be separated by a semicolon
molecule1 = x, y, z:C1 C2 C3 C4 C5 C6 C7 C8 C9 C10 C11 C12 C13 C14

; Base geometry of the second molecule. Can be used as given or potential
; overlaps can be searched by the script. Syntax is the same as molecule1
molecule2 = 1+x, 1+y, +z:C1 C2 C3 C4 C5 C6 C7 C8 C9 C10 C11 C12 C13 C14

[geometry]
; direction1 is determined by placing the the mean plane of molecule1 into the
; xy plane. Direction2 can be chosen to get additional information about
; offsets in the xy plane. Before the colon is a directional vector in the xy
; plane. After the colon the script expects two sets of positions separated by
; a comma. All values within a set will be averaged. Three options are available
; for the definition
; position(0.1 0.2 0.3) is a position in cartesian coordinates
; centroid(1)/centroid(2) are the centroids of molecule1 and 2 respectively
; atom(C2) is an atom position of an atom which is available with the symmetry
; x, y, z in either molecule1 or molecule2

direction2 = 1 0:centroid(1), atom(C2) atom(C3)

; This chooses what is the divisor of the determined overlap molecule. Options
; are: Molecule1, Molecule2: determined area of Molecule1 or two. As the
; plane of projection is determined from Molecule1, this is the more sensible
; choice. Both: The reference polygon is the Union of projected molecule1 and
; projected molecule2 polygons
reference = Molecule1

[search]
; Used to automatically find the symmetry modified molecule2 with the largest
; overlap to molecule1 within the specified parameters. Translation symmetry
; is included in the search

; Do a search for molecule 2(yes) or take positions as given (no)
molecule2_search = yes

; no interatomic distance between molecule can be under this value. Used to
; exclude self overlap and disorder
min_distance = 0.7

; minimum interatomic distance of any two atoms within the two molecules must
; be smaller than this value
max_distance = 5.0

; Maximum angle between mean planes of molecule1 and found molecule2, use to
; exclude T-shaped interactions
max_angle_diff = 45.0

; Minimum overlap to be considered for the search
min_overlap_ratio = 0.01

; Number of cells in x,y and z direction to be searched for tranlational
; symmetry
n_neighbour_cells = 3

[figure]
; Colour of the overlapping region between the two molecules
overlap_colour = #001242

; Colour of the remaining area used as divisor but not in overlap
reference_colour = #7EA8BE

; Colour of the edges of both polygons/molecules
edge_colour = #a0a0a0

; Filename of the graphics output file. See matplotlib for all file options.
; Useful ones should be .png, .pdf and .svg
output_file = overlap.png

; output dpi for the graphics used. Used for png export
output_dpi = 300"""

def ciflike_to_dict(
    cif_fo: str,
    return_descr: Union[str, int, None] = None,
    resolve_esd: bool =True
) -> Union[Dict[str, Dict], Dict]:
    """Function to read in cif or cif-like (e.g. fcf) files. Can return all
    structures contained in the file. To return only the data of one structure
    use return_descr to select by dataset name or index.

    Parameters
    ----------
    filename : str
        Path to a cif or cif-like file
    return_descr : Union[str, int, None], optional
        Can be used to only return a specific dataset from the cif file
        if a string is given as argument, the dataset with that name is
        returned. An integer will return the dataset by index (i.e. 0 will give
        the first dataset in the file), None will return a dict with all
        datasets, with the dataset name as key, by default None
    resolve_esd : bool, optional
        If this argument is set to true, will split arguments, which have an
        esd into two arguments: arg, arg_esd, False will return a string in this
        case (i.e. '12(3)'), by default True

    Returns
    -------
    cif_content: Union[Dict[str, Dict], Dict]
        Returns a dictionary of dataset_name, dataset_dict pairs or a single
        dataset as OrderedDict. Within the dataset all entries are given as
        key, value pairs with the key being the entry in the cif_file without
        the preceding underscore. The routine will try to cast the value into
        a float, int or string, depending on whether digits, digits and a dot
        or other characters are present. All loops are given as a list of pandas
        DataFrame objects under the 'loops' keyword.

    Raises
    ------
    e
        If an exception occurs it is raised after printing the line in the cif
        file for debugging purposes
    ValueError
        The return_descr was in an invalid type
    """
    PATTERN = re.compile(r'''((?:[^ "']|"[^"]*"|'[^']*')+)''')
    lines = [line for line in cif_fo.readlines()]
    datablocks = OrderedDict()
    current_loop_lines = []
    current_loop_titles = []
    # If there is data before the first data entry store it as preblock
    current_block = 'preblock'
    in_loop = False
    in_loop_titles = False
    in_multiline = False
    multiline_title = 'InvalidTitle' # This should never be used
    multiline_entries = []
    current_line_collect = []
    try:
        for index, raw_line in enumerate(lines):
            line = raw_line.strip().lstrip()
            if len(line.strip()) == 0 or line.startswith('#'):
                # empty or comment line
                continue
            if in_loop and not in_loop_titles and (line.startswith('_') or line.startswith('loop_')):
                # The current loop has ended append entries as new DataFrame
                in_loop = False
                if len(current_loop_lines) > 0:
                    new_df = pd.DataFrame(current_loop_lines)
                    for key in new_df:
                        new_df[key] = pd.to_numeric(new_df[key], errors='ignore')
                    if resolve_esd:
                        for column in new_df.columns:
                            if new_df[column].dtype != 'O':
                                continue
                            concatenate = ''.join(new_df[column])
                            if  re.search(r'[\(\)]', concatenate) is not None and re.search(r'[^\d^\.^\(^\)\-\+]', concatenate) is None:
                                values, errors = np.array([split_error(val) for val in new_df[column]]).T
                                new_df[column] = values
                                new_df[column+'_esd'] = errors
                    datablocks[current_block]['loops'].append(new_df)
                # empty all stored entries
                current_loop_lines = []
                current_loop_titles = []
                current_line_collect = []
            if line.startswith('data_'):
                # New data block
                current_block = line[5:]
                if current_block not in datablocks:
                    datablocks[current_block] = OrderedDict([('loops', [])])
            elif line.startswith('loop_'):
                # a new loop / table starts
                in_loop = True
                in_loop_titles = True
            elif in_loop and in_loop_titles and line.startswith('_'):
                # This line is a title entry within a loop
                current_loop_titles.append(line[1:])
            elif in_loop:
                # This line contains data within a loop
                in_loop_titles = False
                line_split = [item.strip() for item in PATTERN.split(line) if item != '' and not item.isspace()]
                line_split = [item[1:-1] if "'" in item else item for item in line_split]
                current_line_collect += line_split
                if len(current_line_collect) == len(current_loop_titles):
                    current_loop_lines.append(OrderedDict())
                    for index2, item in enumerate(current_line_collect):
                        current_loop_lines[-1][current_loop_titles[index2]] = item
                    current_line_collect = []
            elif line.startswith('_'):
                # we are not in a loop -> single line or multiline string entry
                line_split = [item.strip() for item in PATTERN.split(line) if item != '' and not item.isspace()]
                line_split = [item[1:-1] if "'" in item else item for item in line_split]
                if len(line_split) > 1:
                    if resolve_esd:
                        test = line_split[1]
                        if len(test) == 0:
                            datablocks[current_block][line_split[0][1:]] = None
                        elif (re.search(r'[^\d]', test) is None):
                            datablocks[current_block][line_split[0][1:]] = int(test)
                        elif re.search(r'[^\d^\.]', test) is None and re.search(r'\d', test) is not None:
                            datablocks[current_block][line_split[0][1:]] = float(test)
                        elif re.search(r'[\(\)]', test) is not None and re.search(r'[^\d^\.^\(^\)\-\+]', test) is None:
                            val, error = split_error(test)
                            datablocks[current_block][line_split[0][1:]] = val
                            datablocks[current_block][line_split[0][1:] + '_esd'] = error
                        elif test.startswith('-'):
                            # This accounts for negative values without also catching dates
                            if (re.search(r'[^\d]', test[1:]) is None):
                                datablocks[current_block][line_split[0][1:]] = int(test)
                            elif re.search(r'[^\-^\d^\.]', test[1:]) is None and re.search(r'\d', test[1:]) is not None:
                                datablocks[current_block][line_split[0][1:]] = float(test)
                            else:
                                datablocks[current_block][line_split[0][1:]] = line_split[1]
                        elif test == '?':
                            datablocks[current_block][line_split[0][1:]] = None
                        else:
                            datablocks[current_block][line_split[0][1:]] = line_split[1]
                    else:
                        datablocks[current_block][line_split[0][1:]] = line_split[1]
                else:
                    multiline_title = line_split[0][1:]
            elif line.startswith(';') and in_multiline:
                datablocks[current_block][multiline_title] = '\n'.join(multiline_entries)
                multiline_entries = []
                in_multiline = False
            elif line.startswith(';') and not in_multiline:
                in_multiline = True
            elif in_multiline:
                multiline_entries.append(line)
    except Exception as e:
        print('Error in Line {index}')
        print(line)
        raise e

    # We might have a final loop
    if in_loop:
        in_loop = False
        if len(current_loop_lines) > 0:
            new_df = pd.DataFrame(current_loop_lines)
            for key in new_df:
                new_df[key] = pd.to_numeric(new_df[key], errors='ignore')
            if resolve_esd:
                for column in new_df.columns:
                    if new_df[column].dtype != 'O':
                        continue
                    concatenate = ''.join(new_df[column])
                    if  re.search(r'[\(\)]', concatenate) is not None and re.search(r'[^\d^\.^\(^\)\-\+]', concatenate) is None:
                        values, errors = np.array([split_error(val) for val in new_df[column]]).T
                        new_df[column] = values
                        new_df[column+'_esd'] = errors
            datablocks[current_block]['loops'].append(new_df)
    if return_descr is None:
        return datablocks
    elif type(return_descr) is int:
        return datablocks[list(datablocks.keys())[return_descr]]
    elif type(return_descr) is str:
        return datablocks[return_descr]
    else:
        raise ValueError('Invalid return_descr value. Must be either None, index as int or name as str')

def split_error(string: str) -> Union[Tuple[float, float], Tuple[int, int]]:
    """Helper function to split a string containing a value with error in
    brackets to a value-esd pair

    Parameters
    ----------
    string : str
        Input string containing the value to be split

    Returns
    -------
    Union[Tuple[float, float], Tuple[int, int]]
        Pair of floats if a '.' was present in string, otherwise a pair of ints
        containing the value and its esd
    """
    int_search = re.search(r'([\-\d]*)\((\d*)\)', string)
    search = re.search(r'(\-{0,1})([\d]*)\.(\d*)\((\d*)\)', string)
    if search is not None:
        # we have found a float
        sign, before_dot, after_dot, err = search.groups()
        if sign == '-':
            return -1 * (int(before_dot) + int(after_dot) * 10**(-len(after_dot))), int(err) * 10**(-len(after_dot))
        else:
            return int(before_dot) + int(after_dot) * 10**(-len(after_dot)), int(err) * 10**(-len(after_dot))
    elif int_search is not None:
        # we have found an int
        value, error = int_search.groups()
        return int(value), int(error)
    else:
        # no error found
        return float(string), 0.0


def cell_constants_to_M(a, b, c, alpha, beta, gamma):
    """
    Convert cell constants to a matrix with lattice vectors as lines.

    Parameters
    ----------
    a : float
        Cell length along a-axis.
    b : float
        Cell length along b-axis.
    c : float
        Cell length along c-axis.
    alpha : float
        Cell angle alpha in degrees.
    beta : float
        Cell angle beta in degrees.
    gamma : float
        Cell angle gamma in degrees.

    Returns
    -------
    np.ndarray
        3x3 matrix with lattice vectors.
    """
    alpha = alpha / 180.0 * np.pi
    beta = beta / 180.0 * np.pi
    gamma = gamma / 180.0 * np.pi
    M = np.array(
        [
            [
                a,
                0,
                0
            ],
            [
                b * np.cos(gamma),
                b * np.sin(gamma),
                0
            ],
            [
                c * np.cos(beta),
                c * (np.cos(alpha) - np.cos(gamma) * np.cos(beta)) / np.sin(gamma),
                c / np.sin(gamma) * np.sqrt(1.0 - np.cos(alpha)**2 - np.cos(beta)**2
                                            - np.cos(gamma)**2
                                            + 2 * np.cos(alpha) * np.cos(beta) * np.cos(gamma))
            ]
        ]
    )
    return M.T

def symm_to_matrix_vector(instruction: str) -> Tuple[np.ndarray, np.ndarray]:
    """
    Convert a symmetry instruction to a symmetry matrix and translation vector.

    Parameters
    ----------
    instruction : str
        Instruction such as '-x, -y, 0.5+z'.

    Returns
    -------
    Tuple[np.ndarray, np.ndarray]
        Tuple containing symmetry matrix and translation vector.
    """
    instruction_strings = [val.replace(' ', '').upper() for val in instruction.split(',')]
    matrix = np.zeros((3,3), dtype=np.float64)
    vector = np.zeros(3, dtype=np.float64)
    for xyz, element in enumerate(instruction_strings):
        # search for fraction in a/b notation
        fraction1 = re.search(r'(-{0,1}\d{1,3})/(\d{1,3})(?![XYZ])', element)
        # search for fraction in 0.0 notation
        fraction2 = re.search(r'(-{0,1}\d{0,1}\.\d{1,4})(?![XYZ])', element)
        # search for whole numbers
        fraction3 = re.search(r'(-{0,1}\d)(?![XYZ])', element)
        if fraction1:
            vector[xyz] = float(fraction1.group(1)) / float(fraction1.group(2))
        elif fraction2:
            vector[xyz] = float(fraction2.group(1))
        elif fraction3:
            vector[xyz] = float(fraction3.group(1))

        symm = re.findall(r'-{0,1}[\d\.]{0,8}[XYZ]', element)
        for xyz_match in symm:
            if len(xyz_match) == 1:
                sign = 1
            elif xyz_match[0] == '-':
                sign = -1
            else:
                sign = float(xyz_match[:-1])
            if xyz_match[-1] == 'X':
                matrix[xyz][0] = sign
            if xyz_match[-1] == 'Y':
                matrix[xyz][1] = sign
            if xyz_match[-1] == 'Z':
                matrix[xyz][2] = sign
    return matrix, vector

def symm_mat_vec2str(symm_mat: np.ndarray, symm_vec: np.ndarray) -> str:
    """
    Convert symmetry matrix and vector to a string representation.

    Parameters
    ----------
    symm_mat : np.ndarray
        Symmetry matrix.
    symm_vec : np.ndarray
        Symmetry vector.

    Returns
    -------
    str
        String representation of the symmetry matrix and vector.
    """
    symm_string = ''
    for symm_parts, add in zip(symm_mat, symm_vec):
        symm_string_add = str(fractions.Fraction(add).limit_denominator(50))
        if symm_string_add != '0':
            symm_string += symm_string_add
        for symm_part, symbol in zip(symm_parts, ['X', 'Y', 'Z']):
            if abs(symm_part) < 1e-10:
                continue
            if abs(1 - abs(symm_part)) < 1e-10:
                if symm_part > 0:
                    symm_string += f'+{symbol}'
                else:
                    symm_string += f'-{symbol}'
            else:
                fraction = fractions.Fraction(symm_part).limit_denominator(50)
                if str(fraction).startswith('-'):
                    symm_string += f'{str(fraction)}*{symbol}'
                else:
                    symm_string += f'+{str(fraction)}*{symbol}'
        symm_string += ','
    return symm_string[:-1]

def mean_plane2(points: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Calculate mean plane of a set of points.

    Parameters
    ----------
    points : np.ndarray
        Set of points in cartesian coordinates.

    Returns
    -------
    Tuple[np.ndarray, np.ndarray]
        Tuple containing normal vector and centroid of the plane.
    """

    points = points.T
    center = points.mean(axis=0)
    centered = points - center
    A = np.concatenate((centered[:,:2], np.ones((points.shape[0], 1))), axis=1)
    b = centered[:,2, np.newaxis]
    vector = np.linalg.inv(np.dot(A.T, A)) @ A.T @ b
    return vector[:,0] / np.linalg.norm(vector), center[:, None]

def mean_plane(points: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Calculate mean plane of a set of points.

    Parameters
    ----------
    points : np.ndarray
        Set of points in cartesian coordinates.

    Returns
    -------
    Tuple[np.ndarray, np.ndarray]
        Tuple containing normal vector and centroid of the plane.
    """
    center = points.mean(axis=0)
    centered = points - center
    try:
        _, eigenvectors = np.linalg.eigh(np.einsum('ab, ac -> bc', centered, centered))
    except:
        return mean_plane2(points)
    return eigenvectors[:, 0], center

def atan2(y: float, x: float) -> float:
    """
    Calculate atan2 of y and x.

    Parameters
    ----------
    y : float
        Y-coordinate.
    x : float
        X-coordinate.

    Returns
    -------
    float
        atan2 of y and x.
    """
    if x > 0:
        return np.arctan(np.array(y / x))
    elif x < 0 and y >= 0:
        return np.arctan(np.array(y / x)) + np.pi
    elif x < 0 and y < 0:
        return np.arctan(y / x) - np.pi
    elif x == 0 and y > 0:
        return np.pi / 2
    elif x == 0 and y < 0:
        return - np.pi / 2
    elif x == 0 and y == 0:
        return np.nan

def angle_between(vector1: np.ndarray, vector2: np.ndarray) -> float:
    """
    Calculate angle between two vectors.

    Parameters
    ----------
    vector1 : np.ndarray
        First vector.
    vector2 : np.ndarray
        Second vector.

    Returns
    -------
    float
        Angle between the two vectors in radians.
    """
    vector1 = vector1 / np.linalg.norm(vector1)
    vector2 = vector2 / np.linalg.norm(vector2)
    cosang = np.dot(vector1, vector2)
    sinang = np.linalg.norm(np.cross(vector1, vector2))
    return atan2(sinang, cosang)

def to_common_symm_str(symm_str: str) -> str:
    """
    Converts a symmetry string to a common symmetry matrix and vector string
    representation. Used to ensure that the user has flexibility for entering
    the string but comparisons still work

    Parameters
    ----------
    symm_str : str
        Symmetry string to convert.

    Returns
    -------
    str
        Common symmetry matrix and vector string representation.
    """
    symm_mat, symm_vec = symm_to_matrix_vector(symm_str)
    return symm_mat_vec2str(symm_mat, symm_vec)


def update_internal_label(internal_label: str, new_symm: str) -> str:
    """
    Updates the internal label of an atom with a new symmetry operation.

    Parameters
    ----------
    internal_label : str
        Internal label of an atom.
    new_symm : str
        New symmetry operation to apply.

    Returns
    -------
    str
        Updated internal label.
    """
    new_symm_mat, new_symm_vec = symm_to_matrix_vector(new_symm)
    old_symm, atom_label = internal_label.split(':')
    old_symm_mat, old_symm_vec = symm_to_matrix_vector(old_symm)
    transform_mat = new_symm_mat @ old_symm_mat
    transform_vec = new_symm_mat @ old_symm_vec + new_symm_vec
    return symm_mat_vec2str(transform_mat, transform_vec) + ':' + atom_label


def resolve_atom_components(
    components: List[str],
    cart_positions: np.ndarray,
    internal_labels: List[str]
) -> np.ndarray:
    """
    Resolves the atom instruction(s) within the direction definition

    Parameters
    ----------
    components : List[str]
        List of atom component instructions.
    cart_positions : np.ndarray
        Cartesian positions of atoms.
    internal_labels : List[str]
        Internal labels of atoms.

    Returns
    -------
    np.ndarray
        Resolved atom components.
    """
    if len(components) == 0:
        return np.empty((0,3))
    blank_instructions = iter(atom_component[5:].split(':') for atom_component in components)
    atom_names = [f'+X,+Y,+Z:{instr[0]}' if len(instr) == 1 else f'{to_common_symm_str(instr[0])}:{instr[1]}' for instr in blank_instructions]
    atom_indexes = np.array([internal_labels.index(name) for name in atom_names])
    return cart_positions[atom_indexes]

def resolve_centroid_components(components: List[str], centroids: np.ndarray) -> np.ndarray:
    """
    Resolves centroid instruction(s) within the direction definition

    Parameters
    ----------
    components : List[str]
        List of centroid component instructions.
    centroids : np.ndarray
        Centroids of structures.

    Returns
    -------
    np.ndarray
        Resolved centroid components.
    """
    if len(components) == 0:
        return np.empty((0,3))
    indexes = (int(comp[9:]) -1 for comp in components)
    return np.array([centroids[index] for index in indexes])

def resolve_position_components(components: List[str]) -> np.ndarray:
    """
    Resolves position instruction(s) within the direction definition.

    Parameters
    ----------
    components : List[str]
        List of position component instructions.

    Returns
    -------
    np.ndarray
        Resolved position components.
    """
    if len(components) == 0:
        return np.empty((0,3))
    return np.array([
        [float(val) for val in comp[7:].split()] for comp in components
    ])

def resolve_component_strings(
    from_string: str,
    cart_positions: np.ndarray,
    internal_labels: List[str],
    centroids: np.ndarray
) -> np.ndarray:
    """
    Resolves component strings to Cartesian coordinates for the direction
    definition.

    Parameters
    ----------
    from_string : str
        String with instructions for the components.
    cart_positions : np.ndarray
        Cartesian positions of atoms.
    internal_labels : List[str]
        Internal labels of atoms.
    centroids : np.ndarray
        Centroids of structures.

    Returns
    -------
    np.ndarray
        Resolved component strings.
    """
    components = [comp.strip() for comp in from_string.strip().split(')')]
    assert components[-1].strip() == '', 'Missing closing parentheses at the end of direction definition?'

    atom_components = [comp for comp in components if comp.startswith('atom(')]
    remaining = [comp for comp in components if not comp.startswith('atom(')]

    centroid_components = [comp for comp in remaining if comp.startswith('centroid(')]
    remaining = [comp for comp in remaining if not comp.startswith('centroid(')]

    position_components = [comp for comp in remaining if comp.startswith('position(')]
    remaining = [comp for comp in remaining if not comp.startswith('position(')]

    assert len(remaining) == 1, 'Found items that cannot be interpreted in remaining: ' + ')'.join(remaining)


    combined = np.concatenate((
        resolve_atom_components(atom_components, cart_positions, internal_labels),
        resolve_centroid_components(centroid_components, centroids),
        resolve_position_components(position_components)
    ))

    return np.mean(combined, axis=0)

def calc_overlaps(
    coords_cent_rot1: np.ndarray,
    coords_cent_rot2: np.ndarray
) -> Dict[str, float]:
    """
    Calculates overlaps between two sets of coordinates.

    Parameters
    ----------
    coords_cent_rot1 : np.ndarray
        First set of coordinates.
    coords_cent_rot2 : np.ndarray
        Second set of coordinates.

    Returns
    -------
    Dict[str, float]
        Dictionary of overlap values.
    """
    inplane1 = coords_cent_rot1[:, :2].copy()
    inplane2 = coords_cent_rot2[:, :2].copy()
    poly1 = shapely.Polygon(inplane1)
    poly2 = shapely.Polygon(inplane2)
    intersection = poly1.intersection(poly2)

    complete_area = complete_area = poly1.area + poly2.area - intersection.area

    return {
        'overlap / molecule1': intersection.area / poly1.area,
        'overlap / molecule2': intersection.area / poly2.area,
        'overlap / both': intersection.area / complete_area,
        'centroid: delta xyz': np.mean(coords_cent_rot2, axis=0) - np.mean(coords_cent_rot1, axis=0)
    }

def calc_rmat_overlap(cart_positions: np.ndarray) -> np.ndarray:
    """
    Calculates the rotation matrix to put the atoms as close into the xy plane
    as closely as possible. This is the geometry used for calculating the
    overlap.

    Parameters
    ----------
    cart_positions : np.ndarray
        Cartesian positions of atoms.

    Returns
    -------
    np.ndarray
        Rotation matrix for overlap.
    """
    axis1, centre1 = mean_plane(cart_positions)
    #coords_centred1 = cart_positions1 - centre1[None,:]

    #rotate the coordinate system so that the normal of the plane points to z
    z_unit = np.array([0, 0, 1])
    rot = np.cross(axis1, z_unit)
    rot /= np.linalg.norm(rot)
    cos_phi = np.dot(axis1, z_unit)
    n_cross = np.array(
        [[0.0, -rot[2], rot[1]],
         [rot[2], 0.0, -rot[0]],
         [-rot[1], rot[0], 0.0]]
    )
    rmat_overlap = (cos_phi * np.eye(3)
                    + np.sqrt(1 - cos_phi**2) * n_cross
                    + (1 - cos_phi) * np.outer(rot, rot))

    return rmat_overlap

def search_for_symmetry(
    symmetry_table: dict,
    cart_positions1: np.ndarray,
    cart_positions2: np.ndarray,
    M: np.ndarray,
    reference: str,
    search_values: dict
) -> Tuple[np.ndarray, str]:
    """
    Searches for symmetry operations (including translational symmetry), which
    maximises the overlap between the two sets of positions.

    Parameters
    ----------
    symmetry_table : dict
        Table of symmetry operations in cif symop format.
    cart_positions1 : np.ndarray
        First set of Cartesian positions.
    cart_positions2 : np.ndarray
        Second set of Cartesian positions.
    M : np.ndarray
        Transformation matrix.
    reference : str
        Reference for overlap calculations.
    search_values : dict
        Dictionary of search values.

    Returns
    -------
    Tuple[np.ndarray, str]
        Tuple of the best matching positions and the corresponding symmetry operation.
    """
    rmat_overlap = calc_rmat_overlap(cart_positions1)

    n_cells = search_values['n_neighbour_cells']
    x_offsets = np.arange(-n_cells, n_cells + 1, 1)
    y_offsets = np.arange(-n_cells, n_cells + 1, 1)
    z_offsets = np.arange(-n_cells, n_cells + 1, 1)
    cell_symms = [symm_to_matrix_vector(symm) for symm in symmetry_table['space_group_symop_operation_xyz']]
    axis1, centre1 = mean_plane(cart_positions1)

    overlap_max = -99999.9

    for x_offset, y_offset, z_offset, (symm_matrix, trans_vector) in product(x_offsets, y_offsets, z_offsets, cell_symms):
        offset_cart = M @ (np.array([x_offset, y_offset, z_offset]) + trans_vector)
        cart_positions2_symm = np.einsum('xy, zy -> zx', symm_matrix, cart_positions2) + offset_cart[None,:]

        min_distance = np.min(np.linalg.norm(cart_positions1[None, :] - cart_positions2_symm[:, None], axis=-1))
        if min_distance < search_values['min_distance'] or min_distance > search_values['max_distance']:
            continue

        axis2, centre2 = mean_plane(cart_positions2_symm)
        plane_plane_angle = np.pi / 2 - np.abs(np.pi / 2 - angle_between(axis1, axis2))
        if plane_plane_angle > np.deg2rad(search_values['max_angle_diff']):
            continue

        rotated = np.einsum('xy, zy -> zx', rmat_overlap, cart_positions1)
        rotated2 = np.einsum('xy, zy -> zx', rmat_overlap, cart_positions2_symm)

        this_overlap = calc_overlaps(rotated, rotated2)['overlap / ' + reference.lower()]

        if this_overlap > overlap_max and this_overlap > search_values['min_overlap_ratio']:
            best_positions2 = cart_positions2_symm
            overlap_max = this_overlap
            best_symm = (symm_matrix, np.array([x_offset, y_offset, z_offset]) + trans_vector)
    if overlap_max == -99999.9:
        raise Exception('No overlap value found with the given search_values, make sure there is an actual overlapping molecule within the range you specified')
    return best_positions2, symm_mat_vec2str(*best_symm)


def draw_polygons(coords_cent_rot1: np.ndarray, coords_cent_rot2: np.ndarray, reference: str, draw: Dict[str, str]) -> None:
    """
    Draws polygons of the overlapping molecules.

    Parameters
    ----------
    coords_cent_rot1 : np.ndarray
        First set of coordinates.
    coords_cent_rot2 : np.ndarray
        Second set of coordinates.
    reference : str
        Reference for overlap calculations.
    draw : Dict[str, str]
        Dictionary with draw instructions.
    """
    import matplotlib.pyplot as plt
    inplane1 = coords_cent_rot1[:, :2].copy()
    inplane2 = coords_cent_rot2[:, :2].copy()
    poly1 = shapely.Polygon(inplane1)
    poly2 = shapely.Polygon(inplane2)
    intersection = poly1.intersection(poly2)

    fig, ax = plt.subplots(1)
    ax.set_aspect('equal')
    concat1 = np.concatenate((inplane1[None, -1,:], inplane1))
    concat2 = np.concatenate((inplane2[None, -1,:], inplane2))

    if reference.lower() == 'molecule1':
        ax.fill(concat1[:,0], concat1[:,1],
                facecolor=draw['reference_colour'], edgecolor='#00000000')
    elif reference.lower() == 'molecule2':
        ax.fill(concat2[:,0], concat2[:,1],
                facecolor=draw['reference_colour'], edgecolor='#00000000')
    elif reference.lower() == 'both':
        ax.fill(concat1[:,0], concat1[:,1],
                facecolor=draw['reference_colour'], edgecolor='#00000000')
        ax.fill(concat2[:,0], concat2[:,1],
                facecolor=draw['reference_colour'], edgecolor='#00000000')
    else:
        raise NotImplementedError('Reference for draw not implemented. Use Molecule1, Mulecule2 or Both')

    inter_x, inter_y = intersection.exterior.xy
    ax.fill(inter_x, inter_y, facecolor=draw['overlap_colour'], edgecolor='#00000000')

    for start, end in zip(concat1[:-1], concat1[1:]):
        ax.plot((start[0], end[0]), (start[1], end[1]), c=draw['edge_colour'])

    for start, end in zip(concat2[:-1], concat2[1:]):
        ax.plot((start[0], end[0]), (start[1], end[1]), c=draw['edge_colour'])

    ax.axis('off');
    log.info('')
    if draw['output_file'] != 'web.display':
        fig.savefig(draw['output_file'], dpi=draw['output_dpi'], bbox_inches='tight')
        log.info(f"saved graphics to file: {draw['output_file']}")
    return fig

def options_from_config(
    filename: str
) -> Tuple[str, str, str, str, str, str, Dict[str, float], Dict[str, str]]:
    """
    Reads options from a config file.

    Parameters
    ----------
    filename : str
        Path to the config file.

    Returns
    -------
    Tuple[str, str, str, str, str, str, Dict[str, float], Dict[str, str]]
        Tuple with cif_path, dataset_name, atoms_string1, atoms_string2,
        direction_string, reference, search_values, draw.
    """
    parser = ConfigParser()
    parser.read(filenames=[filename])
    cif_path = parser['cif']['path']
    dataset_name = parser['cif']['dataset_label']
    atoms_string1 = parser['cif']['molecule1']
    atoms_string2 = parser['cif']['molecule2']

    direction_string = parser['geometry']['direction2']
    reference = parser['geometry']['reference']

    search_values = {
        'search_molecule_2': parser['search'].getboolean('molecule2_search'),
        'min_distance': parser['search'].getfloat('min_distance'),
        'max_distance': parser['search'].getfloat('max_distance'),
        'max_angle_diff': parser['search'].getfloat('max_angle_diff'),
        'min_overlap_ratio': parser['search'].getfloat('min_overlap_ratio'),
        'n_neighbour_cells': parser['search'].getint('n_neighbour_cells')
    }

    draw = {
        'overlap_colour': parser['figure']['overlap_colour'],
        'reference_colour': parser['figure']['reference_colour'],
        'edge_colour': parser['figure']['edge_colour'],
        'output_file': parser['figure']['output_file'],
        'output_dpi': parser['figure'].getfloat('output_dpi')
    }
    return cif_path, dataset_name, atoms_string1, atoms_string2, direction_string, reference, search_values, draw



def cif2values(cif_fo, dataset_label, atoms_string1, atoms_string2, reference, search_values):
    cif_dict = ciflike_to_dict(cif_fo)
    try:
        cif = cif_dict[dataset_label]
    except KeyError:
        key = list(cif_dict.keys())[int(dataset_label)]
        cif = cif_dict[key]
    cell = np.array([cif['cell_length_a'],
                     cif['cell_length_b'],
                     cif['cell_length_c'],
                     cif['cell_angle_alpha'],
                     cif['cell_angle_beta'],
                     cif['cell_angle_gamma']])
    M = cell_constants_to_M(*cell)

    atom_site_table = next(loop for loop in cif['loops'] if 'atom_site_fract_x' in loop.columns)
    atom_site_table = atom_site_table.set_index('atom_site_label')
    coord_columns = tuple('atom_site_fract_' + coord for coord in ('x', 'y', 'z'))
    #all_atom_labels = list(atom_site_table['atom_site_label'])

    symms1 = [group.split(':') for group in atoms_string1.strip().split(';')]

    symmetry_table = next(table for table in cif['loops'] if 'space_group_symop_operation_xyz' in table.columns or 'symmetry_equiv_pos_as_xyz' in table.columns).copy()
    symmetry_table = symmetry_table.rename({'symmetry_equiv_pos_as_xyz': 'space_group_symop_operation_xyz'}, axis=1)

    cart_positions1 = np.zeros((0, 3))
    internal_labels1 = []

    for symm, atom_names in symms1:
        symm_mat, symm_vec = symm_to_matrix_vector(symm)
        common_symm_str = symm_mat_vec2str(symm_mat, symm_vec)
        atom_labels = atom_names.strip().split(' ')
        atom_pos = atom_site_table.loc[atom_labels, coord_columns].values
        symm_pos = np.einsum('xy, zy -> zx', symm_mat, atom_pos) + symm_vec[None, :]
        cart_pos = np.einsum('xy, zy -> zx', M, symm_pos)
        cart_positions1 = np.concatenate((cart_positions1, cart_pos))
        internal_labels1 += [common_symm_str  + ':' + label for label in atom_labels]

    symms2 = [group.split(':') for group in atoms_string2.strip().split(';')]
    cart_positions2 = np.zeros((0, 3))
    internal_labels2 = []

    for symm, atom_names in symms2:
        symm_mat, symm_vec = symm_to_matrix_vector(symm)
        common_symm_str = symm_mat_vec2str(symm_mat, symm_vec)
        atom_labels = atom_names.strip().split(' ')
        atom_pos = atom_site_table.loc[atom_labels, coord_columns].values
        symm_pos = np.einsum('xy, zy -> zx', symm_mat, atom_pos) + symm_vec[None, :]
        cart_pos = np.einsum('xy, zy -> zx', M, symm_pos)
        cart_positions2 = np.concatenate((cart_positions2, cart_pos))
        internal_labels2 += [common_symm_str  + ':' + label for label in atom_labels]

    if search_values['search_molecule_2']:
        cart_positions2, new_symm = search_for_symmetry(symmetry_table, cart_positions1, cart_positions2, M, reference, search_values)
        log.info('Found symmetry transformation for molecule 2 with the highest overlap')
        log.info(f'Operation: {new_symm}')
        internal_labels2 = [update_internal_label(internal_label, new_symm) for internal_label in internal_labels2]

    return cart_positions1, cart_positions2, internal_labels1, internal_labels2


def overlap_calculation(cif_fo, dataset_name: str, atoms_string1: str, atoms_string2: str, direction_string: str, reference: str, search_values: Dict[str, Union[bool, float]], draw: Dict[str, Union[str, int]]) -> None:
    """
    Calculates the overlap between two planar moieties in a crystal structure.

    Parameters
    ----------
    cif_of : FileObject
        FileObject to the Crystallographic Information File (CIF).
    dataset_name : str
        Label of the dataset within the CIF file. This can be an integer or string.
        Integer indices (0 for the first dataset, 1 for the second, etc.) can be used, but
        string interpretation will be tried first.
    atoms_string1 : str
        String containing atom labels and symmetry operations for the first reference molecule.
        This molecule determines the reference plane.
    atoms_string2 : str
        String containing atom labels and symmetry operations for the second molecule. The
        script can search for potential overlaps with this molecule.
    direction_string : str
        String specifying the direction of the overlap, composed of a directional vector in the
        XY plane, followed by two sets of positions (centroid or atom) separated by a comma.
        All values within a set will be averaged.
    reference : str
        Specifies the reference molecule for calculating overlaps. Options include 'molecule1',
        'molecule2' (for the respective molecules' areas), or 'both' (for the union of both
        molecules' projected areas).
    search_values : Dict[str, Union[bool, float]]
        Dictionary of parameters used in the search for the largest overlap between molecules.
        Includes whether to conduct a search ('molecule2_search'), minimum interatomic distance
        ('min_distance'), maximum interatomic distance ('max_distance'), maximum angle difference
        ('max_angle_diff'), minimum overlap ratio ('min_overlap_ratio'), and number of neighboring
        cells to search for translational symmetry ('n_neighbour_cells').
    draw : Dict[str, Union[str, int]]
        Dictionary with drawing instructions. Includes colors for overlap, remaining area, and
        polygon edges ('overlap_colour', 'reference_colour', 'edge_colour'), filename for output
        ('output_file'), and output dpi ('output_dpi').

    Returns
    -------
    None
    """
    cart_positions1, cart_positions2, internal_labels1, internal_labels2= cif2values(cif_fo, dataset_name, atoms_string1, atoms_string2, reference, search_values)
    log.info(f'Read atoms from dataset {dataset_name}')

    centroids = (np.mean(cart_positions1, axis=0),
                 np.mean(cart_positions2, axis=0))

    to_string, from_string = direction_string.split(':')

    to_vector_inp = np.array([float(val.strip()) for val in to_string.strip().split(' ')])
    to_vector_inp /= np.linalg.norm(to_vector_inp)
    to_vector = np.zeros(3)
    to_vector[:2] = to_vector_inp

    from_string1, from_string2 = (val.strip() for val in from_string.split(','))

    internal_labels = internal_labels1 + internal_labels2
    cart_positions_con = np.concatenate((cart_positions1, cart_positions2), axis=0)

    from_vector1 = resolve_component_strings(from_string1, cart_positions_con, internal_labels, centroids)

    from_vector2 = resolve_component_strings(from_string2, cart_positions_con, internal_labels, centroids)

    from_vector = from_vector2 - from_vector1
    from_vector /= np.linalg.norm(from_vector)

    axis1, centre1 = mean_plane(cart_positions1)
    axis2, centre2 = mean_plane(cart_positions2)

    coords_centred1 = cart_positions1 - centre1[None,:]
    coords_centred2 = cart_positions2 - centre1[None,:]

    r_mat1 = calc_rmat_overlap(cart_positions1)
    from_vec_rot = r_mat1 @ from_vector

    # vector is put into xy plane
    from_vec_rot[2] = 0

    rot = np.cross(from_vec_rot, to_vector)
    rot /= np.linalg.norm(rot)
    cos_phi = np.dot(from_vec_rot, to_vector)
    n_cross = np.array(
        [[0.0, -rot[2], rot[1]],
         [rot[2], 0.0, -rot[0]],
         [-rot[1], rot[0], 0.0]]
    )
    r_mat2 = (cos_phi * np.eye(3)
             + np.sqrt(1 - cos_phi**2) * n_cross
             + (1 - cos_phi) * np.outer(rot, rot))
    r_mat = r_mat2 @ r_mat1
    coords_cent_rot1 = np.einsum('xy, zy -> zx', r_mat, coords_centred1)
    coords_cent_rot2 = np.einsum('xy, zy -> zx', r_mat, coords_centred2)
    log.info('')
    log.info('Atoms of molecule1 in rotated and centred geometry. Z-component not yet set to 0.')
    log.info(f'{"  symm:label":18s} {"    X":8s} {"    Y":8s} {"    Z":8s}')
    for atom_name, coord in zip(internal_labels1, coords_cent_rot1):
        log.info(f'{atom_name:18s} {coord[0]: 8.5f} {coord[1]: 8.5f} {coord[2]: 8.5f}')
    centroid1 = np.mean(coords_cent_rot1, axis=0) # this should be 0 0 0
    log.info(f'{"centroid":18s} {centroid1[0]: 8.5f} {centroid1[1]: 8.5f} {centroid1[2]: 8.5f}')

    log.info('')

    log.info('Atoms of molecule2 in rotated and centred geometry. Z-component not yet set to 0.')
    log.info(f'{"  symm:label":18s} {"    X":8s} {"    Y":8s} {"    Z":8s}')
    for atom_name, coord in zip(internal_labels2, coords_cent_rot2):
        log.info(f'{atom_name:18s} {coord[0]: 8.5f} {coord[1]: 8.5f} {coord[2]: 8.5f}')
    centroid2 = np.mean(coords_cent_rot2, axis=0)
    log.info(f'{"centroid":18s} {centroid2[0]: 8.5f} {centroid2[1]: 8.5f} {centroid2[2]: 8.5f}')
    log.info('')

    plane_plane_angle = np.pi / 2 - np.abs(np.pi / 2 - angle_between(axis1, axis2))

    log.info(f'angle between normals of planes (deg): {np.rad2deg(plane_plane_angle): 5.2f}')
    log.info('')
    overlaps = calc_overlaps(coords_cent_rot1, coords_cent_rot2)

    log.info(f'Overlap is determined in reference to: {reference}')
    chosen_overlap = overlaps['overlap / ' + reference.lower()]
    log.info(f'Determined overlap in %: {100.0*chosen_overlap:6.3f}')

    fig = draw_polygons(coords_cent_rot1, coords_cent_rot2, reference, draw)
    return 100.0*chosen_overlap, fig

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog='Overlap Calculator',
        description='This program calculates the overlap between different (almost) planar moietys in a crystal structure.'
    )
    parser.add_argument('CONFIG_FILE', help='Filename of the ini file for the configureation. If if does not exist an example file will be generated at the given path')
    args = parser.parse_args()
    fh = logging.FileHandler('overlap.log', mode='w', encoding='utf-8')
    fh.setLevel(logging.INFO)
    fh.setFormatter(formatter)
    log.addHandler(fh)

    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    ch.setFormatter(formatter)
    log.addHandler(ch)

    config_file = args.CONFIG_FILE

    if not os.path.exists(config_file):
        with open(config_file, 'w') as fo:
            fo.write(config_template)
    else:
        cif_path, dataset_name, atoms_string1, atoms_string2, direction_string, reference, search_values, draw = options_from_config(config_file)

        with open(cif_path, 'r') as cif_fo:
            _ = overlap_calculation(cif_fo, dataset_name, atoms_string1, atoms_string2, direction_string, reference, search_values, draw)





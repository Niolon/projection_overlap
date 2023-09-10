import streamlit as st
from io import StringIO
from overlap import options_from_config, overlap_calculation, log, formatter
import logging
import textwrap


st.title('Overlap calculator beta')

uploaded_cif = st.file_uploader("Upload a cif file")
with st.expander("See explanation"):
    st.markdown("A valid cif file containing a data block, with the _cell_length and _cell_angle entries, as well as a loop with atom_site entries for label and fract_x,y,z.")

dataset_name = st.text_input('Used dataset within cif', '0')
with st.expander("See explanation"):
    st.markdown(textwrap.dedent("""
        label of the dataset within the cif file. An integer number can be used to
        select the first (0) or second(1) dataset. Interpretation as a string name
        will always be tried first (in case there is a dataset named 1)
    """).strip())

conditions = [uploaded_cif is not None]

atoms_string1 = st.text_input('Atoms forming Polygon1', 'x, y, z: atom_names')
with st.expander("See explanation"):
    st.markdown(textwrap.dedent("""
        First of your reference molecules, this is the molecule that determines the
        reference plane.
        Beginning is a symmetry code in cif convention then a colon and then atom names.
        For the atom names, the ordering matters! Edges of the polygon are between
        neighbouring atoms and the final and first polygon. If multiple symmetry codes
        are needed, these can be separated by a semicolon.
    """.strip()))

atoms_string2 = st.text_input('Atoms forming Polygon2', 'x, y, z: atom_names')
with st.expander("See explanation"):
    st.markdown(textwrap.dedent("""
        Base geometry of the second molecule. Can be used as given or potential
        overlaps can be searched by the script. Syntax is the same as in Atoms forming Polygon1
    """).strip())
conditions.append(atoms_string1 != 'x, y, z: atom_names')
conditions.append(atoms_string2 != 'x, y, z: atom_names')

search_mol2 = st.checkbox('Try symmetry equivalent positions for molecule 2')
with st.expander("See explanation"):
    st.markdown(textwrap.dedent("""
        If active try to find the molecule 2 with the the maximum overlap within
        the given search parameters by applying the unit cell symmetry operations,
        as well as translation to neighbouring cells.
    """).strip())

if search_mol2:
    min_distance = st.number_input(
        'Minimal distance between overlapping moieties.',
        min_value=0.0,
        max_value=20.0,
        value=0.7
    )

    max_distance = st.number_input(
        'Maximal distance between overlapping moieties.',
        min_value=1.0,
        max_value=22.0,
        value=5.0
    )

    max_angle_diff= st.number_input(
        'Maximal angle between molecular planes in degree.',
        min_value=0.0,
        max_value=90,
        value=45.0
    )

    min_overlap = st.number_input(
        'Minimal overlap to be considered in the the search in percent.',
        min_value=0.0,
        max_value=100.0,
        value=0.
    )


direction_string = st.text_input(
    'Defined reference direction',
    '1 0: centroid(1), atom(atom_name1) atom(atom_name2)'
)

with st.expander("See explanation"):
    st.markdown(textwrap.dedent("""
        While this setting does not affect the calculation of the overlap itself,
        the reference direction can be chosen to get additional information about
        offsets in the xy plane, as well as customise the rotation in the produced plot.
        Before the colon is a directional vector in the xy
        plane. After the colon the script expects two sets of positions separated by
        a comma. All values within a set will be averaged. Three options are available
        for the definition.

        position(0.1 0.2 0.3) is a position in cartesian coordinates.

        centroid(1)/centroid(2) are the centroids of molecule1 and 2 respectively.

        atom(C2) is an atom position of an atom which is available with the symmetry x, y, z in either Molecule 1 or
        Molecule 2. If a symmetry equivalent position is needed, it can only be input as position().
    """).strip())

conditions.append(direction_string != '1 0: centroid(1), atom(atom_name1) atom(atom_name2)')


if all(conditions):
    stringio = StringIO(uploaded_cif.getvalue().decode("utf-8"))

    config_file = 'overlap_config.ini'
    #_, _, _, _, _, reference, search_values, draw = options_from_config(config_file)

    reference = 'molecule1'
    search_values = {
        'search_molecule_2': search_mol2,
        'min_distance': min_distance,
        'max_distance': max_distance,
        'max_angle_diff': max_angle_diff,
        'min_overlap_ratio': min_overlap,
        'n_neighbour_cells': 3
    }
    draw = {
        'overlap_colour': '#001242',
        'reference_colour': '#7EA8BE',
        'edge_colour': '#a0a0a0',
        'output_file': 'web.display',
        'output_dpi': 300.0
    }

    log_stream = StringIO()
    ch = logging.StreamHandler(log_stream)
    ch.setLevel(logging.INFO)
    formatter = logging.Formatter('%(message)s')
    ch.setFormatter(formatter)
    log.addHandler(ch)

    draw['output_file'] = 'web.display'
    try:
        overlap, fig = overlap_calculation(stringio, dataset_name, atoms_string1, atoms_string2, direction_string, reference, search_values, draw)
    except Exception as e:
        st.write('Could not process your input. Please recheck and consult the explanations')
        raise e

    st.pyplot(fig)

    st.code((log_stream.getvalue()), language=None)


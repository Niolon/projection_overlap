import streamlit as st
from io import StringIO, BytesIO
from overlap import options_from_config, overlap_calculation, log, formatter
import logging
import base64
import textwrap
import zipfile

def render_svg(svg):
    """Renders the given svg string."""
    b64 = base64.b64encode(svg.encode('utf-8')).decode("utf-8")
    html = r'<img src="data:image/svg+xml;base64,%s"/>' % b64
    st.write(html, unsafe_allow_html=True)

st.title('Overlap calculator beta')

st.markdown(textwrap.dedent("""
    *by Paul Niklas Ruth (paul.n.ruth@durham.ac.uk)*

    This script calculates the overlap between two moieties as the quotient of
    projected areas. This means the mean plane through a first moiety is determined.
    Subsequently the atomic positions of both moieties are projected onto that plane
    and two polygons are constructed from the two investigated moieties. The overlap
    can now be calculated as the quotient of the intersection area and the area of the
    first moiety polygon.

    The algorithm is a reimplementation of the one described in my PhD thesis (see references).

    ### References:

    If you use this website in you research cite the following two references.

    - P. N. Ruth, Method Development for Benchmarking Key Interactions in Quantum Crystallography, p. 39, http://dx.doi.org/10.53846/goediss-9798

    - P. N. Ruth, Overlap calculator beta 0.1.1, DOI: 10.5281/zenodo.8332656

    ### Program

    The program will runand output the result once all necessary information has been filled in and will adapt interactively (Thanks to Streamlit).

""").strip())

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

search_mol2 = st.checkbox('Search best symmetry equivalent position for molecule 2')
with st.expander("See explanation"):
    st.markdown(textwrap.dedent("""
        If active try to find the molecule 2 with the the maximum overlap within
        the given search parameters by applying the unit cell symmetry operations,
        as well as translation to neighbouring cells. The input parameters still
        need to form a closed polygon with the given symmetry elements and
        atom names.
    """).strip())

if search_mol2:
    min_distance = st.number_input(
        'Minimal distance between two atoms in overlapping moieties.',
        min_value=0.0,
        max_value=20.0,
        value=0.7,
        step=0.1
    )
    with st.expander("See explanation"):
        st.markdown(textwrap.dedent("""
            If any two atoms between the two input moieties are closer than this value
            in Angstrom, this overlap will not be considered. Prevents
            self-overlap and to a degree enables investigation of disordered
            moieties.
        """).strip())

    max_distance = st.number_input(
        'Maximal distance between two atoms in overlapping moieties.',
        min_value=1.0,
        max_value=22.0,
        value=5.0,
        step=0.1
    )

    with st.expander("See explanation"):
        st.markdown(textwrap.dedent("""
            If no two atoms between the two input moieties are closer than this value
            in Angstrom, this overlap will not be considered. A smaller value will
            speed up the investigation and prevent overlap calculations between non-adjacent
            molecules.
        """).strip())

    max_angle_diff= st.number_input(
        'Maximal angle between molecular planes in degree.',
        min_value=0.0,
        max_value=90.0,
        value=45.0,
        step=1.0
    )

    with st.expander("See explanation"):
        st.markdown(textwrap.dedent("""
            If the angle between the normal vectors of the mean planes of the two moieties
            is larger than this value in degree, the overlap will not be considered.
            This prevents the calculation of overlaps in T-shaped arrangements.
        """).strip())

    min_overlap = st.number_input(
        'Minimal overlap to be considered in the the search in percent.',
        min_value=0.0,
        max_value=100.0,
        value=1.,
        step=1.0
    )
    with st.expander("See explanation"):
        st.markdown(textwrap.dedent("""
            If the calculated overlap between moieties is smaller than this value in percent
            the overlap is not considered.
        """).strip())

    n_neighbour_cells = st.number_input(
        'Number of translations to check into neighbouring cells',
        min_value=0,
        max_value=5,
        value=1
    )

    with st.expander("See explanation"):
        st.markdown(textwrap.dedent("""
            Number of neighbouring cells to be considered in searching for the symmetry of
            the second moiety. This value might be needed if a lot of your fractional
            coordinates are located outside of the base unit cell (so the coordinates
            are 1.523 or something similar). In that case you might want to recentre your
            molecule in the unit cell. Large values will slow down the script and might
            incur weight times.
        """).strip())
else:
    min_distance = None
    max_distance = None
    max_angle_diff = None
    min_overlap = 0.0
    n_neighbour_cells = None


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

        ### Example
        For a given naphtalene molecule with the following atom names.
    """).strip())

    with open('images/reference_explanation.svg', 'r') as fobj:
        render_svg(fobj.read())

    st.markdown(textwrap.dedent("""
        We would like the reference vector to point from the centroid of this molecule to the centre
        centre of one of the outermost bonds to orient this direction horizontally.

        The input reference direction would be: "1 0: centroid(1), atom(C2) atom(C3)"

        or alternatively "1 0: centroid(1), atom(C7) atom(C8)"

        Note:

        The determination becomes more robust the longer the distance from the centroid actually is.
        A less robust alternative therefore would be: "0 1: centroid(1), atom(C10)".
    """).strip())

conditions.append(direction_string != '1 0: centroid(1), atom(atom_name1) atom(atom_name2)')

with st.expander("Plotting options"):
    col1, col2, col3 = st.columns(3)
    with col1:
        overlap_colour = st.color_picker('Colour of overlap area', value='#001242')
    with col2:
        reference_colour = st.color_picker('Colour of reference moiety area', value='#7EA8BE')
    with col3:
        edge_colour = st.color_picker('Colour of polygon edges', value='#a0a0a0')
    edge_width = st.number_input(
        'Line width of the drawn polygon edges',
        min_value=0.0,
        max_value=50.0,
        value=3.0,
        step=0.1
    )


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
        'min_overlap_ratio': min_overlap / 100,
        'n_neighbour_cells': n_neighbour_cells
    }
    draw = {
        'overlap_colour': overlap_colour,
        'reference_colour': reference_colour,
        'edge_colour': edge_colour,
        'line_width': edge_width,
        'output_file': 'web.display',
        'output_dpi': 300.0
    }

    log_stream = StringIO()
    ch = logging.StreamHandler(log_stream)
    ch.setLevel(logging.INFO)
    formatter = logging.Formatter('%(message)s')
    ch.setFormatter(formatter)
    log.addHandler(ch)

    try:
        overlap, fig = overlap_calculation(stringio, dataset_name, atoms_string1, atoms_string2, direction_string, reference, search_values, draw)
    except Exception as e:
        st.write('Could not process your input. Please recheck and consult the explanations')
        raise e

    st.pyplot(fig)

    txt_output = log_stream.getvalue()
    st.code((txt_output), language=None)

    zip_bool = st.checkbox(
        'Create Output zip (takes significantly longer to process/reload)'
    )
    if zip_bool:
        references_output = textwrap.dedent("""
            P. N. Ruth, Method Development for Benchmarking Key Interactions in Quantum Crystallography, p. 39, http://dx.doi.org/10.53846/goediss-9798

            P. N. Ruth, Overlap calculator beta 0.1.1, DOI: 10.5281/zenodo.8332656
        """).strip()

        fig_buffer = BytesIO()
        fig.savefig(fig_buffer, format='png')
        #fig_buffer.seek(0)

        zip_buffer = BytesIO()
        with zipfile.ZipFile(zip_buffer, 'a', zipfile.ZIP_DEFLATED, False) as zip_file:
            zip_file.writestr('overlap_log.txt', txt_output)
            zip_file.writestr('references.txt', references_output)
            zip_file.writestr('overlap.png', fig_buffer.getvalue())
        st.download_button(
            'Download as zip file',
            zip_buffer.getvalue(),
            'overlap.zip',
            mime='application/zip'
        )




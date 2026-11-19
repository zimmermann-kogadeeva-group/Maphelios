#!/usr/bin/env python3
import argparse
import base64
from pathlib import Path

import streamlit as st
from page1 import get_main_inputs, run_pipeline
from page2_longread import longread_results
from page2_sanger import sanger_results
from page2_shortread import shortread_results
from PIL import Image

from maphelios import __version__

st.set_page_config(page_title="Maphelios", layout="wide")

if "stage" not in st.session_state:
    st.session_state.stage = 0

if "example_active" not in st.session_state:
    st.session_state.example_active = False

LOGO_PATH = Path(__file__).parent / "maphelios_logo.png"
SIDEBAR_LOGO_WIDTH = 180


@st.cache_data
def logo_html(width=None):
    """
    HTML for the logo that follows the app's text colour.

    The logo is a single colour on a transparent background, so it is used
    as a CSS mask over a block filled with ``currentColor``. The browser then
    draws it in whatever text colour the active theme uses, which stays
    correct when the theme is changed (unlike ``st.context.theme``, which
    can be stale on first load and after switching theme).

    ``width`` is the display width in pixels (aspect ratio is preserved);
    defaults to the image's native width.
    """
    with Image.open(LOGO_PATH) as img:
        native_width, native_height = img.size
    width = width or native_width
    height = round(native_height * width / native_width)
    data = base64.b64encode(LOGO_PATH.read_bytes()).decode()
    mask = f"url(data:image/png;base64,{data}) no-repeat center / contain"
    return (
        f'<div role="img" aria-label="Maphelios" style="'
        f"width:{width}px; height:{height}px; max-width:100%; "
        f"background-color:currentColor; "
        f"-webkit-mask:{mask}; mask:{mask};"
        f'"></div>'
    )


def submit(*args):
    with st.spinner("Processing..."):
        all_res = dict()
        for name, param_set in args[0].items():
            res = run_pipeline(**param_set)
            if res is not None:
                all_res[name] = res

        if len(all_res):
            st.session_state.results = all_res
            st.session_state.stage = 1


def main():
    parser = argparse.ArgumentParser("Maphelios")
    parser.add_argument("--workdir", action="store_true")
    args = parser.parse_args()

    col1, col2 = st.columns((6, 1), gap="large", vertical_alignment="bottom")

    if st.session_state.stage == 0:
        st.html(logo_html())
        st.caption(f"*:grey[v{__version__}]*")
        all_inputs = get_main_inputs(args.workdir)
        st.button("Submit", on_click=submit, args=[all_inputs], width="stretch")

    if st.session_state.stage >= 1:
        # Rendered before the results pages add their sidebar widgets, so it
        # sits at the top of the sidebar
        st.sidebar.html(logo_html(width=SIDEBAR_LOGO_WIDTH))
        if st.session_state.seq_type == "blast":
            sanger_results()
        elif st.session_state.seq_type == "minimap2":
            longread_results()
        elif st.session_state.seq_type == "bowtie2":
            shortread_results()


if __name__ == "__main__":
    main()

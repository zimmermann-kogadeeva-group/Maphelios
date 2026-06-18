#!/usr/bin/env python3
import argparse

import streamlit as st
from page1 import get_main_inputs, run_pipeline
from page2 import show_results

from maphelios import __version__

st.set_page_config(layout="wide")

if "stage" not in st.session_state:
    st.session_state.stage = 0


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
    st.title(f"Maphelios")
    st.caption(f"*:grey[v{__version__}]*")

    if st.session_state.stage == 0:
        all_inputs = get_main_inputs(args.workdir)
        st.button(
            "Submit", on_click=submit, args=[all_inputs], use_container_width=True
        )

    if st.session_state.stage >= 1:
        show_results()


if __name__ == "__main__":
    main()

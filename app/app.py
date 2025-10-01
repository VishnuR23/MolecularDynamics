import json
import streamlit as st
import plotly.graph_objects as go

st.set_page_config(page_title="MolDynamicsNet", layout="wide")

st.title("MolDynamicsNet — Binding Energy Dashboard")
pred_file = st.text_input("Prediction JSON path", "runs/demo_pred.json")

if st.button("Load Prediction"):
    with open(pred_file) as f:
        pred = json.load(f)
    st.json(pred)
    curve = pred.get("curve", [])
    fig = go.Figure()
    fig.add_trace(go.Scatter(y=curve, mode='lines+markers', name='Binding curve (toy)'))
    st.plotly_chart(fig, use_container_width=True)

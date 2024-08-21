# utils.py

import base64
import streamlit as st
from pathlib import Path
import requests

def load_image(image_path: str) -> str:
    """Load and encode an image file to base64."""
    try:
        with open(image_path, "rb") as img_file:
            return base64.b64encode(img_file.read()).decode('utf-8')
    except FileNotFoundError:
        st.error(f"Image not found: {image_path}")
        return ""
    
def load_image_from_url(url):
    response = requests.get(url)
    if response.status_code == 200:
        return base64.b64encode(response.content).decode('utf-8')
    else:
        return None
    

def load_image(image_path: str) -> str:
    """Load and encode an image file to base64."""
    try:
        if image_path.startswith("http://") or image_path.startswith("https://"):
            response = requests.get(image_path)
            if response.status_code == 200:
                return base64.b64encode(response.content).decode('utf-8')
            else:
                st.error(f"Unable to load image from URL: {image_path}")
                return ""
        else:
            with open(image_path, "rb") as img_file:
                return base64.b64encode(img_file.read()).decode('utf-8')
    except FileNotFoundError:
        st.error(f"Image not found: {image_path}")
        return ""

def add_custom_css():
    """Add custom CSS to improve the UI."""
    st.markdown("""
        <style>
        body { background-color: #f5f5f5; font-family: Arial, sans-serif; }
        h1, h2, h3, h4, h5, h6 { color: #333; font-weight: bold; }
        .sidebar .sidebar-content { background-color: #37474f; color: #fff; }
        .sidebar .sidebar-content h1 { color: #3840c9; }
        .stRadio > label { font-weight: bold; color: #3840c9; }
        .stFileUploader { border: 2px dashed #9e9e9e; padding: 20px; background-color: #fafafa; }
        .stButton > button {
            background-color: #3840c9; color: white; font-size: 16px;
            padding: 10px 20px; border: none; border-radius: 5px; cursor: pointer;
        }
        .stButton > button:hover { background-color: #4c4cb0; }
        .stDataFrame { margin-top: 20px; border: 1px solid #ddd; border-radius: 5px; }
        .stAlert.success { background-color: #d4edda; color: #155724; }
        .stAlert.error { background-color: #f8d7da; color: #721c24; }
        img { max-width: 100%; height: auto; border-radius: 10px; margin-bottom: 20px; }
        </style>
    """, unsafe_allow_html=True)


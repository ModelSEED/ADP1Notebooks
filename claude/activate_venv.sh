#!/bin/bash
# Quick activation script for ADP1Research virtual environment
source "$(dirname "${BASH_SOURCE[0]}")/.venv/bin/activate"
echo "✅ Virtual environment activated"
echo "📍 Python: $(which python)"
echo "🐍 Version: $(python --version)"
echo ""
echo "Available commands:"
echo "  jupyter notebook  - Start Jupyter Notebook"
echo "  jupyter lab       - Start JupyterLab"
echo "  python            - Python interpreter"
echo "  deactivate        - Exit virtual environment"

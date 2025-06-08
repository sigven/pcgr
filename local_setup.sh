#!/bin/bash
set -e

echo "🚀 DGG Rules Somatic - Local Setup Script"
echo "=========================================="

echo "📋 Checking dependencies..."

if ! command -v git &> /dev/null; then
    echo "❌ git is required but not installed"
    exit 1
fi
echo "✅ git found"

if ! command -v python3 &> /dev/null; then
    echo "❌ python3 is required but not installed"
    exit 1
fi
echo "✅ python3 found"

if ! command -v curl &> /dev/null; then
    echo "❌ curl is required but not installed"
    exit 1
fi
echo "✅ curl found"

if ! command -v poetry &> /dev/null; then
    echo "📦 Installing Poetry..."
    curl -sSL https://install.python-poetry.org | python3 -
    export PATH="$HOME/.local/bin:$PATH"
else
    echo "✅ Poetry found"
fi

if [ ! -f "pyproject.toml" ]; then
    echo "📥 Cloning repository..."
    git clone https://github.com/LauferVA/dgg_rules_somatic.git
    cd dgg_rules_somatic
fi

echo "📦 Installing Python dependencies..."
poetry install --no-root

echo "🗄️ Setting up database..."
poetry run python -c "
import sys
import os
sys.path.insert(0, os.getcwd())
from api.database import create_tables
create_tables()
"

echo "📚 Building knowledge base..."
cd "$(pwd)" && PYTHONPATH="$(pwd)" poetry run python etl/build_reportable.py

echo "✅ Setup complete!"
echo ""
echo "🎯 Next steps:"
echo "  1. Activate environment: poetry shell"
echo "  2. Start API: uvicorn api.main:app --reload"
echo "  3. Start demo: streamlit run demo_ui.py"

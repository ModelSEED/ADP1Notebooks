# Notebooks

This directory contains Jupyter notebooks for exploratory data analysis, experiments, and documentation.

## Structure

- **exploratory/**: Notebooks for initial data exploration and analysis
- **experiments/**: Notebooks for testing hypotheses and running experiments
- **documentation/**: Notebooks that serve as documentation or tutorials

## Usage

1. Install Jupyter in your development environment:

   ```bash
   uv add --group dev jupyter notebook
   ```

2. Start the Jupyter server:

   ```bash
   jupyter notebook
   ```

3. Create new notebooks in the appropriate subdirectory based on their purpose.

## Best Practices

- Use descriptive names for notebooks (e.g., `01_data_exploration.ipynb`)
- Include a clear title and description at the top of each notebook
- Clean up outputs before committing to version control
- Consider using `vibe.md` for brainstorming and planning before creating notebooks

## Template Files

- `vibe.md`: Template for exploratory prompts and brainstorming
- Use this to document ideas, hypotheses, and experimental plans

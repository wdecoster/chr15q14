# GWAS Manhattan Dash Application

This project is a Dash application designed to generate Manhattan plots for Genome-Wide Association Studies (GWAS). The application allows users to filter data based on specified criteria and visualize the results interactively.

## Project Structure

```
gwas-manhattan-dash
├── src
│   ├── app.py                # Main entry point of the Dash application
│   ├── components
│   │   ├── __init__.py       # Empty initializer for the components package
│   │   ├── filters.py         # Implementation of filter components
│   │   └── manhattan_plot.py  # Implementation of the Manhattan plot component
│   ├── data
│   │   └── __init__.py       # Empty initializer for the data package
│   ├── utils
│   │   ├── __init__.py       # Empty initializer for the utils package
│   │   ├── data_loader.py     # Functions to load and preprocess GWAS data
│   │   └── plot_helpers.py     # Helper functions for generating the Manhattan plot
│   └── static
│       └── styles.css         # CSS styles for the Dash application
├── assets                      # Directory for additional assets (images, JS files, etc.)
├── requirements.txt            # Lists dependencies required for the project
└── README.md                   # Documentation for the project
```

## Installation

To set up the project, clone the repository and install the required dependencies:

```bash
git clone <repository-url>
cd gwas-manhattan-dash
pip install -r requirements.txt
```

## Running the Application

Note that the current implementation has the path to the data file hardcoded. You will need to modify the path in `src/utils/data_loader.py` to point to your local GWAS data file, which has some expectations regarding its format and columns, requiring customization based on your specific dataset. To run the Dash application, execute the following command:

```bash
python src/app.py
```

Once the application is running, open your web browser and navigate to `http://127.0.0.1:8050` to access the interface.

## Features

- Interactive Manhattan plot generation based on GWAS data.
- Filtering options using sliders and checkboxes for specified columns.
- User-friendly layout for selecting filters and visualizing results.

## Contributing

Contributions are welcome! Please feel free to submit a pull request or open an issue for any suggestions or improvements.

## License

This project is licensed under the MIT License. See the LICENSE file for more details.
# Digital notebook

Digital notebook is a collection of Python scripts and associated Jupyter notebooks that are meant to help with tasks such as biological sequence curation, alignment, and phylogenetic tree building.

You can either call the Python scripts directly (they sit in the src/ folder) or you can use or create new Jupyter notebooks. 

There are examples of vanilla workflows to be used with your own data in the /workflows folder.

This is a repository that contains the files so they can be used in local deployment - however the files will also soon be available to run entirely online.

## Dependencies for local deployment

1. Python 3
2. Jupyter Notebook
3. Python modules as described in install/requirements.txt

## Installation for local deployment

1. Clone this directory
2. [Install Jupyter Notebook](https://jupyter.readthedocs.io/en/latest/install.html "Install Jupyter Notebook")
3. Install the modules needed to your Python installation `pip3 install -r requirements.txt`
4. Start Jupyter Notebook from the command line in this directory using `jupyter notebook` command
 

## Usage

1. Add the files you want to use to the files directory (you don't have to, but it makes navigation of files easier)
2. Open the Jupyter Notebook you wish to use and load your files by providing the full filepath in one of the first cells in the notebook
3. Run the cells to execute the Python code contained within them


## Contributing

1. Fork it!
2. Create your feature branch: `git checkout -b my-new-feature`
3. Commit your changes: `git commit -am 'Add some feature'`
4. Push to the branch: `git push origin my-new-feature`
5. Submit a pull request :D


## Credits

Gabe Foley

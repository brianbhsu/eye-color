# Make a Code Genome App

## Basic Steps

1. Write code to analyze a VCF file
2. Run the Code Genome App using `run_genome_app.py`
3. Create a `README` for the Code Genome App

## Walkthrough

Let's make a Code Genome App to list your chromosomes and count the number of variants you have in each chromosome.

### 1. Write code to analyze a VCF file

* The `list_chromosomes` function in the `list_chromosomes.py` file opens the VCF file and counts the number of
variants in each chromosome.

* Results are written to the `genome-app-name.output.ga` file in our
[GA file format](https://github.com/Guardiome/create-genome-app/blob/master/simple-genome-app/tools/README.md#ga-file-format).
    * The output MUST be written in our GA file format.

You have just created a Code Genome App. Now you are ready to run the List Chromosomes Genome App!

**You can use ANY code and input files to make your Code Genome App**

### 2. Run the Code Genome App

Modify the `run_genome_app.py` file to call your code.

```sh
from list_chromosomes import list_chromosomes
list_chromosomes()
```

Run the List Chromosomes Genome App with the following command in the `tools` folder:

```sh
python run_genome_app.py
```

It will create a `genome-app-name.output.ga` file inside the output folder.

### 3. Create a README for the Code Genome App

Create a `README` file for the Code Genome App by either modifying this
[file](https://github.com/Guardiome/create-genome-app/blob/master/code-genome-app/README.md) or making your own README.
The README should be in the genome-app-name folder.

The `README` should contain basic information about the Code Genome App in general regardless of the result. Here is
a list of information about the Code Genome App that you may want to add:

* Genome App Name
* Description
* Biology
* Genetics
* Risk Factors

The README created for each Code Genome App will be shown with the result in Guardiome AI.

<div>
    <img src="./simple-genome-app/media/guardiome-logo.png" align="center" width=400 height=150>
</div>

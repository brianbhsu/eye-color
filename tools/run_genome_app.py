"""
Genome App entry point.
File for running Genome App from 1) Genome AI & 2) command line.
"""


def run_genome_app():
    """
    Required function for Guardiome to run this Genome App.
    :return: None
    """

    # Import any function from any file
    from detect_eye_color import detect_eye_color

    # Run any function
    detect_eye_color()


# This Genome App can also be ran from comman line: $ python run_genome_app.py
if __name__ == '__main__':

    run_genome_app()

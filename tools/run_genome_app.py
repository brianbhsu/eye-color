"""
Genome App entry point.
File for running Genome App from either:
    1) Genome AI or
    2) command line ($ python run_genome_app.py).
"""


def run_genome_app():
    """
    Required function for Genome AI to run this Genome App. This Genome App is
        responsible for producing either:
            1) <genome-app-repository>/output/output.json or
            2) <genome-app-repository>/output/output.g2p.
    Arguments:
        None
    Returns:
        None
    """

    from detect_eye_color import detect_eye_color

    detect_eye_color()


if __name__ == '__main__':

    run_genome_app()

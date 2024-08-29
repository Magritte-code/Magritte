from vanZadelhoff_1_1D import create_model as vanZadelhoff_1D_setup
from vanZadelhoff_1_1D import run_model as vanZadelhoff_1D_run
from vanZadelhoff_1_3D_healpix import create_model as vanZadelhoff_3D_setup
from vanZadelhoff_1_3D_healpix import run_model as vanZadelhoff_3D_run

from import_phantom_3D_model import import_phantom
from reduce_phantom_3D_model import reduce_phantom
from run_phantom_3D_model_reduced import run_model as run_phantom

from read_shuffled_lines_one_line_approx import create_model as shuffled_lines_setup
from read_shuffled_lines_one_line_approx import run_model as shuffled_lines_run

import pytest

# @contextmanager
# def does_not_raise():
#     yield

class TestNumeric:

    #incremental testing; see conftest.py or pytest.ini
    @pytest.mark.incremental
    class TestShuffledLines:
        def test_shuffledlines_setup(self):
            shuffled_lines_setup()

        def test_shuffledlines_run(self):
            assert shuffled_lines_run()


    @pytest.mark.incremental
    class TestVanZadelhoff1D:
        def test_vanZadelhoff1D_setup(self):
            vanZadelhoff_1D_setup('a')

        def test_vanZadelhoff1D_run(self):
            assert vanZadelhoff_1D_run('a', nosave=True)

    #for testing whether a 3D model can succesfully converge to the correct solution
    @pytest.mark.incremental
    class TestVanZadelhoff3D:
        def test_vanZadelhoff3D_setup(self):
            vanZadelhoff_3D_setup('a')

        def test_vanZadelhoff3D_run(self):
            assert vanZadelhoff_3D_run('a', nosave=True)

    #for testing whether a phantom model can be correctly imported, reduced and run.
    @pytest.mark.incremental
    class TestPhantom:
        #note: import has some extra filtering: also removes points with zero or negative abundances
        #testing whether the setup runs
        def test_import_phantom(self):
            import_phantom()

        #reduction will probably take the most time of these tests
        def test_reduce_phantom(self):
            reduce_phantom()

        #currently only checks whether it is able to run succesfully, does not yet check whether the results make sense
        def test_run_phantom_model(self):
            assert run_phantom()



    # @pytest.mark.incremental
    # class TestReducePhantom:
    #     def test_density_distribution1D_setup(self):
    #         density_dist_setup('a')
    #
    #     def test_density_distribution1D_run(self):
    #         assert density_dist_run('a', nosave=True)
    #
    # @pytest.mark.incremental
    # class TestRunPhantom:
    #     def test_velocity_gradient_1D_setup(self):
    #         velocity_gradient_setup()
    #
    #     def test_velocity_gradient_1D_run(self):
    #         assert velocity_gradient_run(nosave=True)

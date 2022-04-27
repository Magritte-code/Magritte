from import_phantom_3D_model import import_phantom
from reduce_phantom_3D_model import reduce_phantom
from run_phantom_3D_model_reduced import run_model as run_phantom
# from all_constant_single_ray import run_model as all_constant_run
# from density_distribution_1D import create_model as density_dist_setup
# from density_distribution_1D import run_model as density_dist_run
# from constant_velocity_gradient_1D import create_model as velocity_gradient_setup
# from constant_velocity_gradient_1D import run_model as velocity_gradient_run

import pytest

# @contextmanager
# def does_not_raise():
#     yield

class TestNumeric:
    #incremental testing; see conftest.py or pytest.ini
    @pytest.mark.incremental
    class TestImportPhantom:
        #note: import has some extra filtering: also removes points with zero or negative abundances
        #testing whether the setup runs
        def test_import_phantom(self):
            import_phantom()

    #reduction will probably take the most time of these tests
    @pytest.mark.incremental
    class TestReducePhantom:
        def test_reduce_phantom(self):
            reduce_phantom()

    #currently only checks whether it is able to run succesfully, does not yet check whether the results make sense
    @pytest.mark.incremental
    class TestRunPhantomModel:
        def test_run_phantom_model(self):
            run_phantom()



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

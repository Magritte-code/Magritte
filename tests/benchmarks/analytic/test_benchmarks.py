from all_constant_single_ray import create_model as all_constant_setup
from all_constant_single_ray import run_model as all_constant_run
from density_distribution_1D import create_model as density_dist_setup
from density_distribution_1D import run_model as density_dist_run
from constant_velocity_gradient_1D import create_model as velocity_gradient_setup
from constant_velocity_gradient_1D import run_model as velocity_gradient_run

import pytest

# @contextmanager
# def does_not_raise():
#     yield

class TestAnalytic:
    #incremental testing; see conftest.py or pytest.ini
    @pytest.mark.incremental
    class TestAllConstant:
        #testing whether the setup runs
        def test_all_constant_setup(self):
            all_constant_setup()

        #testing whether we get the correct result
        def test_all_constant_run(self):
            assert all_constant_run(nosave=True)

    @pytest.mark.incremental
    class TestDensityDistribution1D:
        def test_density_distribution1D_setup(self):
            density_dist_setup('a')

        def test_density_distribution1D_run(self):
            assert density_dist_run('a', nosave=True)

    @pytest.mark.incremental
    class TestVelocityGradient1D:
        def test_velocity_gradient_1D_setup(self):
            velocity_gradient_setup()

        def test_velocity_gradient_1D_run(self):
            assert velocity_gradient_run(nosave=True)

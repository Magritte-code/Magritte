from all_constant_single_ray import create_model as all_constant_setup
from all_constant_single_ray import run_model as all_constant_run
from all_zero_single_ray import create_model as all_zero_setup
from all_zero_single_ray import run_model as all_zero_run
from density_distribution_1D import create_model as density_dist_setup
from density_distribution_1D import run_model as density_dist_run
from density_distribution_1D_image import create_model as density_dist_image_setup
from density_distribution_1D_image import run_model as density_dist_image_run
from constant_velocity_gradient_1D import create_model as velocity_gradient_setup
from constant_velocity_gradient_1D import run_model as velocity_gradient_run
from constant_velocity_gradient_1D_image import create_model as velocity_gradient_image_setup
from constant_velocity_gradient_1D_image import run_model as velocity_gradient_image_run

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
    class TestAllZero:
        #testing whether the setup runs
        def test_all_zero_setup(self):
            all_zero_setup()

        #testing whether we get the correct result
        def test_all_zero_run(self):
            assert all_zero_run(nosave=True)


    @pytest.mark.incremental
    class TestDensityDistribution1D:
        def test_density_distribution1D_setup(self):
            density_dist_setup('a')

        def test_density_distribution1D_run(self):
            assert density_dist_run('a', nosave=True, use_widgets=False)

    @pytest.mark.incremental
    class TestDensityDistribution1DImage:
        def test_density_distribution1D_image_setup(self):
            density_dist_image_setup('a')

        def test_density_distribution1D_image_run(self):
            assert density_dist_image_run('a', nosave=True, use_widgets=False)

    @pytest.mark.incremental
    class TestVelocityGradient1D:
        def test_velocity_gradient_1D_setup(self):
            velocity_gradient_setup()

        class TestVelocityGradient1DBenchmarks:
            def test_velocity_gradient_1D_run_bench1(self):
                assert velocity_gradient_run(nosave=True, benchindex=1, use_widgets=False)
            def test_velocity_gradient_1D_run_bench2(self):
                assert velocity_gradient_run(nosave=True, benchindex=2, use_widgets=False)
            def test_velocity_gradient_1D_run_bench3(self):
                assert velocity_gradient_run(nosave=True, benchindex=3, use_widgets=False)

    @pytest.mark.incremental
    class TestVelocityGradient1DImage:
        def test_velocity_gradient_1D_image_setup(self):
            velocity_gradient_image_setup()

        class TestVelocityGradient1DImageBenchmarks:
            def test_velocity_gradient_1D_image_run_bench1(self):
                assert velocity_gradient_image_run(nosave=True, benchindex=1, use_widgets=False)
            def test_velocity_gradient_1D_image_run_bench2(self):
                assert velocity_gradient_image_run(nosave=True, benchindex=2, use_widgets=False)
            def test_velocity_gradient_1D_image_run_bench3(self):
                assert velocity_gradient_image_run(nosave=True, benchindex=3, use_widgets=False)



#TODO ADD ALL ANALYTIC BENCHMARKS, get some manner of performance somewhere else

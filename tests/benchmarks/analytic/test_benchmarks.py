from reject_nonexistent_species import create_model as reject_nonexistent_species_setup
from all_constant_single_ray import create_model as all_constant_setup
from all_constant_single_ray import run_model as all_constant_run
from all_constant_single_ray_testuv import create_model as all_constant_testuv_setup
from all_constant_single_ray_testuv import run_model as all_constant_testuv_run
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
from constant_velocity_gradient_1D_new_imager import create_model as velocity_gradient_new_imager_setup
from constant_velocity_gradient_1D_new_imager import run_model as velocity_gradient_new_imager_run



import pytest

# @contextmanager
# def does_not_raise():
#     yield

#In order to check whether our Io module works correctly
class TestInput:
    def test_read_nonexistent_file(self):
        import magritte.core as magritte
        with pytest.raises(Exception):
            bogusFile = "dlpbmfhqbzeszvroptsqwupklhbtvkzkmebzdduhetthdjakgz.hdf5"
            model = magritte.Model(bogusFile) #reading bogus files should raise an exception
            #note: if we segfault, pytest will also complain, so this should protect against segfaults

    def test_read_nonexistent_species(self):
        with pytest.raises(Exception):
            reject_nonexistent_species_setup()

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
    class TestUVFeautrier:
        def test_uv_feautrier_setup(self):
            all_constant_testuv_setup()

        def test_uv_feautrier_run(self):
            assert all_constant_testuv_run(nosave=True)

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


#To check whether the new imager works on 1D models
class TestNewImager:
    @pytest.mark.incremental
    class TestVelocityGradient1DNewImager:
        def test_velocity_gradient_1D_new_imager_setup(self):
            velocity_gradient_new_imager_setup()

        class TestVelocityGradient1DNewImagerBenchmarks:
            def test_velocity_gradient_1D_new_imager_run_bench1(self):
                assert velocity_gradient_new_imager_run(nosave=True, benchindex=1, use_widgets=False)
            def test_velocity_gradient_1D_new_imager_run_bench2(self):
                assert velocity_gradient_new_imager_run(nosave=True, benchindex=2, use_widgets=False)
            def test_velocity_gradient_1D_new_imager_run_bench3(self):
                assert velocity_gradient_new_imager_run(nosave=True, benchindex=3, use_widgets=False)

    @pytest.mark.incremental
    class TestVelocityGradient3DNewImager:
        def test_velocity_gradient_3D_new_imager_setup(self):
            velocity_gradient_new_imager_setup()

        def test_velocity_gradient_3D_new_imager_run(self):
            assert velocity_gradient_new_imager_run(nosave=True, use_widgets=False)


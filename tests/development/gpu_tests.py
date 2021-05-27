import numpy         as np
import magritte.core as magritte

model = magritte.Model('../models/density_distribution_VZa_single_ray.hdf5')

# model.compute_spectral_discretisation ()
# model.compute_inverse_line_widths     ()
# model.compute_LTE_level_populations   ()

magritte.list_accelerators()


model.set()

print(np.array(model.add()))

# print('Computing radiation field...')
# model.compute_radiation_field_feautrier_order_2 ()
# print('Done!')
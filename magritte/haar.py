import numpy as np
import numba


# Characterisation of the 3D Haar wavelets
hwav = np.ones((7,2,2,2), dtype=np.int64)

hwav[0,0,0,1] = -1
hwav[0,0,1,1] = -1
hwav[0,1,0,1] = -1
hwav[0,1,1,1] = -1

hwav[1,0,1,0] = -1
hwav[1,0,1,1] = -1
hwav[1,1,1,0] = -1
hwav[1,1,1,1] = -1

hwav[2,1,0,0] = -1
hwav[2,1,0,1] = -1
hwav[2,1,1,0] = -1
hwav[2,1,1,1] = -1

hwav[3,1,0,0] = -1
hwav[3,1,0,1] = -1
hwav[3,0,1,0] = -1
hwav[3,0,1,1] = -1

hwav[4,0,0,1] = -1
hwav[4,1,0,1] = -1
hwav[4,0,1,0] = -1
hwav[4,1,1,0] = -1

hwav[5,1,0,0] = -1
hwav[5,1,1,0] = -1
hwav[5,0,0,1] = -1
hwav[5,0,1,1] = -1

hwav[6,1,0,0] = -1
hwav[6,0,1,0] = -1
hwav[6,0,0,1] = -1
hwav[6,1,1,1] = -1


class Haar():

    def __init__(self, points, q):
        """
        Constructor for a Haar object.

        Parameters
        ----------
        points    : n by 3 array
        q         : Maximum level of refinement (side of the cube is 2**q, number of volume elements 8**q)
        """
        # Store the maximum (finenst) level
        self.q = q

        # Compute size of the box
        self.xyz_min = np.min(points, axis=0)
        self.xyz_max = np.max(points, axis=0)

        # Set a tolerance (which should be smaller than the resolution of the cube)
        tol = 1.0e-3 / 2**self.q

        # Compute xyz size of the box (1 tol'th larger)
        self.xyz_L = (1.0 + tol) * (self.xyz_max - self.xyz_min)

        # Normalize the point coordinates to be in [0,1]
        self.points_normed = (points - self.xyz_min) / self.xyz_L

        # Count the number of points in each octree grid cell
        self.num = Haar.get_number_matrices(self.points_normed, self.q)

        # Compute cube (ix, iy, iz) indices in the linear representation
        self.ids = [np.array(Haar.__cub__(np.arange(8**k, dtype=np.int64)), dtype=np.int64).T for k in range(q)]


        # # Extract the cube indices (Is)
        # self.Is = Haar.__lin__(self.indices[:,0], self.indices[:,1], self.indices[:,2])

        # # Get a convenient ordering by ordering cube indices (Is)
        # self.order = np.argsort(self.Is)

        # # Get the unique indices, a map from the point indices to these unique indices, and a count of each unique index
        # self.Is_unique, self.Is_omap, self.Is_count = np.unique(self.Is[self.order], return_inverse=True, return_counts=True)


    def map_data(self, data, interpolate=True):
        """
        Map data to the cube.
        """
        # Integrate (i.e. sum) the point data over each cell
        dat = Haar.integrate_data(self.points_normed, self.q, data)
        # Divide by the numberof points in each cell to get the average
        for k in range(self.q):
            dat[k] = np.divide(dat[k], self.num[k], out=dat[k], where=(self.num[k]!=0))
        # Interpolate empty cells in each cube
        if interpolate:
            for k in range(1, self.q):
                dat[k] = Haar.interpolate(self.ids[k], self.num[k], dat[k], dat[k-1], k)
        # Return a hierarchical representation of the data
        return dat


    @staticmethod
    @numba.njit(parallel=True)
    def get_number_matrices(points_normed, q):
        """
        Compute the number of points that live in each cell at each level.
        """
        # Create number matrices
        num = [np.zeros((2**k, 2**k, 2**k), dtype=np.int64) for k in range(q)]
        # For all levels
        for k in range(q):
            # Compute the indices of the points in the octree grid
            indices = (points_normed * 2**k).astype(np.int64)
            # Count the number of points at every index
            for ix, iy, iz in indices:
                num[k][ix, iy, iz] += 1
        # Return the number matrices
        return num


    @staticmethod
    @numba.njit(parallel=True)
    def integrate_data(points_normed, q, data):
        """
        Sum the data in each cell at each level.
        """
        # Create number matrices
        dat = [np.zeros((2**k, 2**k, 2**k), dtype=data.dtype) for k in range(q)]
        # For all levels
        for k in range(q):
            # Compute the indices of the points in the octree grid
            indices = (points_normed * 2**k).astype(np.int64)
            # Count the number of points at every index
            for i, (ix, iy, iz) in enumerate(indices):
                dat[k][ix, iy, iz] += data[i]
        # Return the integrated data
        return dat


    @staticmethod
    @numba.njit(parallel=True)
    def interpolate(ids, num, dat, dat_up, k):
        # For all cells in the cube at this level
        for i in numba.prange(8**k):
            # Get the cube triple index
            ix, iy, iz = ids[i]
            # If the cell at this level is empty take the value from the above level
            if (num[ix, iy, iz] == 0):
                dat[ix, iy, iz] = dat_up[ix//2, iy//2, iz//2]
        # Return data cube
        return dat

    def generate_points(self, avg, wav, threshold=1.0e-3):
        """
        Generate points sampling the avg and wav data. Also generates corresponding boundary points.
        """
        points = []
        bulk = []
        #start with 8 edge points of the cube as boundary
        boundary = [
            [self.xyz_min + np.array([0,0,0]) * self.xyz_L],
            [self.xyz_min + np.array([0,0,1]) * self.xyz_L],
            [self.xyz_min + np.array([0,1,0]) * self.xyz_L],
            [self.xyz_min + np.array([0,1,1]) * self.xyz_L],
            [self.xyz_min + np.array([1,0,0]) * self.xyz_L],
            [self.xyz_min + np.array([1,0,1]) * self.xyz_L],
            [self.xyz_min + np.array([1,1,0]) * self.xyz_L],
            [self.xyz_min + np.array([1,1,1]) * self.xyz_L]]

        for k in range(self.q-1):
            pts = np.argwhere(np.sum(wav[k]/avg[k] > threshold, axis=0))
            bulk_pts = (pts + 0.5) * (self.xyz_L / 2**k) + self.xyz_min
            bulk.append(bulk_pts)
            # Get subset of points on the edge, add projected versions to boundary
            MIN_IDX = 0
            MAX_IDX = 2**k - 1
            ## get x=const plane
            x_min_pts = pts[np.where(pts[:,0]==MIN_IDX)]#TODO MIGHT CONTAIN ERROR
            x_max_pts = pts[np.where(pts[:,0]==MAX_IDX)]
            # projecting onto the position of the bounds
            x_min_pts = (x_min_pts + 0.5) * (self.xyz_L / 2**k) + self.xyz_min
            x_max_pts = (x_max_pts + 0.5) * (self.xyz_L / 2**k) + self.xyz_min
            x_min_pts[:, 0] = self.xyz_min[0]
            x_max_pts[:, 0] = self.xyz_min[0] + self.xyz_L[0]
            boundary.append(x_min_pts)
            boundary.append(x_max_pts)
            ## y = const plane
            y_min_pts = pts[np.where(pts[:,1]==MIN_IDX)]#TODO MIGHT CONTAIN ERROR
            y_max_pts = pts[np.where(pts[:,1]==MAX_IDX)]
            # projecting onto the bounds
            y_min_pts = (y_min_pts + 0.5) * (self.xyz_L / 2**k) + self.xyz_min
            y_max_pts = (y_max_pts + 0.5) * (self.xyz_L / 2**k) + self.xyz_min
            y_min_pts[:, 1] = self.xyz_min[1]
            y_max_pts[:, 1] = self.xyz_min[1] + self.xyz_L[1]
            boundary.append(y_min_pts)
            boundary.append(y_max_pts)
            ## z = const plane
            z_min_pts = pts[np.where(pts[:,2]==MIN_IDX)]#TODO MIGHT CONTAIN ERROR
            z_max_pts = pts[np.where(pts[:,2]==MAX_IDX)]
            # projecting onto the bounds
            z_min_pts = (z_min_pts + 0.5) * (self.xyz_L / 2**k) + self.xyz_min
            z_max_pts = (z_max_pts + 0.5) * (self.xyz_L / 2**k) + self.xyz_min
            z_min_pts[:, 2] = self.xyz_min[2]
            z_max_pts[:, 2] = self.xyz_min[2] + self.xyz_L[2]
            boundary.append(z_min_pts)
            boundary.append(z_max_pts)

        bulk = np.concatenate(bulk)
        boundary = np.concatenate(boundary)

        # get number of non-boundary points to return
        nb_bulk_points = len(bulk)
        nb_boundary_points = len(boundary)
        print("number of bulk points: ", nb_bulk_points)
        print("number of boundary points: ", nb_boundary_points)

        # Return a numpy array of points
        return (np.concatenate((boundary, bulk)), nb_boundary_points)

    # def map_data_sparse(self, data):
    #     """
    #     Map data to the cube using a sparse representation.
    #     The resulting data is indexed by Is_unique.
    #     """
    #     return Haar.__map_data__(data, self.order, self.Is_unique, self.Is_omap, self.Is_count)


    # @staticmethod
    # @numba.njit(parallel=True)
    # def __map_data_sparse__(data, order, Is_unique, Is_omap, Is_count):
    #     """
    #     Jitted version of the map_data function.
    #     """
    #     data_mapped = np.zeros((Is_unique.shape[0],) + data.shape[1:], dtype=data.dtype)
    #     data        = data[order]
    #
    #     for i in numba.prange(data.shape[0]):
    #         data_mapped[Is_omap[i]] = data_mapped[Is_omap[i]] + data[i]
    #
    #     for i in numba.prange(data_mapped.shape[0]):
    #         data_mapped[i] = data_mapped[i] / Is_count[i]
    #
    #     return data_mapped


    @staticmethod
    @numba.njit(parallel=True)
    def __lin__(i, j, k):
        """
        If the binary representation of i, j, and k are given by:
            bin(i) = (..., i8, i4, i2, i1)
            bin(j) = (..., j8, j4, j2, j1)
            bin(k) = (..., k8, k4, k2, k1)
        this function returns the number, r, for which the binary repressentation is given by:
            bin(r) = (..., k8, j8, i8, k4, j4, i4, k2, j2, i2, k1, j1, i1)
        inducing an hierarchical ordering.
        """

        r = 0

        j = j << 1
        k = k << 2

        r = r + (i & 2**0)
        r = r + (j & 2**1)
        r = r + (k & 2**2)

        # Yes, this can be rewritten as a for-loop
        # Done this way since I am not sure if the compiler unrolls loops.
        # Should be tested!

        i = i << 2
        j = j << 2
        k = k << 2

        r = r + (i & 2**3)
        r = r + (j & 2**4)
        r = r + (k & 2**5)

        i = i << 2
        j = j << 2
        k = k << 2

        r = r + (i & 2**6)
        r = r + (j & 2**7)
        r = r + (k & 2**8)

        i = i << 2
        j = j << 2
        k = k << 2

        r = r + (i & 2**9)
        r = r + (j & 2**10)
        r = r + (k & 2**11)

        i = i << 2
        j = j << 2
        k = k << 2

        r = r + (i & 2**12)
        r = r + (j & 2**13)
        r = r + (k & 2**14)

        i = i << 2
        j = j << 2
        k = k << 2

        r = r + (i & 2**15)
        r = r + (j & 2**16)
        r = r + (k & 2**17)

        i = i << 2
        j = j << 2
        k = k << 2

        r = r + (i & 2**18)
        r = r + (j & 2**19)
        r = r + (k & 2**20)

        i = i << 2
        j = j << 2
        k = k << 2

        r = r + (i & 2**21)
        r = r + (j & 2**22)
        r = r + (k & 2**23)

        i = i << 2
        j = j << 2
        k = k << 2

        r = r + (i & 2**24)
        r = r + (j & 2**25)
        r = r + (k & 2**26)

        i = i << 2
        j = j << 2
        k = k << 2

        r = r + (i & 2**27)
        r = r + (j & 2**28)
        r = r + (k & 2**29)

        i = i << 2
        j = j << 2
        k = k << 2

        r = r + (i & 2**30)
        r = r + (j & 2**31)
        r = r + (k & 2**32)

        return r


    @staticmethod
    @numba.njit(parallel=True)
    def __cub__(r):
        """
        If the binary representation of r is given by:
            bin(r) = (..., k8, j8, i8, k4, j4, i4, k2, j2, i2, k1, j1, i1)
        this function returns the numbers, i, j, and k, for which the binary repressentations are given by:
            bin(i) = (..., i8, i4, i2, i1)
            bin(j) = (..., j8, j4, j2, j1)
            bin(k) = (..., k8, k4, k2, k1)
        inducing an hierarchical ordering.
        """

        i = 0
        j = 0
        k = 0

        # Yes, this can be rewritten as a for-loop
        # Done this way since I am not sure if the compiler unrolls loops.
        # Should be tested!

        i = i + (r & 2**0)
        r = r >> 1
        j = j + (r & 2**0)
        r = r >> 1
        k = k + (r & 2**0)

        i = i + (r & 2**1)
        r = r >> 1
        j = j + (r & 2**1)
        r = r >> 1
        k = k + (r & 2**1)

        i = i + (r & 2**2)
        r = r >> 1
        j = j + (r & 2**2)
        r = r >> 1
        k = k + (r & 2**2)

        i = i + (r & 2**3)
        r = r >> 1
        j = j + (r & 2**3)
        r = r >> 1
        k = k + (r & 2**3)

        i = i + (r & 2**4)
        r = r >> 1
        j = j + (r & 2**4)
        r = r >> 1
        k = k + (r & 2**4)

        i = i + (r & 2**5)
        r = r >> 1
        j = j + (r & 2**5)
        r = r >> 1
        k = k + (r & 2**5)

        i = i + (r & 2**6)
        r = r >> 1
        j = j + (r & 2**6)
        r = r >> 1
        k = k + (r & 2**6)

        i = i + (r & 2**7)
        r = r >> 1
        j = j + (r & 2**7)
        r = r >> 1
        k = k + (r & 2**7)

        i = i + (r & 2**8)
        r = r >> 1
        j = j + (r & 2**8)
        r = r >> 1
        k = k + (r & 2**8)

        i = i + (r & 2**9)
        r = r >> 1
        j = j + (r & 2**9)
        r = r >> 1
        k = k + (r & 2**9)

        i = i + (r & 2**10)
        r = r >> 1
        j = j + (r & 2**10)
        r = r >> 1
        k = k + (r & 2**10)

        i = i + (r & 2**11)
        r = r >> 1
        j = j + (r & 2**11)
        r = r >> 1
        k = k + (r & 2**11)

        i = i + (r & 2**12)
        r = r >> 1
        j = j + (r & 2**12)
        r = r >> 1
        k = k + (r & 2**12)

        i = i + (r & 2**13)
        r = r >> 1
        j = j + (r & 2**13)
        r = r >> 1
        k = k + (r & 2**13)

        i = i + (r & 2**14)
        r = r >> 1
        j = j + (r & 2**14)
        r = r >> 1
        k = k + (r & 2**14)

        i = i + (r & 2**15)
        r = r >> 1
        j = j + (r & 2**15)
        r = r >> 1
        k = k + (r & 2**15)

        return i, j, k


    def get_Haar_wavelets(self, data_cube):
        """
        3D Haar transform of data.
        """
        averages = [np.zeros((   2**k, 2**k, 2**k), dtype=data_cube[0].dtype) for k in range(self.q  )]
        wavelets = [np.zeros((7, 2**k, 2**k, 2**k), dtype=data_cube[0].dtype) for k in range(self.q-1)]

        averages[-1] = data_cube[-1]

        for k in reversed(range(self.q-1)):
            averages[k], wavelets[k] = Haar.__get_Haar_wavelets__(self.ids[k], self.ids[k+1], averages[k], averages[k+1], wavelets[k], k)

        return averages, wavelets

    @staticmethod
    @numba.njit(parallel=True)
    def __get_Haar_wavelets__(ids, ids_down, avg, dat, wav, k):
        """
        3D Haar transform on data given i.
        """
        # For all cells in the cube at this level
        for p in numba.prange(8**k):
            P = 8*p
            # Get the cube triple index
            i0 = ids_down[P  ][0], ids_down[P  ][1], ids_down[P  ][2]
            i1 = ids_down[P+1][0], ids_down[P+1][1], ids_down[P+1][2]
            i2 = ids_down[P+2][0], ids_down[P+2][1], ids_down[P+2][2]
            i3 = ids_down[P+3][0], ids_down[P+3][1], ids_down[P+3][2]
            i4 = ids_down[P+4][0], ids_down[P+4][1], ids_down[P+4][2]
            i5 = ids_down[P+5][0], ids_down[P+5][1], ids_down[P+5][2]
            i6 = ids_down[P+6][0], ids_down[P+6][1], ids_down[P+6][2]
            i7 = ids_down[P+7][0], ids_down[P+7][1], ids_down[P+7][2]

            j = ids[p][0], ids[p][1], ids[p][2]

            avg   [j] = 0.125 * (dat[i0] + dat[i1] + dat[i2] + dat[i3] + dat[i4] + dat[i5] + dat[i6] + dat[i7])
            wav[0][j] = 0.125 * (dat[i0] + dat[i1] + dat[i2] + dat[i3] - dat[i4] - dat[i5] - dat[i6] - dat[i7])
            wav[1][j] = 0.125 * (dat[i0] + dat[i1] - dat[i2] - dat[i3] + dat[i4] + dat[i5] - dat[i6] - dat[i7])
            wav[2][j] = 0.125 * (dat[i0] - dat[i1] + dat[i2] - dat[i3] + dat[i4] - dat[i5] + dat[i6] - dat[i7])
            wav[3][j] = 0.125 * (dat[i0] - dat[i1] - dat[i2] + dat[i3] + dat[i4] - dat[i5] - dat[i6] + dat[i7])
            wav[4][j] = 0.125 * (dat[i0] + dat[i1] - dat[i2] - dat[i3] - dat[i4] - dat[i5] + dat[i6] + dat[i7])
            wav[5][j] = 0.125 * (dat[i0] - dat[i1] + dat[i2] - dat[i3] - dat[i4] + dat[i5] - dat[i6] + dat[i7])
            wav[6][j] = 0.125 * (dat[i0] - dat[i1] - dat[i2] + dat[i3] - dat[i4] + dat[i5] + dat[i6] - dat[i7])

        return avg, wav


    def wavelet_approx(self, averages, wavelets):
        """
        Wavelet approximation of the data.
        """
        approx = [np.zeros((2**k, 2**k, 2**k), dtype=averages[0].dtype) for k in range(self.q)]

        approx[0] = averages[0]

        for k in range(1, self.q):
            approx[k] = Haar.__wavelet_approx__(self.ids[k-1], self.ids[k], wavelets[k-1], approx[k-1], approx[k], k)

        return approx


    @staticmethod
    @numba.njit(parallel=True)
    def __wavelet_approx__(ids_up, ids, wav, app_up, app, k):
        """
        3D Haar wavelet approximation on data given i.
        """
        # For all cells in the cube at this level
        for p in numba.prange(8**(k-1)):
            P = 8*p
            # Get the cube triple index
            i0 = ids[P  ][0], ids[P  ][1], ids[P  ][2]
            i1 = ids[P+1][0], ids[P+1][1], ids[P+1][2]
            i2 = ids[P+2][0], ids[P+2][1], ids[P+2][2]
            i3 = ids[P+3][0], ids[P+3][1], ids[P+3][2]
            i4 = ids[P+4][0], ids[P+4][1], ids[P+4][2]
            i5 = ids[P+5][0], ids[P+5][1], ids[P+5][2]
            i6 = ids[P+6][0], ids[P+6][1], ids[P+6][2]
            i7 = ids[P+7][0], ids[P+7][1], ids[P+7][2]

            j = ids_up[p][0], ids_up[p][1], ids_up[p][2]

            app[i0] = app_up[j] + wav[0][j] + wav[1][j] + wav[2][j] + wav[3][j] + wav[4][j] + wav[5][j] + wav[6][j]
            app[i1] = app_up[j] + wav[0][j] + wav[1][j] - wav[2][j] - wav[3][j] + wav[4][j] - wav[5][j] - wav[6][j]
            app[i2] = app_up[j] + wav[0][j] - wav[1][j] + wav[2][j] - wav[3][j] - wav[4][j] + wav[5][j] - wav[6][j]
            app[i3] = app_up[j] + wav[0][j] - wav[1][j] - wav[2][j] + wav[3][j] - wav[4][j] - wav[5][j] + wav[6][j]
            app[i4] = app_up[j] - wav[0][j] + wav[1][j] + wav[2][j] + wav[3][j] - wav[4][j] - wav[5][j] - wav[6][j]
            app[i5] = app_up[j] - wav[0][j] + wav[1][j] - wav[2][j] - wav[3][j] - wav[4][j] + wav[5][j] + wav[6][j]
            app[i6] = app_up[j] - wav[0][j] - wav[1][j] + wav[2][j] - wav[3][j] + wav[4][j] - wav[5][j] + wav[6][j]
            app[i7] = app_up[j] - wav[0][j] - wav[1][j] - wav[2][j] + wav[3][j] + wav[4][j] + wav[5][j] - wav[6][j]

        return app


    def wavelet_approx_threshold(self, averages, wavelets, thres_frac=0.01):
        """
        Wavelet approximation of the data.
        """
        approx = [np.zeros((2**k, 2**k, 2**k), dtype=averages[0].dtype) for k in range(self.q)]
        ns     = [0                                                     for _ in range(self.q)]

        approx[0] = averages[0]
        ns    [0] = 1

        for k in range(1, self.q):
            approx[k], ns[k] = Haar.__wavelet_approx_threshold__(self.ids[k-1], self.ids[k], wavelets[k-1], approx[k-1], approx[k], k, thres_frac=thres_frac)

        return approx, ns


    @staticmethod
    @numba.njit(parallel=True)
    def __wavelet_approx_threshold__(ids_up, ids, wav, app_up, app, k, thres_frac=0.01):
        """
        3D Haar wavelet approximation on data given i.
        """

        # Note that n will not computed correctly when parallel=True
        n = 0

        # For all cells in the cube at this level
        for p in numba.prange(8**(k-1)):
            P = 8*p
            # Get the cube triple index
            i0 = ids[P  ][0], ids[P  ][1], ids[P  ][2]
            i1 = ids[P+1][0], ids[P+1][1], ids[P+1][2]
            i2 = ids[P+2][0], ids[P+2][1], ids[P+2][2]
            i3 = ids[P+3][0], ids[P+3][1], ids[P+3][2]
            i4 = ids[P+4][0], ids[P+4][1], ids[P+4][2]
            i5 = ids[P+5][0], ids[P+5][1], ids[P+5][2]
            i6 = ids[P+6][0], ids[P+6][1], ids[P+6][2]
            i7 = ids[P+7][0], ids[P+7][1], ids[P+7][2]

            j = ids_up[p][0], ids_up[p][1], ids_up[p][2]

            thres = thres_frac * app_up[j]

            app[i0] = app_up[j]
            app[i1] = app_up[j]
            app[i2] = app_up[j]
            app[i3] = app_up[j]
            app[i4] = app_up[j]
            app[i5] = app_up[j]
            app[i6] = app_up[j]
            app[i7] = app_up[j]

            # No, this cannot be rewritten as a for-loop!
            # Mind the additions and subtractions!

            if (np.abs(wav[0][j]) > thres):
                n += 1
                app[i0] += wav[0][j]
                app[i1] += wav[0][j]
                app[i2] += wav[0][j]
                app[i3] += wav[0][j]
                app[i4] -= wav[0][j]
                app[i5] -= wav[0][j]
                app[i6] -= wav[0][j]
                app[i7] -= wav[0][j]

            if (np.abs(wav[1][j]) > thres):
                n += 1
                app[i0] += wav[1][j]
                app[i1] += wav[1][j]
                app[i2] -= wav[1][j]
                app[i3] -= wav[1][j]
                app[i4] += wav[1][j]
                app[i5] += wav[1][j]
                app[i6] -= wav[1][j]
                app[i7] -= wav[1][j]

            if (np.abs(wav[2][j]) > thres):
                n += 1
                app[i0] += wav[2][j]
                app[i1] -= wav[2][j]
                app[i2] += wav[2][j]
                app[i3] -= wav[2][j]
                app[i4] += wav[2][j]
                app[i5] -= wav[2][j]
                app[i6] += wav[2][j]
                app[i7] -= wav[2][j]

            if (np.abs(wav[3][j]) > thres):
                n += 1
                app[i0] += wav[3][j]
                app[i1] -= wav[3][j]
                app[i2] -= wav[3][j]
                app[i3] += wav[3][j]
                app[i4] += wav[3][j]
                app[i5] -= wav[3][j]
                app[i6] -= wav[3][j]
                app[i7] += wav[3][j]

            if (np.abs(wav[4][j]) > thres):
                n += 1
                app[i0] += wav[4][j]
                app[i1] += wav[4][j]
                app[i2] -= wav[4][j]
                app[i3] -= wav[4][j]
                app[i4] -= wav[4][j]
                app[i5] -= wav[4][j]
                app[i6] += wav[4][j]
                app[i7] += wav[4][j]

            if (np.abs(wav[5][j]) > thres):
                n += 1
                app[i0] += wav[5][j]
                app[i1] -= wav[5][j]
                app[i2] += wav[5][j]
                app[i3] -= wav[5][j]
                app[i4] -= wav[5][j]
                app[i5] += wav[5][j]
                app[i6] -= wav[5][j]
                app[i7] += wav[5][j]

            if (np.abs(wav[6][j]) > thres):
                n += 1
                app[i0] += wav[6][j]
                app[i1] -= wav[6][j]
                app[i2] -= wav[6][j]
                app[i3] += wav[6][j]
                app[i4] -= wav[6][j]
                app[i5] += wav[6][j]
                app[i6] += wav[6][j]
                app[i7] -= wav[6][j]

        return app, n


    def wavelet_approx_2(self, averages, wavelets):
        """
        Wavelet approximation of the data.
        """
        approx = np.zeros((2**(self.q-1), 2**(self.q-1), 2**(self.q-1)))

        return Haar.__wavelet_approx_2__(self.ids[self.q-1], averages[0], wavelets, approx, self.q)


    @staticmethod
    # @numba.njit(parallel=True)
    def __wavelet_approx_2__(ids_q, averages_0, wavelets, approx, q):
        """
        3D Haar transform on data given i.
        """

        # For all cells in the cube at the finest level
        for p in numba.prange(8**(q-1)):

            ix, iy, iz = ids_q[p]
            jx, jy, jz = ix, iy, iz

            approx[jx, jy, jz] = averages_0[0, 0, 0]

            # For the other levels levels
            for k in range(q-2, -1, -1):

                iix = ix & 1   # Extract first bit
                iiy = iy & 1   # Extract first bit
                iiz = iz & 1   # Extract first bit

                ix = ix >> 1   # Discard first bit
                iy = iy >> 1   # Discard first bit
                iz = iz >> 1   # Discard first bit

                # ((1<<(k+1))-1)) is k ones in bits

                tl = 1<<k

                dx = ( ( ids_q[p][0] & (tl-1) ) + 0.5 ) / tl
                dy = ( ( ids_q[p][1] & (tl-1) ) + 0.5 ) / tl
                dz = ( ( ids_q[p][2] & (tl-1) ) + 0.5 ) / tl

                print(k, iix,iiy,iiz, dx,dy,dz)


                # print(ids_q[p], ix,iy,iz, iix,iiy,iiz)

                for w in range(7):
                    approx[jx, jy, jz] += hwav[w, iix, iiy, iiz] * wavelets[k][w, ix, iy, iz]


        return approx


    @staticmethod
    def remesh(positions, data, q = 9, threshold = 1e-3):
        '''
        Wrapper for remeshing a point cloud (positions) based on the variations in data by thresholding the wavelet coefficients.
        The maximum depth for the wavelets is given by q.
        The returned vector of positions starts with the nb_boundary boundary points.
        '''
        # Initialize haar object, defining the rectangular wavelet grid
        haar = Haar(positions, q=q)
        # Map the data to every location
        map_data = haar.map_data(data)
        # Generate the wavelets from the mapped data, obtaining averages and coefficients for the wavelets
        data_avg, data_wav = haar.get_Haar_wavelets(map_data)
        # By thresholding the wavelet coefficients, a sparse grid of points can be obtained
        # This also generates a boundary, putting these in the nb_boundary first positions of the vector haar_points
        haar_positions, nb_boundary = haar.generate_points(data_avg, data_wav, threshold=threshold)

        return (haar_positions, nb_boundary)

import time

def remesh_recursive(positions, data, q=9, threshold= 1e-2, hullorder = 3):
    '''
    Remeshing method by simply comparing the maximal variation of the data.
    '''

    # division_order = 1#TODO: improve script such that this can be better/automatically parallelized
    # unidir_splits = 2**division_order
    # n_threads = 8**division_order#should be power of 8, as we divide the domain into 8 equal pieces # == unidir_splits**3


    # temp_positions = np.zeros((n_threads, 8*len(positions), 3))
    # bdy_positions = np.zeros((8*len(positions), 3))#also allocate way too large size for bdy positions
    new_positions = np.zeros((8*len(positions), 3))#should be large enough to contain the new positions
    remesh_nb_points = 0
    # print("remesh: ", remesh_nb_points)
    # bdy_nb_points = 0

    xyz_min = np.min(positions, axis=0)
    xyz_max = np.max(positions, axis=0)

    # xyzsplits = np.linspace(xyz_min, xyz_max, unidir_splits, endpoint=False)
    # splitids = np.arange(unidir_splits)
    # xyzdelta = (xyz_max - xyz_min) / unidir_splits
    #
    # print(xyzsplits)
    #
    # xsplits = xyzsplits[:,0]
    # ysplits = xyzsplits[:,1]
    # zsplits = xyzsplits[:,2]
    #
    #
    #
    # for xid in numba.prange(unidir_splits):
    #     for yid in numba.prange(unidir_splits):
    #         for zid in numba.prange(unidir_splits):
    #             tid = unidir_splits**2*xid+unidir_splits*yid+zid
    #             box_mincoord = np.array([xsplits[xid], ysplits[yid], zsplits[zid]])
    #             box_maxcoord = np.array([xsplits[xid], ysplits[yid], zsplits[zid]]) + xyzdelta
    #             boxindices   = np.sum((positions>=box_mincoord) & (positions<=box_maxcoord), axis=1) == 3
    #             # print(box1indices)
    #             boxpositions = positions[boxindices,:]
    #             boxdata      = data[boxindices]
    #             print(box_mincoord)
    #             print(box_maxcoord)
    #             _,_,remesh_nb_points[tid] = get_mean_remesh(boxpositions, boxdata, division_order, q, threshold, box_mincoord, box_maxcoord, temp_positions[tid,:,:], remesh_nb_points[tid])
    #
    #             # new_positions[tid,:,:].resize((remesh_nb_points[tid],3))
    #
    # #lazy concatenating is probably the simplest option
    # new_positions = np.zeros((0,3))
    # # print(new_positions.shape)
    # # print(remesh_nb_points)
    # for thread in range(n_threads):
    #     # print(remesh_nb_points[thread])
    #     # print(temp_positions[thread, 0:remesh_nb_points[thread], :])
    #     slice = temp_positions[thread, 0:remesh_nb_points[thread], :]
    #     # print(slice.shape)s
    #     new_positions = np.concatenate((new_positions, slice), axis=0)
    # sum_remesh_nb_points = np.sum(remesh_nb_points)
    # # new_positions = new_positions.reshape(sum_remesh_nb_points, 3)
    #
    # print(new_positions)

    # _,_,remesh_nb_points = get_mean_remesh(positions, data, 0, q, threshold, xyz_min, xyz_max, new_positions, remesh_nb_points)
    remesh_nb_points = get_mean_remesh_ver2(positions, data, 0, q, threshold, xyz_min, xyz_max, new_positions, remesh_nb_points)
    print("new interior points: ", remesh_nb_points)
    # #shrink positions to actual ones
    new_positions.resize((remesh_nb_points,3))

    hull = create_cubic_uniform_hull(xyz_min, xyz_max, order=hullorder)
    nb_boundary = hull.shape[0]
    print("number boundary points: ", nb_boundary)
    new_positions = np.concatenate((hull, new_positions), axis = 0)

    #TODO: add boundary points on the 8 points of the encompassing cube AND ALSO on the smaller cube-ish thing?

    # (new_positions, nb_boundary) = remesh_recursive_numba(positions, data, xyz_min, xyz_max, q=9, threshold=threshold, hullorder=hullorder)
    return (new_positions, nb_boundary)#, nb_boundary_points

# # @numba.njit(parallel=True)
# def remesh_recursive_numba(positions, data, xyz_min, xyz_max, q=9, threshold= 1e-2, hullorder = 3):
#     division_order = 1#TODO: improve script such that this can be better/automatically parallelized
#     unidir_splits = 2**division_order
#     n_threads = 8**division_order#should be power of 8, as we divide the domain into 8 equal pieces # == unidir_splits**3
#
#     print("temp_positions")
#
#     temp_positions = np.zeros((len(positions)*n_threads, 3))
#     # bdy_positions = np.zeros((8*len(positions), 3))#also allocate way too large size for bdy positions
#     # new_positions.resize(len(positions)*2)#should be large enough to contain the new positions
#     # remesh_nb_points = np.zeros((n_threads), dtype=int)
#     #separating the starting positions of the flattened array
#     start_indices = np.arange(n_threads)*len(positions)
#     remesh_nb_points = np.arange(n_threads)*len(positions)
#     print("remesh: ", remesh_nb_points)
#     # bdy_nb_points = 0
#
#     # xyz_min = np.min(positions, axis=0)
#     # xyz_max = np.max(positions, axis=0)
#
#     # xyzsplits = np.linspace(xyz_min, xyz_max, unidir_splits, endpoint=False)
#     xsplits = np.linspace(xyz_min[0], xyz_max[0], unidir_splits+1)[:-1]
#     ysplits = np.linspace(xyz_min[1], xyz_max[1], unidir_splits+1)[:-1]
#     zsplits = np.linspace(xyz_min[2], xyz_max[2], unidir_splits+1)[:-1]
#     splitids = np.arange(unidir_splits)
#     xyzdelta = (xyz_max - xyz_min) / unidir_splits
#
#     print("xsplits: ", xsplits)
#     #
#     # xsplits = xyzsplits[:,0]
#     # ysplits = xyzsplits[:,1]
#     # zsplits = xyzsplits[:,2]
#
#
#
#     for xid in numba.prange(unidir_splits):
#     # for xid in numba.prange(unidir_splits):
#         for yid in numba.prange(unidir_splits):
#         # for yid in numba.prange(unidir_splits):
#             for zid in numba.prange(unidir_splits):
#             # for zid in numba.prange(unidir_splits):
#                 tid = unidir_splits**2*xid+unidir_splits*yid+zid
#                 box_mincoord = np.array([xsplits[xid], ysplits[yid], zsplits[zid]])
#                 box_maxcoord = np.array([xsplits[xid], ysplits[yid], zsplits[zid]]) + xyzdelta
#                 boxindices   = np.sum((positions>=box_mincoord) & (positions<=box_maxcoord), axis=1) == 3
#                 # print(box1indices)
#                 boxpositions = positions[boxindices,:]
#                 boxdata      = data[boxindices]
#                 # print(box_mincoord)
#                 # print(box_maxcoord)
#                 _,_,remesh_nb_points[tid] = get_mean_remesh(boxpositions, boxdata, division_order, q, threshold, box_mincoord, box_maxcoord, temp_positions, remesh_nb_points[tid])
#
#                 # new_positions[tid,:,:].resize((remesh_nb_points[tid],3))
#
#     #lazy concatenating is probably the simplest option
#     new_positions = np.zeros((0,3))
#
#     print("newpositions")
#     # print(new_positions.shape)
#     # print(remesh_nb_points)
#     for thread in range(n_threads):
#         # print(remesh_nb_points[thread])
#         # print(temp_positions[thread, 0:remesh_nb_points[thread], :])
#         print("start: ", start_indices[thread])
#         print("end: ", remesh_nb_points[thread])
#         slice = temp_positions[start_indices[thread]:remesh_nb_points[thread], :]
#         # print(slice.shape)s
#         new_positions = np.concatenate((new_positions, slice), axis=0)
#     sum_remesh_nb_points = np.sum(remesh_nb_points)
#     # new_positions = new_positions.reshape(sum_remesh_nb_points, 3)
#
#     print(new_positions)
#
#     # print("new interior points: ", remesh_nb_points)
#     # #shrink positions to actual ones
#     # new_positions.resize((remesh_nb_points,3))
#
#     hull = create_cubic_uniform_hull(xyz_min, xyz_max, order=hullorder)
#     nb_boundary = hull.shape[0]
#     print("number boundary points: ", nb_boundary)
#     new_positions = np.concatenate((hull, new_positions), axis = 0)
#
#     #TODO: add boundary points on the 8 points of the encompassing cube AND ALSO on the smaller cube-ish thing?
#     return (new_positions, nb_boundary)#, nb_boundary_points

@numba.njit(cache=True)
def create_cubic_uniform_hull(xyz_min, xyz_max, order=3):
    # hull = np.zeros((1000, 3))#TODO: actually fill in correct value
    nx, ny, nz = (2**order+1, 2**order+1, 2**order+1)
    x_vector = np.linspace(xyz_min[0], xyz_max[0], nx)
    y_vector = np.linspace(xyz_min[1], xyz_max[1], ny)
    z_vector = np.linspace(xyz_min[2], xyz_max[2], nz)

    #x plane does not yet intersect with other planes
    xmin_plane = grid3D(np.array([xyz_min[0]]), y_vector, z_vector)
    xmax_plane = grid3D(np.array([xyz_max[0]]), y_vector, z_vector)

    #y plane intersects with x plane, so using reduced vectors for x coordinate
    ymin_plane = grid3D(x_vector[1:nx-1], np.array([xyz_min[1]]), z_vector)
    ymax_plane = grid3D(x_vector[1:nx-1], np.array([xyz_max[1]]), z_vector)

    #z plane also intersects with x plane
    zmin_plane = grid3D(x_vector[1:nx-1], y_vector[1:ny-1], np.array([xyz_min[2]]))
    zmax_plane = grid3D(x_vector[1:nx-1], y_vector[1:ny-1], np.array([xyz_max[2]]))

    #At the edges, the hull will contain duplicate points. These need to be removed
    # hull = np.unique(np.concatenate((xmin_plane, xmax_plane, ymin_plane, ymax_plane, zmin_plane, zmax_plane), axis = 0), axis = 0)
    hull = np.concatenate((xmin_plane, xmax_plane, ymin_plane, ymax_plane, zmin_plane, zmax_plane), axis = 0)

    # print(hull)
    return hull

@numba.njit(cache=True)
def grid3D(x, y, z):
    xyz = np.empty(shape=(x.size*y.size*z.size, 3))
    idx = 0
    for k in range(x.size):
        for j in range(y.size):
            for i in range(z.size):
                xyz[idx] = [x[k], y[j], z[i]]
                idx+=1
    return xyz


# @numba.jit()
# @numba.njit(cache=True, fastmath=True) #if I rewrite the .all stuff with sum < 3
# @numba.njit(cache=True, parallel=False)
@numba.njit(cache=True)
def get_mean_remesh(positions, data, depth, max_depth, threshold, min_coord, max_coord, remesh_points, remesh_nb_points):
    '''
    TODO EXPLAIN EVERY VARIABLE
    bdy_... denotes whether to add a boundary on (when adding points)
    '''

    # print("here0")
    if len(data)==0:
        # print("returning, no data")
        return (0, 0, remesh_nb_points)
    if len(data)==1 or depth==max_depth:
        # print("returning, 1 point or max depth")
        # print("len(data): ", len(data))
        # print("data: ", data)
        return (np.sum(data)/len(data), len(data), remesh_nb_points)

    delta_coord = (max_coord - min_coord) / 2.0
    #defining coordinates for sub-boxes
    def_mincoord = min_coord
    def_maxcoord = min_coord + delta_coord
    #TODO: use max by starting from true maximum (otherwise rounding errors!)

    box1_mincoord = def_mincoord + np.array([0,0,0]) * delta_coord
    box1_maxcoord = def_maxcoord + np.array([0,0,0]) * delta_coord
    box2_mincoord = def_mincoord + np.array([0,0,1]) * delta_coord
    box2_maxcoord = def_maxcoord + np.array([0,0,1]) * delta_coord
    box3_mincoord = def_mincoord + np.array([0,1,0]) * delta_coord
    box3_maxcoord = def_maxcoord + np.array([0,1,0]) * delta_coord
    box4_mincoord = def_mincoord + np.array([0,1,1]) * delta_coord
    box4_maxcoord = def_maxcoord + np.array([0,1,1]) * delta_coord
    box5_mincoord = def_mincoord + np.array([1,0,0]) * delta_coord
    box5_maxcoord = def_maxcoord + np.array([1,0,0]) * delta_coord
    box6_mincoord = def_mincoord + np.array([1,0,1]) * delta_coord
    box6_maxcoord = def_maxcoord + np.array([1,0,1]) * delta_coord
    box7_mincoord = def_mincoord + np.array([1,1,0]) * delta_coord
    box7_maxcoord = def_maxcoord + np.array([1,1,0]) * delta_coord
    box8_mincoord = def_mincoord + np.array([1,1,1]) * delta_coord
    box8_maxcoord = def_maxcoord + np.array([1,1,1]) * delta_coord
    # print("here")

    # box1indices   = np.all((positions>=box1_mincoord) & (positions<=box1_maxcoord), axis=1)
    box1indices   = np.sum((positions>=box1_mincoord) & (positions<=box1_maxcoord), axis=1) == 3
    # print(box1indices)
    box1positions = positions[box1indices,:]
    box1data      = data[box1indices]
    # print("box1: ", box1positions)
    # print("box 1 data: ", box1data)
    box1mean, box1N, remesh_nb_points = get_mean_remesh(box1positions, box1data, depth+1, max_depth, threshold, box1_mincoord, box1_maxcoord, remesh_points, remesh_nb_points)

    # box2indices   = np.all((positions>=box2_mincoord) & (positions<=box2_maxcoord), axis=1)/
    box2indices   = np.sum((positions>=box2_mincoord) & (positions<=box2_maxcoord), axis=1) == 3
    # print("2min: ", box2_mincoord)
    # print("2max: ", box2_maxcoord)
    box2positions = positions[box2indices,:]
    box2data      = data[box2indices]
    # print(box2indices)
    # print("box2: ", box2positions)
    # print("box 2 data: ", box2data)
    box2mean, box2N, remesh_nb_points = get_mean_remesh(box2positions, box2data, depth+1, max_depth, threshold, box2_mincoord, box2_maxcoord, remesh_points, remesh_nb_points)

    # box3indices   = np.all((positions>=box3_mincoord) & (positions<=box3_maxcoord), axis=1)
    box3indices   = np.sum((positions>=box3_mincoord) & (positions<=box3_maxcoord), axis=1) == 3
    box3positions = positions[box3indices,:]
    box3data      = data[box3indices]
    # print("box3: ", box3positions)
    # print("box 3 data: ", box3data)
    box3mean, box3N, remesh_nb_points = get_mean_remesh(box3positions, box3data, depth+1, max_depth, threshold, box3_mincoord, box3_maxcoord, remesh_points, remesh_nb_points)

    # box4indices   = np.all((positions>=box4_mincoord) & (positions<=box4_maxcoord), axis=1)
    box4indices   = np.sum((positions>=box4_mincoord) & (positions<=box4_maxcoord), axis=1) == 3
    box4positions = positions[box4indices,:]
    box4data      = data[box4indices]
    # print("box4: ", box4positions)
    # print("box 4 data: ", box4data)
    box4mean, box4N, remesh_nb_points = get_mean_remesh(box4positions, box4data, depth+1, max_depth, threshold, box4_mincoord, box4_maxcoord, remesh_points, remesh_nb_points)

    # box5indices   = np.all((positions>=box5_mincoord) & (positions<=box5_maxcoord), axis=1)
    box5indices   = np.sum((positions>=box5_mincoord) & (positions<=box5_maxcoord), axis=1) == 3
    box5positions = positions[box5indices,:]
    box5data      = data[box5indices]
    # print("box5: ", box5positions)
    # print("box 5 data: ", box5data)
    box5mean, box5N, remesh_nb_points = get_mean_remesh(box5positions, box5data, depth+1, max_depth, threshold, box5_mincoord, box5_maxcoord, remesh_points, remesh_nb_points)

    # box6indices   = np.all((positions>=box6_mincoord) & (positions<=box6_maxcoord), axis=1)
    box6indices   = np.sum((positions>=box6_mincoord) & (positions<=box6_maxcoord), axis=1) == 3
    box6positions = positions[box6indices,:]
    box6data      = data[box6indices]
    # print("box6: ", box6positions)
    # print("box 6 data: ", box6data)
    box6mean, box6N ,remesh_nb_points = get_mean_remesh(box6positions, box6data, depth+1, max_depth, threshold, box6_mincoord, box6_maxcoord, remesh_points, remesh_nb_points)

    # box7indices   = np.all((positions>=box7_mincoord) & (positions<=box7_maxcoord), axis=1)
    box7indices   = np.sum((positions>=box7_mincoord) & (positions<=box7_maxcoord), axis=1) == 3
    box7positions = positions[box7indices,:]
    box7data      = data[box7indices]
    # print("box7: ", box7positions)
    # print("box 7 data: ", box7data)
    box7mean, box7N, remesh_nb_points = get_mean_remesh(box7positions, box7data, depth+1, max_depth, threshold, box7_mincoord, box7_maxcoord, remesh_points, remesh_nb_points)

    # box8indices   = np.all((positions>=box8_mincoord) & (positions<=box8_maxcoord), axis=1)
    box8indices   = np.sum((positions>=box8_mincoord) & (positions<=box8_maxcoord), axis=1) == 3
    box8positions = positions[box8indices,:]
    box8data      = data[box8indices]
    # print("box8: ", box8positions)
    # print("box 8 data: ", box8data)
    box8mean, box8N, remesh_nb_points = get_mean_remesh(box8positions, box8data, depth+1, max_depth, threshold, box8_mincoord, box8_maxcoord, remesh_points, remesh_nb_points)

    # print(data)
    totN = box1N+box2N+box3N+box4N+box5N+box6N+box7N+box8N
    # print("totN: ", totN)
    mean = (box1mean*box1N + box2mean*box2N + box3mean*box3N + box4mean*box4N + box5mean*box5N + box6mean*box6N + box7mean*box7N + box8mean*box8N)/totN
    # sum = (box1mean*box1N + box2mean*box2N + box3mean*box3N + box4mean*box4N + box5mean*box5N + box6mean*box6N + box7mean*box7N + box8mean*box8N)
    # print("sum: ", sum)
    # print("mean: ", mean)
    #Note: variance might be the wrong word: this is actually the variance if we replace the values in each box by their average
    variance = (box1mean**2*box1N+box2mean**2*box2N+box3mean**2*box3N+box4mean**2*box4N+box5mean**2*box5N+box6mean**2*box6N+box7mean**2*box7N+box8mean**2*box8N)/totN - mean**2
    stddev = np.sqrt(variance)#err, this is not correct in the slightest; give variance another name; units stddev are correct though
    # print("stddev: ", stddev)

    #add boxes to remeshed list if variance is deemed important enough/for now, I just add the current point
    if (stddev/mean > threshold):
        # remesh_points[remesh_nb_points, :] = def_maxcoord#coordinate of middle point
        # remesh_nb_points+=1
        remesh_points[remesh_nb_points, :] = (box1_mincoord + box1_maxcoord) / 2.0
        remesh_points[remesh_nb_points+1, :] = (box2_mincoord + box2_maxcoord) / 2.0
        remesh_points[remesh_nb_points+2, :] = (box3_mincoord + box3_maxcoord) / 2.0
        remesh_points[remesh_nb_points+3, :] = (box4_mincoord + box4_maxcoord) / 2.0
        remesh_points[remesh_nb_points+4, :] = (box5_mincoord + box5_maxcoord) / 2.0
        remesh_points[remesh_nb_points+5, :] = (box6_mincoord + box6_maxcoord) / 2.0
        remesh_points[remesh_nb_points+6, :] = (box7_mincoord + box7_maxcoord) / 2.0
        remesh_points[remesh_nb_points+7, :] = (box8_mincoord + box8_maxcoord) / 2.0
        remesh_nb_points+=8
        # print("adding points to array")


    # this is an attempt at using decently defined quantities, instead of some undefined variance; this does not seem to work
    # if (np.abs(box1mean-mean)/mean > threshold):
    #     remesh_points[remesh_nb_points, :] = (box1_mincoord + box1_maxcoord) / 2.0
    #     remesh_nb_points+=1
    #
    # if (np.abs(box2mean-mean)/mean > threshold):
    #     remesh_points[remesh_nb_points, :] = (box2_mincoord + box2_maxcoord) / 2.0
    #     remesh_nb_points+=1
    #
    # if (np.abs(box3mean-mean)/mean > threshold):
    #     remesh_points[remesh_nb_points, :] = (box3_mincoord + box3_maxcoord) / 2.0
    #     remesh_nb_points+=1
    #
    # if (np.abs(box4mean-mean)/mean > threshold):
    #     remesh_points[remesh_nb_points, :] = (box4_mincoord + box4_maxcoord) / 2.0
    #     remesh_nb_points+=1
    #
    # if (np.abs(box5mean-mean)/mean > threshold):
    #     remesh_points[remesh_nb_points, :] = (box5_mincoord + box5_maxcoord) / 2.0
    #     remesh_nb_points+=1
    #
    # if (np.abs(box6mean-mean)/mean > threshold):
    #     remesh_points[remesh_nb_points, :] = (box6_mincoord + box6_maxcoord) / 2.0
    #     remesh_nb_points+=1
    #
    # if (np.abs(box7mean-mean)/mean > threshold):
    #     remesh_points[remesh_nb_points, :] = (box7_mincoord + box7_maxcoord) / 2.0
    #     remesh_nb_points+=1
    #
    # if (np.abs(box8mean-mean)/mean > threshold):
    #     remesh_points[remesh_nb_points, :] = (box8_mincoord + box8_maxcoord) / 2.0
    #     remesh_nb_points+=1


    return (mean, totN, remesh_nb_points)




@numba.njit(cache=True)
def get_mean_remesh_ver2(positions, data, depth, max_depth, threshold, min_coord, max_coord, remesh_points, remesh_nb_points):
    '''
    Uses recursion (uses all data to determine whether to recurse on a smaller scale)
    '''

    if len(data)<=1 or depth==max_depth:
        #add this point to the list
        remesh_points[remesh_nb_points, :] = (min_coord + max_coord) / 2.0
        remesh_nb_points+=1
        return remesh_nb_points

    minval, maxval, meanval = np.min(data), np.max(data), np.mean(data)
    # print(maxval, minval, meanval)

    if (maxval-minval)<threshold*meanval:
        #add this point to the list
        remesh_points[remesh_nb_points, :] = (min_coord + max_coord) / 2.0
        remesh_nb_points+=1
        return remesh_nb_points

    else:
        #go and do some more recursive investigation
        delta_coord = (max_coord - min_coord) / 2.0
        #defining coordinates for sub-boxes
        def_mincoord = min_coord
        def_maxcoord = min_coord + delta_coord
        #TODO: use max by starting from true maximum (otherwise rounding errors!)

        box1_mincoord = def_mincoord + np.array([0,0,0]) * delta_coord
        box1_maxcoord = def_maxcoord + np.array([0,0,0]) * delta_coord
        box2_mincoord = def_mincoord + np.array([0,0,1]) * delta_coord
        box2_maxcoord = def_maxcoord + np.array([0,0,1]) * delta_coord
        box3_mincoord = def_mincoord + np.array([0,1,0]) * delta_coord
        box3_maxcoord = def_maxcoord + np.array([0,1,0]) * delta_coord
        box4_mincoord = def_mincoord + np.array([0,1,1]) * delta_coord
        box4_maxcoord = def_maxcoord + np.array([0,1,1]) * delta_coord
        box5_mincoord = def_mincoord + np.array([1,0,0]) * delta_coord
        box5_maxcoord = def_maxcoord + np.array([1,0,0]) * delta_coord
        box6_mincoord = def_mincoord + np.array([1,0,1]) * delta_coord
        box6_maxcoord = def_maxcoord + np.array([1,0,1]) * delta_coord
        box7_mincoord = def_mincoord + np.array([1,1,0]) * delta_coord
        box7_maxcoord = def_maxcoord + np.array([1,1,0]) * delta_coord
        box8_mincoord = def_mincoord + np.array([1,1,1]) * delta_coord
        box8_maxcoord = def_maxcoord + np.array([1,1,1]) * delta_coord
        # print("here")

        # box1indices   = np.all((positions>=box1_mincoord) & (positions<=box1_maxcoord), axis=1)
        box1indices   = np.sum((positions>=box1_mincoord) & (positions<=box1_maxcoord), axis=1) == 3
        # print(box1indices)
        box1positions = positions[box1indices,:]
        box1data      = data[box1indices]
        # print("box1: ", box1positions)
        # print("box 1 data: ", box1data)
        remesh_nb_points = get_mean_remesh_ver2(box1positions, box1data, depth+1, max_depth, threshold, box1_mincoord, box1_maxcoord, remesh_points, remesh_nb_points)

        # box2indices   = np.all((positions>=box2_mincoord) & (positions<=box2_maxcoord), axis=1)/
        box2indices   = np.sum((positions>=box2_mincoord) & (positions<=box2_maxcoord), axis=1) == 3
        # print("2min: ", box2_mincoord)
        # print("2max: ", box2_maxcoord)
        box2positions = positions[box2indices,:]
        box2data      = data[box2indices]
        # print(box2indices)
        # print("box2: ", box2positions)
        # print("box 2 data: ", box2data)
        remesh_nb_points = get_mean_remesh_ver2(box2positions, box2data, depth+1, max_depth, threshold, box2_mincoord, box2_maxcoord, remesh_points, remesh_nb_points)

        # box3indices   = np.all((positions>=box3_mincoord) & (positions<=box3_maxcoord), axis=1)
        box3indices   = np.sum((positions>=box3_mincoord) & (positions<=box3_maxcoord), axis=1) == 3
        box3positions = positions[box3indices,:]
        box3data      = data[box3indices]
        # print("box3: ", box3positions)
        # print("box 3 data: ", box3data)
        remesh_nb_points = get_mean_remesh_ver2(box3positions, box3data, depth+1, max_depth, threshold, box3_mincoord, box3_maxcoord, remesh_points, remesh_nb_points)

        # box4indices   = np.all((positions>=box4_mincoord) & (positions<=box4_maxcoord), axis=1)
        box4indices   = np.sum((positions>=box4_mincoord) & (positions<=box4_maxcoord), axis=1) == 3
        box4positions = positions[box4indices,:]
        box4data      = data[box4indices]
        # print("box4: ", box4positions)
        # print("box 4 data: ", box4data)
        remesh_nb_points = get_mean_remesh_ver2(box4positions, box4data, depth+1, max_depth, threshold, box4_mincoord, box4_maxcoord, remesh_points, remesh_nb_points)

        # box5indices   = np.all((positions>=box5_mincoord) & (positions<=box5_maxcoord), axis=1)
        box5indices   = np.sum((positions>=box5_mincoord) & (positions<=box5_maxcoord), axis=1) == 3
        box5positions = positions[box5indices,:]
        box5data      = data[box5indices]
        # print("box5: ", box5positions)
        # print("box 5 data: ", box5data)
        remesh_nb_points = get_mean_remesh_ver2(box5positions, box5data, depth+1, max_depth, threshold, box5_mincoord, box5_maxcoord, remesh_points, remesh_nb_points)

        # box6indices   = np.all((positions>=box6_mincoord) & (positions<=box6_maxcoord), axis=1)
        box6indices   = np.sum((positions>=box6_mincoord) & (positions<=box6_maxcoord), axis=1) == 3
        box6positions = positions[box6indices,:]
        box6data      = data[box6indices]
        # print("box6: ", box6positions)
        # print("box 6 data: ", box6data)
        remesh_nb_points = get_mean_remesh_ver2(box6positions, box6data, depth+1, max_depth, threshold, box6_mincoord, box6_maxcoord, remesh_points, remesh_nb_points)

        # box7indices   = np.all((positions>=box7_mincoord) & (positions<=box7_maxcoord), axis=1)
        box7indices   = np.sum((positions>=box7_mincoord) & (positions<=box7_maxcoord), axis=1) == 3
        box7positions = positions[box7indices,:]
        box7data      = data[box7indices]
        # print("box7: ", box7positions)
        # print("box 7 data: ", box7data)
        remesh_nb_points = get_mean_remesh_ver2(box7positions, box7data, depth+1, max_depth, threshold, box7_mincoord, box7_maxcoord, remesh_points, remesh_nb_points)

        # box8indices   = np.all((positions>=box8_mincoord) & (positions<=box8_maxcoord), axis=1)
        box8indices   = np.sum((positions>=box8_mincoord) & (positions<=box8_maxcoord), axis=1) == 3
        box8positions = positions[box8indices,:]
        box8data      = data[box8indices]
        # print("box8: ", box8positions)
        # print("box 8 data: ", box8data)
        remesh_nb_points = get_mean_remesh_ver2(box8positions, box8data, depth+1, max_depth, threshold, box8_mincoord, box8_maxcoord, remesh_points, remesh_nb_points)


    return remesh_nb_points

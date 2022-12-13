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
        boundary = []
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

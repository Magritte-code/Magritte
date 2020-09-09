#include <iostream>
using std::cout;
using std::endl;


struct array_gpu
{
    size_t  size;
    double* data;


    void* operator new (size_t len)
    {
        void* ptr;
        cudaMallocManaged (&ptr, len);
        cudaDeviceSynchronize ();
        return ptr;
    }

    void operator delete (void *ptr)
    {
        cudaDeviceSynchronize ();
        cudaFree (ptr);
    }

    array_gpu (const size_t s)
    {
//        array_gpu* this_ptr;
//        cudaMalloc (&this_ptr, sizeof(array_gpu));
//        this = this_ptr;
        size = s;
        cudaMallocManaged (&data, size*sizeof(double));
        cudaDeviceSynchronize ();
    }

    ~array_gpu()
    {
        cudaDeviceSynchronize ();
        cudaFree (data);
    }

};

struct array_cpu
{
    size_t  size;
    double* data;

    array_cpu (const size_t s)
    {
        size = s;
        data = (double*) std::malloc (size*sizeof(double));
    }

    ~array_cpu ()
    {
        std::free (data);
    }
};

void copy (array_cpu& arr_cpu, array_gpu& arr_gpu)
{
    cudaMemcpy(arr_gpu.data,
               arr_cpu.data,
               arr_cpu.size*sizeof(double),
               cudaMemcpyHostToDevice      );
    cudaDeviceSynchronize ();
}

void copy (array_gpu& arr_gpu, array_cpu& arr_cpu)
{
    cudaDeviceSynchronize ();
    cudaMemcpy(arr_cpu.data,
               arr_gpu.data,
               arr_gpu.size*sizeof(double),
               cudaMemcpyDeviceToHost      );
}


__global__ void kernel (array_gpu& arr_gpu)
{
    for (size_t i = 0; i < arr_gpu.size; i++)
    {
        arr_gpu.data[i]++;
    }
}




//int main ()
//{
//    const size_t size = 5;
//    double* arr = (double*) malloc(size*sizeof(double));
//
//    array a = array(size);
//
//    for (size_t i = 0; i < size; i++)
//    {
//        arr[i] = i + 4.0;
//        a.data[i] = arr[i];
//        cout << arr[i] << endl;
//    }
//
//    double* arr_dev;
//    cudaMalloc(&arr_dev, size*sizeof(double));
//
//
//    cudaMemcpy(arr_dev, arr, size*sizeof(double), cudaMemcpyHostToDevice);
//
//
//    kernel<<<1,1>>>(arr_dev);
//    cudaDeviceSynchronize();
//
//    cudaMemcpy(arr, arr_dev, size*sizeof(double), cudaMemcpyDeviceToHost);
//
//    for (size_t i = 0; i < size; i++)
//    {
//       cout << arr[i] << endl;
//    }
//
//
//    free    (arr);
//    cudaFree(arr_dev);
//
//    return (0);
//}


int main ()
{
    const size_t size = 5;

    array_cpu* arr_cpu_ptr = new array_cpu (size);
    array_gpu* arr_gpu_ptr = new array_gpu (size);

    array_cpu arr_cpu = *arr_cpu_ptr;
    array_gpu arr_gpu = *arr_gpu_ptr;

    for (size_t i = 0; i < size; i++)
    {
        arr_cpu.data[i] = i + 4.0;
        cout << arr_cpu.data[i] << endl;
    }

    copy (arr_cpu, arr_gpu);

    cout << "on gpu ---" << endl;
    cout << &arr_gpu     << endl;
    cout <<  arr_gpu_ptr << endl;

    cout << "on cpu ---" << endl;
    cout << &arr_cpu     << endl;
    cout <<  arr_cpu_ptr << endl;

    kernel<<<1,1>>>(*arr_gpu_ptr);
    cudaDeviceSynchronize();

    copy (arr_gpu, arr_cpu);

    for (size_t i = 0; i < size; i++)
    {
        cout << "--- " << arr_cpu.data[i] << endl;
    }

    return (0);
}

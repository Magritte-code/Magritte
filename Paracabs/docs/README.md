#Questions

- Currently the datatypes have to be explicitly constructed using their cusom `new` operators. As a result, we are constantly working with pointers instead of actual objects. Is there a way to circumvent this?


Array Design strategy:

Frederik: separate objects on host and device.

hemi: single object containing both a host and device pointer.


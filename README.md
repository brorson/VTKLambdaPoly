May 2021.

This is an attempt to create a viewer for the lambda poly.  The
desire is to use a vtkImageData as the place to write the
values of the lambda poly, then display them.  Then use VTK
facilities to allow the user to rubber band an area and zoom in.

It turned out to be difficult -- I got stuck in the VTK morass.

--------------------------------------------------------------

April 2025.

I poked at this program and after some work it now
works fine.  I also implemented a feature where the program shows
where the lambda poly's roots are.  The roots are colored black.  

To build:  

cmake -S . -B build  
cd build  
make  

To run:  

cd build  
./LambdaPoly -N 16



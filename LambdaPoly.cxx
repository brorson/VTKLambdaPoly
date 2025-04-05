#include <vtkNamedColors.h>
#include <vtkNew.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkInteractorStyleImage.h>
#include <vtkRenderer.h>
#include <vtkUniformGrid.h>
#include <vtkLookupTable.h>
#include <vtkImageMapToColors.h>
#include <vtkImageActor.h>
#include <vtkImageMapper3D.h>
#include <vtkScalarBarActor.h>
#include <vtkTextProperty.h>
#include <vtkVector.h>
#include <vtkMath.h>
#include <vtkNamedColors.h>

#include <stdlib.h>
#include <unistd.h>
#include <getopt.h>

#include <complex>
#include <vector>
#include <array>
// #include <mpreal.h> 

// Window dimensions
#define NX 600
#define NY 600

// Default number of logistic map iterations.
#define NITER 64

// Amount to grow/shrink when turning mouse wheel.
#define SCALE 1.6

// Value at which to saturate poly (both pos and neg)
#define SAT 1.0

// Type declaration macro
#define MY_CREATE(type, name) \
    type *name = type::New()

//-----------------------------------------------------------------
// Declare fcns defining the Lambda Poly surface in the complex plane.
// Under CUDA these will become __global__ fcns.
void computeLambdaPoly(vtkUniformGrid *imageData, int N);
std::complex<double> f(double x, double y, int N);

// host-side graphics manipulation fcns.
void moveZoom(int i, int j, double zoom);
void moveTranslate(vtkVector<int, 4> p);

//-------------------------------------------------------------
// Create most VTK objects as globals so I can access them from
// everywhere.  Some say it's bad practice to use globals, but
// I say it's easier than trying to pass around pointers to
// objects from my main prog to the callbacks in the custom
// interactor
MY_CREATE(vtkUniformGrid, rImageData);
//MY_CREATE(vtkUniformGrid, iImageData);
MY_CREATE(vtkImageMapToColors, colorComplexPlane);
MY_CREATE(vtkImageActor, imageActor);
MY_CREATE(vtkRenderer, renderer);
MY_CREATE(vtkRenderWindow, renWin);
MY_CREATE(vtkRenderWindowInteractor, iren);
MY_CREATE(vtkLookupTable, lookupTable);
MY_CREATE(vtkScalarBarActor, scalarBar);


// Scalar global vars
int N;          // Number of logistic map iterations (settable).
double w, h;    // Width, height of image in real numbers.
                // I can probably replace this global with
                // quantities derived from the image itself.

//------------------------------------------------------------------
// Started from
// https://kitware.github.io/vtk-examples/site/Cxx/Interaction/MouseEvents/
class customMouseInteractorStyle : public vtkInteractorStyleImage
{
public:
  static customMouseInteractorStyle* New();
  vtkTypeMacro(customMouseInteractorStyle, vtkInteractorStyleImage);

  vtkVector<int, 4> evt;    // Event coords -- down xy, up xy
  double scale;             // Scale to zoom in/out
  bool quit;
  
  void OnMouseWheelForward() override {
    //std::cout << "MouseWheelForward ... ";
    int i = this->Interactor->GetEventPosition()[0];
    int j = this->Interactor->GetEventPosition()[1];
    //std::cout << "[i,j] = [" << i << ", " << j << "]" << std:: endl;
    scale = 1.0/SCALE;
    moveZoom(i, j, scale);
    // Tell pipeline to update
    renderer->ResetCamera();
    renWin->Render();
  }

  void OnMouseWheelBackward() override {
    //std::cout << "MouseWheelBackward ... ";
    int i = this->Interactor->GetEventPosition()[0];
    int j = this->Interactor->GetEventPosition()[1];
    //std::cout << "[i,j] = [" << i << ", " << j << "]" << std:: endl;
    scale = SCALE;    
    moveZoom(i, j, scale);
    // Tell pipeline to update
    renderer->ResetCamera();
    renWin->Render();
  }
  
  void OnMiddleButtonDown() override
  {
    // This returns point in window where button went down.
    //std::cout << " MiddleButtonDown ..." << std::endl;
    int i = this->Interactor->GetEventPosition()[0];
    int j = this->Interactor->GetEventPosition()[1];
    //std::cout << "[i,j] = [" << i << ", " << j << "]" << std:: endl;
    evt[0] = i;
    evt[1] = j;
    vtkInteractorStyleImage::OnMiddleButtonDown();
    // Nothing to do here -- must wait until button pops up.
  }


  void OnMiddleButtonUp() override
  {
    //std::cout << " MiddleButtonUp ..." << std::endl;
    int i = this->Interactor->GetEventPosition()[0];
    int j = this->Interactor->GetEventPosition()[1];
    //std::cout << "[i,j] = [" << i << ", " << j << "]" << std:: endl;
    evt[2] = i;
    evt[3] = j;
    vtkInteractorStyleImage::OnMiddleButtonUp();
    moveTranslate(evt);
    // Tell pipeline to update
    renderer->ResetCamera();
    renWin->Render();
  }

  
  void OnLeftButtonDown() override {
    std::cout << "Left button down ..." << std::endl;
    //vtkInteractorStyleImage::OnLeftButtonDown();
  }

  void OnLeftButtonUp() override {
    std::cout << "Left button up ..." << std::endl;
    //vtkInteractorStyleImage::OnLeftButtonUp();    
  } 

  void OnKeyDown() override {
    //std::cout << "Key down ..." << std::endl;
    std::string key = this->Interactor->GetKeySym();
    std::cout << "Key pressed: " << key << std::endl;
    if (key == "q"){
      quit = true;    
    } else {
      quit = false;    
    }
    this->Interactor->ExitCallback ();
  }

  
  vtkVector<int, 4> getEvt(void) {
    return evt;
  }

  double getScale(void) {
    return scale;
  }

};
vtkStandardNewMacro(customMouseInteractorStyle);
// Instantiate iStyle here, after defining it.
MY_CREATE(customMouseInteractorStyle, iStyle);

//------------------------------------------------------
// Graphics manipulation
void moveZoom(int i, int j, double zoom) {
  // Convert first pixel point to real number 
  double xyz[3];
  int iz = 0;
  double x0, y0, w, h;
  double dx, dy;
  double xmin, xmax, ymin, ymax;
  int dims[2];
  double *zz;  // Z min and max

  /*
  std::cout << "moveZoom, before update, rImageData = " << std::endl;
  rImageData->PrintSelf(std::cout,vtkIndent(2));
  std::cout << "moveZoom, before update, imageActor = " << std::endl;
  imageActor->PrintSelf(std::cout,vtkIndent(2));
  std::cout << "moveZoom, before update, renWin = " << std::endl;
  renWin->PrintSelf(std::cout,vtkIndent(2));
  std::cout << "----------------------------------------------" << std::endl;
  */
  
  // Get old height and width
  rImageData->GetSpacing(xyz);
  w = NX*xyz[0];
  h = NY*xyz[1];

  // Get old image origin
  rImageData->GetOrigin(xyz);
  xmin = xyz[0];
  ymin = xyz[1];
  
  rImageData->TransformIndexToPhysicalPoint (i, j, iz, xyz);
  x0 = xyz[0];  // New center of image
  y0 = xyz[1];  // New center of image

  // Update w and h
  w = w*zoom;
  h = h*zoom;

  // Convert these values to new min, max, and spacing
  xmin = x0 - w/2.0;
  xmax = x0 + w/2.0;  
  ymin = y0 - h/2.0;
  ymax = y0 + h/2.0;  
  dx = w/(NX-1);
  dy = h/(NY-1);
  
  // Now that I have new extents, must update ImageData
  rImageData->SetSpacing(dx, dy, 1.0);
  rImageData->SetOrigin(xmin, ymin, 0.0);  // This sets lower left corner.
  rImageData->AllocateScalars(VTK_DOUBLE, 1);

  // Get new origin, height and width as check
  rImageData->GetOrigin(xyz);
  xmin = xyz[0];
  ymin = xyz[1];
  rImageData->GetSpacing(xyz);
  w = NX*xyz[0];
  h = NY*xyz[1];
  printf("moveZoom: New x0 = %f, y0 = %f, w = %e, h = %e\n", x0, y0, w, h);

  /*
  std::cout << "moveZoom, after update, rImageData = " << std::endl;
  rImageData->PrintSelf(std::cout,vtkIndent(2));
  std::cout << "moveZoom, after update, imageActor = " << std::endl;
  imageActor->PrintSelf(std::cout,vtkIndent(2));
  std::cout << "moveZoom, after update, renWin = " << std::endl;
  renWin->PrintSelf(std::cout,vtkIndent(2));
  printf("New x0 = %f, y0 = %f, w = %e, h = %e\n", x0, y0, w, h);
  std::cout << "----------------------------------------------" << std::endl;
  */
  
  // Now compute lambda poly using new origin and spacing.
  computeLambdaPoly(rImageData, N);  // Compute the whole set.

  // Find min & max of data in rImageData
  /*
  double min = 1;
  double max = 0;
  rImageData->GetDimensions(dims);
  int Nx = dims[0];
  int Ny = dims[1];
  double *z;
  int iy = 0;
  // Set color bar by finding min & max on real number line.
  //for (int iy = 0; iy < Ny; iy++) {
    for (int ix = 0; ix < Nx; ix++) {
      // First figure out what is my [x,y] point in the complex plane
      z = static_cast<double*>(rImageData->GetScalarPointer(ix, iy, iz));      
      // Set color map so limits of color are set by values on visible
      // real line.
      if (*z > 1.0) break;    // Out of bounds.
      else if (*z < min) min = *z;
      else if (*z > max) max = *z;
    }
  //}
  //printf("min = %f, max = %f\n", min, max);
  // Choose larger of distance between 1/2 and min or max.
  double diff;
  if (abs(max-0.5) > abs(0.5-min)) diff = max-0.5;
  else diff = 0.5-min;
  zz = lookupTable->GetTableRange();
  // Only shrink color range, don't grow it.
  if (*zz < 0.5-diff) lookupTable->SetTableRange(0.5-diff, 0.5+diff);
  lookupTable->Build();
  scalarBar->SetLookupTable( lookupTable );
  imageActor->GetMapper()->Update();
  */
  
  return;
}

//------------------------------------------------------
void moveTranslate(vtkVector<int, 4> p) {
  // Convert first pixel point to real number.  p
  // is 4 element vector holding [xdn, ydn, xup, yup].
  int iz = 0;
  double xyz[3];
  double xmin, ymin, x1, y1, x2, y2;
  double dx, dy;
  double w, h;

  // Location of middle button down
  rImageData->TransformIndexToPhysicalPoint (p[0], p[1], iz, xyz);
  x1 = xyz[0];
  y1 = xyz[1];

  // Location of middle button up
  rImageData->TransformIndexToPhysicalPoint (p[2], p[3], iz, xyz);
  x2 = xyz[0];
  y2 = xyz[1];

  printf("moveTranslate, [x1, y1] = [%f, %f] ... [x2, y2] = [%f, %f]\n",
	 x1, y1, x2, y2);

  // Get old image origin
  rImageData->GetOrigin(xyz);
  xmin = xyz[0];
  ymin = xyz[1];

  // Get old height and width
  rImageData->GetSpacing(xyz);
  w = NX*xyz[0];
  h = NY*xyz[1];

  // Amount to translate.
  dx = x2-x1;
  dy = y2-y1;

  // New origin.
  xmin = xmin - 1.3*dx;
  ymin = ymin - 1.3*dy;
  
  // Move origin to new location
  rImageData->SetOrigin(xmin, ymin, 0.0);  // This sets lower left corner.

  // Get new origin, height and width as check
  rImageData->GetOrigin(xyz);
  xmin = xyz[0];
  ymin = xyz[1];
  rImageData->GetSpacing(xyz);
  w = NX*xyz[0];
  h = NY*xyz[1];
  printf("New x0 = %f, y0 = %f, w = %e, h = %e\n",
	 xmin+w/2.0, ymin+h/2.0, w, h);
  
  // Now compute lambda poly using new dimensions
  computeLambdaPoly(rImageData, N);  // Compute the whole set.  
  return;
}


//======================================================================
int main(int argc, char* argv[])
{
  vtkNew<vtkNamedColors> colors;
  double alpha;
  int i, j;      // Center of image in pixel coords.
  int Nx = NX;   // Width of image in pixels
  int Ny = NY;   // Height of image in pixels
  double x0, y0;
  double xmin, xmax, ymin, ymax;
  int c;

  /*
  printf("argc = %d\n", argc);
  for(i = 0 ; i < argc ; ++i) {
    printf("argv[%d] = '%s'\n", i, argv[i]);
  }
  */
  
  // This is initial view window
  x0 = 1.0;
  y0 = 0.0;
  w = 6.0;
  h = 6.0;
  N = NITER;
  // Process command line args (if any)
  static struct option long_options[] =
    {
     {"x",  required_argument, 0, 'x'},
     {"y",  required_argument, 0, 'y'},
     {"w",  required_argument, 0, 'w'},
     {"h",  required_argument, 0, 'h'},
     {"N",  required_argument, 0, 'N'},       
     {0, 0, 0, 0}
    };
  /* getopt_long stores the option index here. */
  int option_index = 0;
  while (1) {
    c = getopt_long (argc, argv, "x:y:w:h:N:", long_options, &option_index);
    //std::cout << "c = " << c << std::endl;
    if (c == -1) break;
    switch (c) {
    case 'x':
      x0 = atof(optarg);
      //std::cout << "x = " << x0 << std::endl;
      break;
    case 'y':
      y0 = atof(optarg);
      //std::cout << "y = " << y0 << std::endl;      
      break;
    case 'w':
      w = atof(optarg);
      //std::cout << "w = " << w << std::endl;      
      break;
    case 'h':
      h = atof(optarg);
      //std::cout << "h = " << h << std::endl;            
      break;
    case 'N':
      N = atoi(optarg);
      //std::cout << "N = " << N << std::endl;            
      break;
    case '?':
      fprintf (stderr,
	       "Unknown option character 0x%x'.\n",
	       optopt);
      return 1;
    default:
      abort ();
    }
  }
  printf("Starting x0 = %f, y0 = %f, w = %e, h = %e\n", x0, y0, w, h);  
  
  //--------------------------------------------------------------------

  // Map the scalar values in the image to colors with a lookup table
  lookupTable->SetNumberOfTableValues(256);
  //lookupTable->SetTableRange(-SAT, SAT);
  //lookupTable->SetTableRange(0, N);
  lookupTable->SetTableRange(0, 1);  
  lookupTable->SetAboveRangeColor(0.0, 0.0, 0.0, 0.0);
  lookupTable->SetNanColor(0.0, 0.0, 0.0, 1.0);
  //lookupTable->SetRampToLinear();
  //lookupTable->SetRampToSQRT();
  lookupTable->SetRampToSCurve();
  lookupTable->Build();

  //----------------------------------------------------------------
  // Colorbar to show off color map
  scalarBar->SetLookupTable( lookupTable );
  scalarBar->SetOrientationToVertical();
  scalarBar->GetLabelTextProperty()->SetColor(0,0,1);
  scalarBar->GetTitleTextProperty()->SetColor(0,0,1);
  
  // Position scalarBar in window
  scalarBar->GetPositionCoordinate()->SetCoordinateSystemToNormalizedViewport();
  scalarBar->SetPosition(0.85, 0.1);
  scalarBar->SetWidth(.10);
  scalarBar->SetHeight(0.8);


  //--------------------------------------------------------
  // Pass the original image and the lookup table to a
  // filter to create a color image.
  cout << "Configure colorComplexPlane ... " << endl;
  colorComplexPlane->SetLookupTable(lookupTable);
  colorComplexPlane->PassAlphaToOutputOn();
  colorComplexPlane->SetInputData(rImageData);  // set to real or imag plane

  // Configure ImageData
  cout << "Configure ImageData ... " << endl;
  xmin = x0 - w/2.0;
  xmax = x0 + w/2.0;  
  ymin = y0 - h/2.0;
  ymax = y0 + h/2.0;  
  double dx = (xmax-xmin)/(Nx-1);
  double dy = (ymax-ymin)/(Ny-1);
  printf("xmin = %f, xmax = %f, ymin = %f, ymax = %f, dx = %f, dy = %f\n",
	 xmin, xmax, ymin, ymax, dx, dy);
  
  rImageData->SetExtent( 0, Nx-1, 0, Ny-1, 0, 0 );  // Set image size in pixels
  rImageData->SetSpacing(dx, dy, 1.0);
  rImageData->SetOrigin(xmin, ymin, 0.0);   // This sets pos of left corner.
  rImageData->AllocateScalars(VTK_DOUBLE, 1); 
  
  /*
  iImageData->SetExtent( 0, Nx-1, 0, Ny-1, 0, 0 );  // Set image size in pixels
  iImageData->SetSpacing(dx, dy, 1.0);
  iImageData->SetOrigin(xmin, ymin, 0.0);   // This sets pos of left corner.
  iImageData->AllocateScalars(VTK_DOUBLE, 1); 
  */

  // Compute initial lambda poly for display
  cout << "Compute initial lambda poly ... " << endl;  
  computeLambdaPoly(rImageData, N);  // Compute the whole set.
  
  // Configure image actor.  Actor has built-in mapper.
  cout << "Configure image actor ... " << endl;    
  imageActor->InterpolateOff();
  imageActor->GetMapper()->SetInputConnection(colorComplexPlane->GetOutputPort());
  
  // Configure renderer
  cout << "Configure renderer ..." << endl;
  renderer->AddActor(imageActor);
  renderer->AddActor(scalarBar);
  renderer->SetBackground(colors->GetColor3d("MidnightBlue").GetData());

  // Configure render window
  cout << "Configure render window ..." << endl;
  renWin->AddRenderer(renderer);
  renWin->SetSize(Nx, Ny); // set window size in pixels
  renWin->SetWindowName("Lambda poly in complex plane");
            
  // Configure interactor and interactor style
  iren->SetRenderWindow(renWin);
  iStyle->SetInteractor(iren);
  iren->SetInteractorStyle(iStyle);
  
  // Start rendering thread
  cout << "Start rendering thread......" << endl;  
  renWin->Render();

  cout << "Initialize interactor......" << endl;  
  iren->Initialize();
  //-------------------------------------------------------------------

  std::cout << "----------------------------------------" << std::endl;

  cout << "Start interactor event loop......" << endl;
  iren->Start();

  // If I get here, it's because the event loop terminated.
  
  if (iStyle->quit == true) {
    std::cout << "User requested quit.  Exiting ..." << std::endl;
    return 0;
  } else {
    std::cout << "Returned from event loop for unknown reasons." << std::endl;
    return -1;
  }

}

//=====================================================================
void computeLambdaPoly(vtkUniformGrid *imageData, int N) {
  int dims[3];
  double xyz[3];
  double x, y;
  std::complex<double> z;
  double zx, zy, zmag;
  double* pixel;  // Temp variable
  int iz = 0;

  std::cout << "Computing lambda poly ... " << std::flush;

  imageData->GetDimensions(dims);
  int Nx = dims[0];
  int Ny = dims[1];

  for (int iy = 0; iy < Ny; iy++) {
    for (int ix = 0; ix < Nx; ix++) {
      // First figure out what is my [x,y] point in the complex plane
      imageData->TransformIndexToPhysicalPoint (ix, iy, iz, xyz);
      x = xyz[0];
      y = xyz[1];

      // Now compute complex value of fcn at this point and stick it into
      // imageData.
      z = f(x,y,N);  // Do iteration N times
      // Get ptr to place to put computed value
      pixel = static_cast<double*>(imageData->GetScalarPointer(ix, iy, iz));
      // Insert computed value.
      zx = std::real(z);
      zy = std::imag(z);
      zmag = std::sqrt(zx*zx + zy*zy);
      if (zmag < 0.05) {
	// Color roots different color
	*pixel = vtkMath::Nan();
      } else {
	//*pixel = zmag;
	*pixel = zx;
	//*pixel = zy;
      }

    }
  }
  
  std::cout << "done!" << std::endl;
  return;
}

//---------------------------------------------------------------
// This fcn iterates a point in the complex plane.
inline std::complex<double> f(double lamr, double lami, int N) {
  
  register std::complex<double> lam(lamr, lami);
  register std::complex<double> x(0.5, 0.0);

  // Do iteration
  for (int i=0; i<N; i++) {
    x = lam*x*(1.0-x);
  }

  // Saturate values
  if (x.real() > SAT) {
    x.real(SAT);
  } else if (x.real() < -SAT) {
    x.real(-SAT);
  }

  if (x.imag() > SAT) {
    x.imag(SAT);
  } else if (x.imag() < -SAT) {
    x.imag(-SAT);
  }

  std::complex<double> z = x;  
  
  return z;
}



#include <vtkImageData.h>
#include <vtkPNGWriter.h>
#include <vtkPointData.h>
#include <vtkCellArray.h>
#include <vtkPoints.h>
#include <vtkUnsignedCharArray.h>
#include <vtkFloatArray.h>
#include <vtkDataSetReader.h>
#include <vtkRectilinearGrid.h>
#include <vtkFloatArray.h>
#include <vtkPolyData.h>
#include <vtkDataSetWriter.h>
#include <vtkTubeFilter.h>
#include <vtkPolyDataNormals.h>
#include <vtkSphereSource.h>

#include <vtkCamera.h>
#include <vtkDataSetMapper.h>
#include <vtkRenderer.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkSmartPointer.h>

#include <math.h>

// ****************************************************************************
//  Function: GetNumberOfPoints
//
//  Arguments:
//     dims: an array of size 3 with the number of points in X, Y, and Z.
//           2D data sets would have Z=1
//
//  Returns:  the number of points in a rectilinear mesh
//
// ****************************************************************************

int GetNumberOfPoints(const int *dims)
{
    // 3D
    //return dims[0]*dims[1]*dims[2];
    // 2D
    return dims[0]*dims[1];
}

// ****************************************************************************
//  Function: GetNumberOfCells
//
//  Arguments:
//
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  the number of cells in a rectilinear mesh
//
// ****************************************************************************

int GetNumberOfCells(const int *dims)
{
    // 3D
    //return (dims[0]-1)*(dims[1]-1)*(dims[2]-1);
    // 2D
    return (dims[0]-1)*(dims[1]-1);
}


// ****************************************************************************
//  Function: GetPointIndex
//
//  Arguments:
//      idx:  the logical index of a point.
//              0 <= idx[0] < dims[0]
//              1 <= idx[1] < dims[1]
//              2 <= idx[2] < dims[2] (or always 0 if 2D)
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  the point index
//
// ****************************************************************************

int GetPointIndex(const int *idx, const int *dims)
{
    // 3D
    //return idx[2]*dims[0]*dims[1]+idx[1]*dims[0]+idx[0];
    // 2D
    return idx[1]*dims[0]+idx[0];
}


// ****************************************************************************
//  Function: GetCellIndex
//
//  Arguments:
//      idx:  the logical index of a cell.
//              0 <= idx[0] < dims[0]-1
//              1 <= idx[1] < dims[1]-1 
//              2 <= idx[2] < dims[2]-1 (or always 0 if 2D)
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  the cell index
//
// ****************************************************************************

int GetCellIndex(const int *idx, const int *dims)
{
    // 3D
    //return idx[2]*(dims[0]-1)*(dims[1]-1)+idx[1]*(dims[0]-1)+idx[0];
    // 2D
    return idx[1]*(dims[0]-1)+idx[0];
}

// ****************************************************************************
//  Function: GetLogicalPointIndex
//
//  Arguments:
//      idx (output):  the logical index of the point.
//              0 <= idx[0] < dims[0]
//              1 <= idx[1] < dims[1] 
//              2 <= idx[2] < dims[2] (or always 0 if 2D)
//      pointId:  a number between 0 and (GetNumberOfPoints(dims)-1).
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  None (argument idx is output)
//
// ****************************************************************************

void GetLogicalPointIndex(int *idx, int pointId, const int *dims)
{
    // 3D
    // idx[0] = pointId%dim[0];
    // idx[1] = (pointId/dims[0])%dims[1];
    // idx[2] = pointId/(dims[0]*dims[1]);

    // 2D
    idx[0] = pointId%dims[0];
    idx[1] = pointId/dims[0];
}


// ****************************************************************************
//  Function: GetLogicalCellIndex
//
//  Arguments:
//      idx (output):  the logical index of the cell index.
//              0 <= idx[0] < dims[0]-1
//              1 <= idx[1] < dims[1]-1 
//              2 <= idx[2] < dims[2]-1 (or always 0 if 2D)
//      cellId:  a number between 0 and (GetNumberOfCells(dims)-1).
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  None (argument idx is output)
//
// ****************************************************************************

void GetLogicalCellIndex(int *idx, int cellId, const int *dims)
{
    // 3D
    // idx[0] = cellId%(dims[0]-1);
    // idx[1] = (cellId/(dims[0]-1))%(dims[1]-1);
    // idx[2] = cellId/((dims[0]-1)*(dims[1]-1));

    // 2D
    idx[0] = cellId%(dims[0]-1);
    idx[1] = cellId/(dims[0]-1);
}



/***************************************************************************
// Function: getCellFromVal
//
// Arguments:
//      A: an array that contains either the X or the Y locations of a rectilinear mesh
//      dim: Length of array A
//      val: the value we're trying to find

//  Returns: the Cell index that contains a certain value
//              -1 if the value isn't in the range of the array values
**************************************************************************/
int
GetCellFromVal(const float *A, const int *dim, float val) 
{
    int returnVal = -1;
    for (int i=0; i < *(dim)-1; i++) {
        if (A[i] <= val && A[i+1] >= val){
            if (val == A[i+1]){
               returnVal = i+1;
               break;
            }
            else {
                returnVal = i;
                break;
            }
        }
    }
    return returnVal;
}

/***************************************************************************
// Function: getT
//
// Arguments:
//      X: float val to interpolate
//      A: float val to interpolate between
//      B: float val to interpolate between

//  Returns: the float value for T
**************************************************************************/


float
GetT(const float X, const float A, const float B) {
    return ((X-A)/(B-A));
}

// ****************************************************************************
//  Function: Interpolate
//
//  Arguments:
//     Fa: F(A)
//     Fb: F(B)
//     X: The value we want to interpolate
//     A, B
//
//
//   Returns: An interpolated value
//
// ****************************************************************************
void
Interpolate(const float Fa, const float Fb, const float X, const float A, const float B, float *returnVal)
{
    float t = GetT(X, A, B);
    *returnVal = (Fa + (t * (Fb-Fa)));
}

// ****************************************************************************
//  Function: EvaluateVectorFieldAtLocation
//
//  Arguments:
//     pt: a two-dimensional location
//     dims: an array of size two.  
//              The first number is the size of the array in argument X, 
//              the second the size of Y.
//     X: an array (size is specified by dims).  
//              This contains the X locations of a rectilinear mesh.
//     Y: an array (size is specified by dims).  
//              This contains the Y locations of a rectilinear mesh.
//     F: a vector field defined on the mesh.  Its size is 2*dims[0]*dims[1].
//        The first value in the field is the x-component for the first point.
//        The second value in the field is the y-component for the first point.
//
//     rv (output): the interpolated field value. (0,0) if the location is out of bounds.
//
// ****************************************************************************

void EvaluateVectorFieldAtLocation(const float *pt, const int *dims, const float *X, 
                              const float *Y, const float *F, float *rv)
{
    rv[0] = 0; // setting the x-component of the velocity
    rv[1] = 0; // setting the y-component of the velocity
    
    
    //Get the Cell that contains point P
    int a = GetCellFromVal(X, dims, pt[0]);
    int b = GetCellFromVal(Y, &dims[1], pt[1]);

    //If the point pt is within the mesh, then continue
    if (a != -1 && b != -1) {
        float scalarValues[8];
        int vertices[8];

        //Get Cell's Vertices v0 v1 v2 v3
        
        vertices[0] = a;
        vertices[1] = b;
        
        vertices[2] = a+1;         //v1[0]
        vertices[3] = b;           //v1[1]
        
        vertices[4] = a;           //v2[0]
        vertices[5] = b+1;

        vertices[6] = a+1;         //v3[0]
        vertices[7] = b+1;

        // Get x- and y- components of velocity at the four vertices of the cell
        for (int i =0; i < 7; i += 2) {
            int pointIndex = GetPointIndex((vertices+i), dims)*2;
            scalarValues[i] = F[pointIndex];
            scalarValues[i+1] = F[pointIndex+1];
        }

        //Perform linear interpolation to location P for X component
        float line1, line2;
        //Interpolate between v0 and v1 to get I1
        Interpolate(scalarValues[0], scalarValues[2], pt[0], X[vertices[0]], X[vertices[2]], &line1);
        //Interpolate between v2 and v3 to get I4
        Interpolate(scalarValues[4], scalarValues[6], pt[0], X[vertices[4]], X[vertices[6]], &line2);
        //Interpolate between v0 and v2 to get I2
        Interpolate(line1, line2, pt[1], Y[vertices[1]], Y[vertices[5]], rv);

        float line4, line5;
        //Perform linear interpolation to location P for Y component
        Interpolate(scalarValues[1], scalarValues[5], pt[1], Y[vertices[1]], Y[vertices[5]], &line4);
        //Interpolate between v2 and v3 to get I4
        Interpolate(scalarValues[3], scalarValues[7], pt[1], Y[vertices[3]], Y[vertices[7]], &line5);
        //Interpolate between v0 and v2 to get I2
        Interpolate(line4, line5, pt[0], X[vertices[0]], X[vertices[2]], (rv+1));
    }
}

// ****************************************************************************
//  Function: Get Speed
//
//  Arguments:
//     rv: two dimensional velocity at a location
//      Essentially performs pythagorean theorem
//
//  Returns:
//     returnVal: the speed at that location
//
// ****************************************************************************
float
GetSpeed(float *rv) {
    float returnVal = 0.0;
    float squareOne =0.0;
    float squareTwo = 0.0;

    squareOne = rv[0]*rv[0];
    squareTwo = rv[1]*rv[1];
    returnVal = sqrt(squareOne+squareTwo);

    return returnVal;
}

// ****************************************************************************
//  Function: AdvanceStep
//
//  Arguments:
//     pos0:  a two dimensional initial position
//     time0: an initial time  
//     h: The size of the Euler step
//     v: magnitudes of two velocity vectors
//                  first is in X, second is in Y
//     
//     pos1 (output): a two dimensional finishing position
//                      first is in X, second is in Y
//     time1 (output): a finishing time
//     
//
// ****************************************************************************
void
AdvanceStep(const float *pos0, const float h, const float *v, float *pos1)
{
    //Get new position in X
    pos1[0] = pos0[0] + (h * v[0]);
    //Get new position in Y
    pos1[1] = pos0[1] + (h * v[1]);
}

// ****************************************************************************
//  Function: AdvectWithEulerStep
//
//  Arguments:
//     pt: the seed location (two-dimensions)
//     dims: an array of size two.  
//              The first number is the size of the array in argument X, 
//              the second the size of Y.
//     X: an array (size is specified by dims).  
//              This contains the X locations of a rectilinear mesh.
//     Y: an array (size is specified by dims).  
//              This contains the Y locations of a rectilinear mesh.
//     F: a vector field defined on the mesh.  Its size is 2*dims[0]*dims[1].
//     h: The size of the Euler step
//     nsteps: The number of Euler steps to take
//     output_locations (output): An array of size 2*(nsteps+1).  It's first entry
//        should be the seed location.  The second entry should be the result
//        of the first advection step.  The final entry should be the result
//        of the final advection step.
//     speeds (output): An array of size (nsteps+1).  It's first entry should be the
//        speed at the seed location.  It's final entry should be the speed
//        at the location of the result of the final advection step.
//        Recall that speed is the magnitude of the velocity.
//
// ****************************************************************************

void
AdvectWithEulerStep(const float *pt, const int *dims, const float *X, 
                    const float *Y, const float *F, 
                    float h, int nsteps, float *output_locations, float *speeds)
{
    float currentVelocity[2];                //vOld is for storing the velocity magnitudes at a certain location Pi
    float newVelocity[2];                //vNew is for storing the velocity magnitudes at a certain location Pi+1
    float currentPosition[2];
    float newPosition[2];         //positionNew is for storing the output locations from the previous step
    int stepsTaken = 1;     
    
//Initialize the seed location
    //Get First Position
    output_locations[0] = pt[0]; // set the x component of the first output location
    output_locations[1] = pt[1]; // set the y component of the first output location

    //Get First Velocities
    EvaluateVectorFieldAtLocation(pt, dims, X, Y, F, currentVelocity);

    //Get First Speed
    speeds[0] = GetSpeed(currentVelocity); 

    while (stepsTaken <= nsteps) {

        //get the new position from the old position and velocities
        AdvanceStep(&output_locations[(stepsTaken-1)*2], h, currentVelocity, newPosition);

        //get the speed from the new velocities and store it
        speeds[stepsTaken] = GetSpeed(currentVelocity);
        
        //get the new velocities from the new position
        EvaluateVectorFieldAtLocation(newPosition, dims, X, Y, F, newVelocity);

        //store the new positions
        output_locations[stepsTaken*2] = newPosition[0];
        output_locations[stepsTaken*2+1] = newPosition[1];

        currentVelocity[0] = newVelocity[0];
        currentVelocity[1] = newVelocity[1];

        //rinse and repeat
        stepsTaken++;
    }

}

// ****************************************************************************
//  Function: CalculateArcLength
//
//  Arguments:
//     locations: an array of 2D locations.
//     nlocations: the number of locations in the array "locations".
//
//  Returns: the arc length, meaning the distance between each successive
//           pair of points
//
// ****************************************************************************

float
CalculateArcLength(const float *output_locations, int nlocations)
{
    float returnVal = 0.0;
    float differences[2];
    float arcLength = 0.0;
    differences[0] = 0.0;
    differences[1] = 0.0;

    if (nlocations > 1) {
        for (int i=0; i < nlocations-1; i++) {
            //The X Value from location i+1 - X val from location [i]
            differences[0] = output_locations[(i+1)*2] - output_locations[i*2];
            // Y value from location i+1 - Y val from location i
            differences[1] = output_locations[((i+1)*2)+1] - output_locations[(i*2)+1];
            // Find the distance between the two locations
            arcLength = sqrt( (differences[0]*differences[0] + differences[1]*differences[1]));
            // Add it to the accumulator
            returnVal += arcLength;
        }
    }
    return returnVal;
}

void
WriteImage(vtkImageData *img, const char *filename)
{
    std::string full_filename = filename;
    full_filename += ".png";
    vtkPNGWriter *writer = vtkPNGWriter::New();
    writer->SetInputData(img);
    writer->SetFileName(full_filename.c_str());
    writer->Write();
    writer->Delete();
}

vtkImageData *
NewImage(int width, int height)
{
    vtkImageData *image = vtkImageData::New();
    image->SetDimensions(width, height, 1);
    //image->SetWholeExtent(0, width-1, 0, height-1, 0, 0);
    //image->SetUpdateExtent(0, width-1, 0, height-1, 0, 0);
    //image->SetNumberOfScalarComponents(3);
    image->AllocateScalars(VTK_UNSIGNED_CHAR, 3);
    //image->AllocateScalars();

    return image;
}

// VTK files are only 3D, so the vector data is all of the form (X,Y,0).
// Remove the 0's since it is counter-intuitive for students who are 
// thinking of this as 2D data.
float *
Convert3DVectorDataTo2DVectorData(const int *dims, const float *F)
{
    float *rv = new float[dims[0]*dims[1]*2];
    int index3D = 0;
    int index2D = 0;
    for (int i = 0 ; i < dims[0] ; i++)
       for (int j = 0 ; j < dims[1] ; j++)
       {
           rv[index2D]   = F[index3D];
           rv[index2D+1] = F[index3D+1];
           index2D += 2;
           index3D += 3;
       }

    return rv;
}

vtkPolyData *
CreateVTKPolyData(int nseeds, int nsteps, float **output_locations, float **speeds)
{
    int numPoints = nseeds*(nsteps+1);
    vtkPoints *pts = vtkPoints::New();
    pts->SetNumberOfPoints(numPoints);
    vtkFloatArray *var = vtkFloatArray::New();
    var->SetName("speed");
    var->SetNumberOfTuples(numPoints);
    int ptIdx = 0;
    vtkCellArray *lines = vtkCellArray::New();
    lines->EstimateSize(numPoints,2);
    for (int i = 0 ; i < nseeds ; i++)
    {
        for (int j = 0 ; j < nsteps+1 ; j++)
        {
            double pt[3];
            pt[0] = output_locations[i][2*j];
            pt[1] = output_locations[i][2*j+1];
            pt[2] = 0.;
            pts->SetPoint(ptIdx, pt);
            var->SetTuple1(ptIdx, speeds[i][j]);
            if (j > 0)
            {
                vtkIdType ids[2] = { ptIdx-1, ptIdx };
                lines->InsertNextCell(2, ids);
            }
            ptIdx++;
        }
    }

    vtkPolyData *pd = vtkPolyData::New();
    pd->SetPoints(pts);
    pd->GetPointData()->AddArray(var);
    pd->GetPointData()->SetActiveScalars("speed");
    pd->SetLines(lines);
    lines->Delete();
    var->Delete();
    pts->Delete();

    return pd;
}

int main()
{
    int  i, j;

    vtkDataSetReader *rdr = vtkDataSetReader::New();
    rdr->SetFileName("proj4_data.vtk");
    rdr->Update();

    int dims[3];
    vtkRectilinearGrid *rgrid = (vtkRectilinearGrid *) rdr->GetOutput();
    rgrid->GetDimensions(dims);

    float *X = (float *) rgrid->GetXCoordinates()->GetVoidPointer(0);
    float *Y = (float *) rgrid->GetYCoordinates()->GetVoidPointer(0);
    float *F_3D = (float *) rgrid->GetPointData()->GetVectors()->GetVoidPointer(0);
    float *F = Convert3DVectorDataTo2DVectorData(dims, F_3D);
    
    
    const int npts = 10;
    float pt[npts][3] =
         {
            {10.1119, 1.22062, 0},
            {8.62376, 13.3839, 0},
            {1.55026, 1.26123, 0},
            {6.9736, 0.653565, 0},
            {2, 2.74117, 0},
            {8.93699, 10.4111, 0},
            {6.08791, -0.533753, 0},
            {10.0543, 1.38024, 0},
            {3.84128, -0.768977, 0},
            {6.66757, 6.0259, 0},
         };


    for (i = 0 ; i < npts ; i++)
    {
       float vec[2];
       EvaluateVectorFieldAtLocation(pt[i], dims, X, Y, F, vec);
       cerr << "Velocity at (" << pt[i][0] <<", "<<pt[i][1] << ") is (" << vec[0] << ", " << vec[1] << ")" << endl;
    }

    float h = 0.01;
    const int nsteps = 5000;
    float **output_locations = new float*[2*(npts+1)];
    float **speeds = new float*[npts+1];
    for (i = 0 ; i < npts ; i++)
    {
       output_locations[i] = new float[(nsteps+1)*2];
       speeds[i] = new float[nsteps];
       AdvectWithEulerStep(pt[i], dims, X, Y, F, h, nsteps, output_locations[i], speeds[i]);
       float length = CalculateArcLength(output_locations[i], nsteps+1);
       cerr << "Arc length for (" << pt[i][0] << ", " << pt[i][1] << ") is " << length << endl;
    }

    vtkPolyData *pd = CreateVTKPolyData(npts, nsteps, output_locations, speeds);

    //This can be useful for debugging
/*
    vtkDataSetWriter *writer = vtkDataSetWriter::New();
    writer->SetFileName("paths.vtk");
    writer->SetInput(pd);
    writer->Write();
 */

    vtkSmartPointer<vtkDataSetMapper> win1Mapper =
      vtkSmartPointer<vtkDataSetMapper>::New();
    win1Mapper->SetInputData(pd);
    win1Mapper->SetScalarRange(0, 0.15);
  
    vtkSmartPointer<vtkActor> win1Actor =
      vtkSmartPointer<vtkActor>::New();
    win1Actor->SetMapper(win1Mapper);
  
    vtkSmartPointer<vtkRenderer> ren1 =
      vtkSmartPointer<vtkRenderer>::New();
  
    vtkSmartPointer<vtkRenderWindow> renWin =
      vtkSmartPointer<vtkRenderWindow>::New();
    renWin->AddRenderer(ren1);
  
    vtkSmartPointer<vtkRenderWindowInteractor> iren =
      vtkSmartPointer<vtkRenderWindowInteractor>::New();
    iren->SetRenderWindow(renWin);
    ren1->AddActor(win1Actor);
    ren1->SetBackground(0.0, 0.0, 0.0);
    renWin->SetSize(800, 800);
  
    ren1->GetActiveCamera()->SetFocalPoint(5,5,0);
    ren1->GetActiveCamera()->SetPosition(5,5,30);
    ren1->GetActiveCamera()->SetViewUp(0,1,0);
    ren1->GetActiveCamera()->SetClippingRange(20, 120);
    ren1->GetActiveCamera()->SetDistance(30);

    // This starts the event loop and invokes an initial render.
    //
    iren->Initialize();
    iren->Start();

    delete [] F;
    pd->Delete();
}

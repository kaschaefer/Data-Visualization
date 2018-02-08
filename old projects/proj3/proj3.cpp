#include <vtkImageData.h>
#include <vtkPNGWriter.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkUnsignedCharArray.h>
#include <vtkFloatArray.h>
#include <vtkDataSetReader.h>
#include <vtkRectilinearGrid.h>
#include <vtkFloatArray.h>


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
getT(const float X, const float A, const float B) {
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
InterpolateForColor(const float Fa, const float Fb, const float X, const float A, const float B, unsigned char *returnVal)
{
    float t = getT(X, A, B);
    *returnVal = (Fa + (t * (Fb-Fa)));
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
    float t = getT(X, A, B);
    *returnVal = (Fa + (t * (Fb-Fa)));
}

// ****************************************************************************
//  Function: EvaluateFieldAtLocation
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
//     F: a scalar field defined on the mesh.  Its size is dims[0]*dims[1].
//
//   Returns: the interpolated field value. 0 if the location is out of bounds.
//
// ****************************************************************************

float
EvaluateFieldAtLocation(const float *pt, const int *dims, 
                        const float *X, const float *Y, const float *F)
{
    float returnVal = 0;
    
    //Get the Cell that contains point P
    int a = GetCellFromVal(X, dims, pt[0]);
    int b = GetCellFromVal(Y, &dims[1], pt[1]);

    //If the point P is within the mesh, then continue
    if (a != -1 && b != -1) {
        float f0 = 0.0;
        float f1 = 0.0;
        float f2 = 0.0;
        float f3 = 0.0;
        int v0[2];
        int v1[2];
        int v2[2];
        int v3[2];
        //Get Cell's Vertices v0 v1 v2 v3
        
        v0[0] = a;
        v0[1] = b;
        
        v1[0] = a+1;
        v1[1] = b;
        
        v2[0] = a;
        v2[1] = b+1;

        v3[0] = a+1;
        v3[1] = b+1;

        // Find F(v0) F(v1) F(v2) F(v3)
        f0 = F[GetPointIndex(v0, dims)];
        f1 = F[GetPointIndex(v1, dims)];
        f2 = F[GetPointIndex(v2, dims)];
        f3 = F[GetPointIndex(v3, dims)];

        //Perform bilinear interpolation to location P

            float line1, line2, line3;
            //Interpolate between v0 and v1 to get I1
            Interpolate(f0, f1, pt[0], X[v0[0]], X[v1[0]], &line1);
            //Interpolate between v2 and v3 to get I4
            Interpolate(f2, f3, pt[0], X[v2[0]], X[v3[0]], &line2);
            //Interpolate between v0 and v2 to get I2
            Interpolate(line1, line2, pt[1], Y[v0[1]], Y[v2[1]], &line3);
            returnVal = line3;

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

// ****************************************************************************
//  Function: ApplyBlueHotColorMap
//
//  Purpose: 
//     Maps a normalized scalar value F (0<=F<=1) to a color using the blue 
//     hot color map.
//
//     The blue hot color map has:
//        F=0: (0,0,128) 
//        F=1: (255,255,255) 
//       and smooth interpolation in between
//
//  Arguments:
//       F (input):     a scalar value between 0 and 1
//       RGB (output):  the location to store the color
//      
// ****************************************************************************

void
ApplyBlueHotColorMap(float F, unsigned char *RGB)
{
    InterpolateForColor(0.0, 255.0, F, 0.0, 1.0, RGB);
    InterpolateForColor(0.0, 255.0, F, 0.0, 1.0, (RGB+1));
    InterpolateForColor(128.0, 255.0, F, 0.0, 1.0, (RGB+2));
}


// ****************************************************************************
//  Function: ApplyDifferenceColorMap
//
//  Purpose: 
//     Maps a normalized scalar value F (0<=F<=1) to a color using a divergent colormap
//
//     The divergent color map has:
//        F=0: (0,0,128) 
//        F=0.5: (255,255,255) 
//        F=1: (128, 0, 0)
//       and smooth interpolation in between
//
//  Arguments:
//       F (input):     a scalar value between 0 and 1
//       RGB (output):  the location to store the color
//      
// ****************************************************************************
void
ApplyDifferenceColorMap(float F, unsigned char *RGB)
{
    if (F < 0.5) {
        InterpolateForColor(0.0, 255.0, F, 0.0, 0.5, RGB);
        InterpolateForColor(0.0, 255.0, F, 0.0, 0.5, (RGB+1));
        InterpolateForColor(128.0, 255.0, F, 0.0, 0.5, (RGB+2));
    }
    else {
        InterpolateForColor(255.0, 128.0, F, 0.5, 1.0, RGB);
        InterpolateForColor(255.0, 0.0, F, 0.5, 1.0, (RGB+1));
        InterpolateForColor(255.0, 0.0, F, 0.5, 1.0, (RGB+2));
    }
}

// ****************************************************************************
//  Function: ApplyBHSVColorMap
//
//  Purpose: 
//     Maps a normalized scalar value F (0<=F<=1) to a color using an HSV rainbow colormap
//
//     The rainbow colormap uses a saturation =1.0, value = 1.0, 
//     and interpolates hue from 0 to 360 degrees 
//
//  Arguments:
//       F (input):     a scalar value between 0 and 1
//       RGB (output):  the location to store the color
//      
// ****************************************************************************
void
ApplyHSVColorMap(float F, unsigned char *RGB)
{
    float hue = 0.0;
    float saturation = 1.0;
    float val = 1.0;
    float p = 0.0;
    float q = 0.0;
    float t = 0.0;
    int i = 0;
    float f = 0.0;

    hue = 360*F;

    hue /= 60.f;
    i = floor(hue);
    f = hue - i;

    p = val * (1 - saturation);
    q = val * (1 - saturation*f);
    t = val * (1 - saturation* (1-f));
    
    switch(i)
    {
        case 0:
            *RGB = val*255.f;
            *(RGB+1) = t*255.f;
            *(RGB+2) = p*255.f;
            break;
        case 1:
            *RGB = q*255.f;
            *(RGB+1) = val*255.f;
            *(RGB+2) = p*255.f;
            break;
        case 2:
            *RGB = p*255.f;
            *(RGB+1) = val*255.f;
            *(RGB+2) = t*255.f;
            break;
        case 3:
            *RGB = p*255.f;
            *(RGB+1) = q*255.f;
            *(RGB+2) = val*255.f;
            break;
        case 4:
            *RGB = t*255.f;
            *(RGB+1) = p*255.f;
            *(RGB+2) = val*255.f;
            break;
        case 5:
            *RGB = val*255.f;
            *(RGB+1) = p*255.f;
            *(RGB+2) = q*255.f;
            break;
    }

}


int main()
{
    int  i, j;

    vtkDataSetReader *rdr = vtkDataSetReader::New();
    rdr->SetFileName("proj3_data.vtk");
    rdr->Update();

    int dims[3];
    vtkRectilinearGrid *rgrid = (vtkRectilinearGrid *) rdr->GetOutput();
    rgrid->GetDimensions(dims);

    float *X = (float *) rgrid->GetXCoordinates()->GetVoidPointer(0);
    float *Y = (float *) rgrid->GetYCoordinates()->GetVoidPointer(0);
    float *F = (float *) rgrid->GetPointData()->GetScalars()->GetVoidPointer(0);
    
    int nx = 500;
    int ny = 500;

    vtkImageData *images[3];
    unsigned char *buffer[3];
    for (i = 0 ; i < 3 ; i++)
    {
        images[i] = NewImage(nx, ny);
        buffer[i] = (unsigned char *) images[i]->GetScalarPointer(0,0,0);
    }

    for (i = 0 ; i < 3*nx*ny ; i++)
        for (j = 0 ; j < 3 ; j++)
            buffer[j][i] = 0;

    for (i = 0 ; i < nx ; i++)
        for (j = 0 ; j < ny ; j++)
        {
            // ITERATE OVER PIXELS
            float aNum = (float)i;
            float anotherNum = (float)j;
            float pt[2];
            pt[0] = aNum/499*18+(-9);
            pt[1] = anotherNum/499*18+(-9);
            float f = EvaluateFieldAtLocation(pt, dims, X, Y, F);
            float normalizedF = getT(f, 1.20, 5.02); 
            // I TAKE OVER HERE
            int offset = 3*(j*nx+i);
            ApplyBlueHotColorMap(normalizedF, buffer[0]+offset);
            ApplyDifferenceColorMap(normalizedF, buffer[1]+offset);
            ApplyHSVColorMap(normalizedF, buffer[2]+offset);
        }

    WriteImage(images[0], "bluehot");
    WriteImage(images[1], "difference");
    WriteImage(images[2], "hsv");
}

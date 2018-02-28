/*=========================================================================

  Program:   Visualization Toolkit
  Module:    SpecularSpheres.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
//
// This examples demonstrates the effect of specular lighting.
//
#include "vtkSmartPointer.h"
#include "vtkSphereSource.h"
#include "vtkPolyDataMapper.h"
#include "vtkActor.h"
#include "vtkInteractorStyle.h"
#include "vtkObjectFactory.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkProperty.h"
#include "vtkCamera.h"
#include "vtkLight.h"
#include "vtkOpenGLPolyDataMapper.h"
#include "vtkJPEGReader.h"
#include "vtkImageData.h"
#include <vtkPNGWriter.h>

#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkPolyDataReader.h>
#include <vtkPoints.h>
#include <vtkUnsignedCharArray.h>
#include <vtkFloatArray.h>
#include <vtkDoubleArray.h>
#include <vtkCellArray.h>
#include <vtkDataSetReader.h>
#include <vtkContourFilter.h>
#include <vtkRectilinearGrid.h>

#include <vtkCamera.h>
#include <vtkDataSetMapper.h>
#include <vtkRenderer.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkSmartPointer.h>

#include "TriangleList.h"
#include "tricase.cpp"
#include <cmath>


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
    return dims[0]*dims[1]*dims[2];
    // 2D
    //return dims[0]*dims[1];
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
    return (dims[0]-1)*(dims[1]-1)*(dims[2]-1);
    // 2D
   // return (dims[0]-1)*(dims[1]-1);
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
    return idx[2]*dims[0]*dims[1]+idx[1]*dims[0]+idx[0];
    // 2D
    //return idx[1]*dims[0]+idx[0];
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
    return idx[2]*(dims[0]-1)*(dims[1]-1)+idx[1]*(dims[0]-1)+idx[0];
    // 2D
    //return idx[1]*(dims[0]-1)+idx[0];
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
    idx[0] = pointId%dims[0];
    idx[1] = (pointId/dims[0])%dims[1];
    idx[2] = pointId/(dims[0]*dims[1]);

    // 2D
    //idx[0] = pointId%dims[0];
    //idx[1] = pointId/dims[0];
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
    idx[0] = cellId%(dims[0]-1);
    idx[1] = (cellId/(dims[0]-1))%(dims[1]-1);
    idx[2] = cellId/((dims[0]-1)*(dims[1]-1));

    // 2D
    //idx[0] = cellId%(dims[0]-1);
    //idx[1] = cellId/(dims[0]-1);
}

/************************************
GenerateEdges function

Takes a reference to 4x2 integer array
Fills the Array with which vertices correspond to which edge
Ex: Edge 0 corresponds to vertices 1 and 0
    Edge 2 corersponds to vertices 2 and 1
Returns Nothing
*************************************/
//Change this??????????
void
GenerateEdges(int arr[][2]) {
    //edge 0
    arr[0][0] = 0;
    arr[0][1] = 1;
    //edge 1
    arr[1][0] = 1;
    arr[1][1] = 3;
    //edge2
    arr[2][0] = 2;
    arr[2][1] = 3;
    //edge3
    arr[3][0] = 0;
    arr[3][1] = 2;
    //edge4
    arr[4][0] = 4;
    arr[4][1] = 5;
    //edge5
    arr[5][0] = 5;
    arr[5][1] = 7;
    //edge6
    arr[6][0] = 6;
    arr[6][1] = 7;
    //edge7
    arr[7][0] = 4;
    arr[7][1] = 6;
    //edge8
    arr[8][0] = 0;
    arr[8][1] = 4;
    //edge9
    arr[9][0] = 1;
    arr[9][1] = 5;
    //edge10
    arr[10][0] = 2;
    arr[10][1] = 6;
    //edge11
    arr[11][0] = 3;
    arr[11][1] = 7;
    return;
}

/************************************
IdentifyCase Function

Takes a Cell
Returns the Integer Corresponding to the Isoline Case
*************************************/
int
IdentifyCase (const float isoValue, float *vertices)
{
    int counter = 0;
    for (int i = 0; i < 8; i++) {
        if (vertices[i] > isoValue){
            counter += pow(2, i);
        }
    }
    return counter;
}

int main()
{
    int  i, j;
    const float isoValue = 3.20;
    int numCells = 0;
    int iCase = 0;

    //Lookup Tables
    int edgeTable[12][2];
    int pointIndex = 0;
    float scalarValsAtVertices[8];

    vtkDataSetReader *rdr = vtkDataSetReader::New();
    rdr->SetFileName("proj7.vtk");
    rdr->Update();

    int dims[3];
    vtkRectilinearGrid *rgrid = (vtkRectilinearGrid *) rdr->GetOutput();
    rgrid->GetDimensions(dims);

    float *X = (float *) rgrid->GetXCoordinates()->GetVoidPointer(0);
    float *Y = (float *) rgrid->GetYCoordinates()->GetVoidPointer(0);
    float *Z = (float *) rgrid->GetZCoordinates()->GetVoidPointer(0);
    float *F = (float *) rgrid->GetPointData()->GetScalars()->GetVoidPointer(0);
    
    // Make a Triangle List

    TriangleList t1;
// YOUR CODE TO GENERATE ISOLINES SHOULD GO HERE!
    //Generate LookUp Table and num segments table
    //Tricase[256][16] is the lookup Table for edges once case has been determined
    numCells = GetNumberOfCells(dims);
    GenerateEdges(edgeTable);
    
    //Stuff to be moved up later
    int idx[3];
    int logicalVertices[8][3];
    for (i=0; i < numCells; i++) {
        //Store logical point indices for V0
        GetLogicalCellIndex(logicalVertices[0], i, dims);
        //V0
        idx[0] = logicalVertices[0][0];
        idx[1] = logicalVertices[0][1];
        idx[2] = logicalVertices[0][2];

        //V1
        logicalVertices[1][0] = idx[0]+1;
        logicalVertices[1][1] = idx[1];
        logicalVertices[1][2] = idx[2];

        //V2
        logicalVertices[2][0] = idx[0];
        logicalVertices[2][1] = idx[1]+1;
        logicalVertices[2][2] = idx[2];

        //V3
        logicalVertices[3][0] = idx[0]+1;
        logicalVertices[3][1] = idx[1]+1;
        logicalVertices[3][2] = idx[2];

        //V4
        logicalVertices[4][0] = idx[0];
        logicalVertices[4][1] = idx[1];
        logicalVertices[4][2] = idx[2]+1;

        //V5
        logicalVertices[5][0] = idx[0]+1;
        logicalVertices[5][1] = idx[1];
        logicalVertices[5][2] = idx[2]+1;

        //V6
        logicalVertices[6][0] = idx[0];
        logicalVertices[6][1] = idx[1]+1;
        logicalVertices[6][2] = idx[2]+1;

        //V7
        logicalVertices[7][0] = idx[0]+1;
        logicalVertices[7][1] = idx[1]+1;
        logicalVertices[7][2] = idx[2]+1;

        //Get Point Index of Each Vertex to Get Scalar Values
        for (j = 0; j < 8; j++) {
            pointIndex = GetPointIndex(logicalVertices[j], dims);
            scalarValsAtVertices[j] = F[pointIndex];
        }
        //Get Case
        iCase = IdentifyCase(isoValue, scalarValsAtVertices);

        //Get Edges for Triangles
        int *triEdges = triCase[iCase];

        //For each triangle (e.g. if *edges != -1)
        while (*triEdges != -1) {
            //Get the 3 edges
            int edgeArray[3];
            edgeArray[0] = *triEdges++;
            edgeArray[1] = *triEdges++;
            edgeArray[2] = *triEdges++;
            float pts[9];
            //For each edge that has a vertex of the triangle on it
            for (j=0; j < 3; j++) {
            //Get bounding Vertices [lookup Table]
                //save current edge in temp var to simplify below code
                int edge = edgeArray[j];
                int lowX = logicalVertices[edgeTable[edge][0]][0];
                int lowY = logicalVertices[edgeTable[edge][0]][1];
                int lowZ = logicalVertices[edgeTable[edge][0]][2];

                int highX = logicalVertices[edgeTable[edge][1]][0];
                int highY = logicalVertices[edgeTable[edge][1]][1];
                int highZ = logicalVertices[edgeTable[edge][1]][2];

                //Get Scalar Values
                float scalarLow = scalarValsAtVertices[edgeTable[edge][0]];
                float scalarHigh = scalarValsAtVertices[edgeTable[edge][1]];

                //Interpolate between them
                //pts[j*3] is the x value
                pts[j*3] = X[lowX] + ((isoValue-scalarLow) / (scalarHigh-scalarLow)) * (X[highX]-X[lowX]);
                //pts[j*3+1] is the y value
                pts[j*3+1] = Y[lowY] + ((isoValue-scalarLow) / (scalarHigh-scalarLow)) * (Y[highY]-Y[lowY]);
                //pts[j*3+2] is the z value
                pts[j*3+2] = Z[lowZ] + ((isoValue-scalarLow) / (scalarHigh-scalarLow)) * (Z[highZ]- Z[lowZ]);
            }

            //Add Triangle --  AddTriangle(float X1, float Y1, float Z1, float X2, float Y2, float Z2, float X3, float Y3, float Z3)
            t1.AddTriangle(pts[0], pts[1], pts[2], pts[3], pts[4], pts[5], pts[6], pts[7], pts[8]);
        }
    }


    vtkPolyData *pd = t1.MakePolyData();

    //This can be useful for debugging
/*
    vtkDataSetWriter *writer = vtkDataSetWriter::New();
    writer->SetFileName("paths.vtk");
    writer->SetInputData(pd);
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

    ren1->GetActiveCamera()->SetFocalPoint(0,0,0);
    ren1->GetActiveCamera()->SetPosition(0,0,50);
    ren1->GetActiveCamera()->SetViewUp(0,1,0);
    ren1->GetActiveCamera()->SetClippingRange(20, 120);
    ren1->GetActiveCamera()->SetDistance(30);

    // This starts the event loop and invokes an initial render.
    //
    iren->Initialize();
    iren->Start();

    pd->Delete();
}
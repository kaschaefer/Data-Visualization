
/************************************
IdentifyCase Function

Takes a Cell
Returns the Integer Corresponding to the Isoline Case
*************************************/
int
IdentifyCase (const float isoValue, float *vertices)
{
    int counter = 0;
    for (int i = 0; i < 4; i++) {
        if (vertices[i] > isoValue){
            counter += pow(2, i);
        }
    }
    return counter;
}

/**************************************
Generate Segments function

Takes a reference to an integer array (length 16)
the length is hardcoded into the program and thus not checked for/protected against array outof bounds
Fills the Array with NumSegments in Each Case
Returns Nothing

**************************************/
void
GenerateSegments(int* arr)
{
    arr[0] = 0;
    arr[1] = 1;
    arr[2] = 1;
    arr[3] = 1;
    arr[4] = 1;
    arr[5] = 1;
    arr[6] = 2;
    arr[7] = 1;
    arr[8] = 1;
    arr[9] = 2;
    arr[10] = 1;
    arr[11] = 1;
    arr[12] = 1;
    arr[13] = 1;
    arr[14] = 1;
    arr[15] = 0;
}

/************************************
GenerateTable function

Takes a reference to 16x4 integer array
Fills the Array with IsoSurface Case Information.
Returns Nothing
*************************************/
void
GenerateTable(int lup[][4]) 
{
    //Case 0
    lup[0][0] = lup[0][1] = lup[0][2] = lup[0][3] = -1;
    //Case 1
    lup[1][0] = 0;
    lup[1][1] = 3;
    lup[1][2] = lup[1][3] = -1;
    //Case 2
    lup[2][0] = 0;
    lup[2][1] = 1;
    lup[2][2] = lup[2][3] = -1;
    //Case 3
    lup[3][0] = 1;
    lup[3][1] = 3;
    lup[3][2] = lup[3][3] = -1;
    //Case 4
    lup[4][0] = 2;
    lup[4][1] = 3;
    lup[4][2] = lup[4][3] = -1;
    //Case 5
    lup[5][0] = 0;
    lup[5][1] = 2;
    lup[5][2] = lup[5][3] = -1;
    //Case 6
    lup[6][0] = 0;
    lup[6][1] = 1;
    lup[6][2] = 2;
    lup[6][3] = 3;
    //Case 7
    lup[7][0] = 1;
    lup[7][1] = 2;
    lup[7][2] = lup[7][3] = -1;
    //Case 8
    lup[8][0] = 1;
    lup[8][1] = 2;
    lup[8][2] = lup[8][3] = -1;
    //Case 9
    lup[9][0] = 0;
    lup[9][1] = 3;
    lup[9][2] = 1;
    lup[9][3] = 2; 
    //Case 10
    lup[10][0] = 0;
    lup[10][1] = 2;
    lup[10][2] = lup[10][3] = -1;
    //Case 11
    lup[11][0] = 2;
    lup[11][1] = 3;
    lup[11][2] = lup[11][3] = -1;
    //Case 12
    lup[12][0] = 1;
    lup[12][1] = 3;
    lup[12][2] = lup[12][3] = -1;
    //Case 13
    lup[13][0] = 0;
    lup[13][1] = 1;
    lup[13][2] = lup[13][3] = -1;
    //Case 14
    lup[14][0] = 0;
    lup[14][1] = 3;
    lup[14][2] = lup[14][3] = -1;
    //Case 15
    lup[15][0] = lup[15][1] = lup[15][2] = lup[15][3] = -1;
}

/************************************
GenerateEdges function

Takes a reference to 4x2 integer array
Fills the Array with which vertices correspond to which edge
Ex: Edge 0 corresponds to vertices 1 and 0
    Edge 2 corersponds to vertices 2 and 1
Returns Nothing
*************************************/

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
    return;
}

//      THIS CODE IS PART OF MAIN()

// YOUR CODE TO GENERATE ISOLINES SHOULD GO HERE!
    //Generate LookUp Table and num segments table
    GenerateTable(lookupTable);
    GenerateSegments(numSegments);
    GenerateEdges(edges);
    numCells = GetNumberOfCells(dims);
    
    //Stuff to be moved up later
    int idx[2];
    int logicalVertices[4][2];
    int pointIndices[4];
    for (i=0; i < numCells; i++) {
        //Store logical point indices for V0
        GetLogicalCellIndex(logicalVertices[0], i, dims);

        //V1
        logicalVertices[1][0] = idx[0]+1;
        logicalVertices[1][1] = idx[1];

        //V2
        logicalVertices[2][0] = idx[0];
        logicalVertices[2][1] = idx[1]+1;

        //V3
        logicalVertices[3][0] = idx[0]+1;
        logicalVertices[3][1] = idx[1]+1;

        //Get Point Indices
        pointIndices[0] = GetPointIndex(logicalVertices[0], dims);
        pointIndices[1] = GetPointIndex(logicalVertices[1], dims);
        pointIndices[2] = GetPointIndex(logicalVertices[2], dims);
        pointIndices[3] = GetPointIndex(logicalVertices[3], dims);

        //Get Scalar Values at Vertices
        scalarValsAtVertices[0] = F[pointIndices[0]];
        scalarValsAtVertices[1] = F[pointIndices[1]];
        scalarValsAtVertices[2] = F[pointIndices[2]];
        scalarValsAtVertices[3] = F[pointIndices[3]];

        //Get Case
        iCase = IdentifyCase(isoValue, scalarValsAtVertices);

        //Get Num Segments
        nSegment = numSegments[iCase];
        if (nSegment != 0) {
            for (j=0; j < nSegment; j++) {
                float pt1[2];
                int edge1 = lookupTable[iCase][2*j];

                //Get Logic
                int lowX = logicalVertices[edges[edge1][0]][0];
                int lowY = logicalVertices[edges[edge1][0]][1];

                int highX = logicalVertices[edges[edge1][1]][0];
                int highY = logicalVertices[edges[edge1][1]][1];

                float scalarLow = scalarValsAtVertices[edges[edge1][0]];
                float scalarHigh = scalarValsAtVertices[edges[edge1][1]];

                
                
                pt1[0] = X[lowX] + ((isoValue-scalarLow) / (scalarHigh-scalarLow)) * (X[highX]-X[lowX]);
                pt1[1] = Y[lowY] + ((isoValue-scalarLow) / (scalarHigh-scalarLow)) * (Y[highY]-Y[lowY]);

                int edge2 = lookupTable[iCase][2*j+1];
                float pt2[2];

                lowX = logicalVertices[edges[edge2][0]][0];
                lowY = logicalVertices[edges[edge2][0]][1];

                highX = logicalVertices[edges[edge2][1]][0];
                highY = logicalVertices[edges[edge2][1]][1];

                scalarLow = scalarValsAtVertices[edges[edge2][0]];
                scalarHigh = scalarValsAtVertices[edges[edge2][1]];

                pt2[0] = X[lowX] + ((isoValue-scalarLow) / (scalarHigh-scalarLow)) * (X[highX]-X[lowX]);
                pt2[1] = Y[lowY] + ((isoValue-scalarLow) / (scalarHigh-scalarLow)) * (Y[highY]-Y[lowY]);
                
                sl.AddSegment(pt1[0], pt1[1], pt2[0], pt2[1]);
            }
        }


    }

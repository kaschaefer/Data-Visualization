#include <vtkDataSetReader.h>
#include <vtkContourFilter.h>
#include <vtkDataSetWriter.h>
#include <vtkSmartPointer.h>
#include <vtkActor.h>
#include <vtkRenderer.h>

int main(int argc, char* argv[]){
    //Get DataSet Reader
    vtkDataSetReader *rdr = vtkDataSetReader::New();
    rdr->SetFileName("proj8.vtk");
    rdr->Update();

    //hardyglobal is a field F
    vtkContourFilter *cf = vtkContourFilter::New();
    cf->SetValue(2, 2.5, 5.0);


    //Contour Filter for IsoSurfacing Renderer 1
        //set active attribute for contour filter?

    //Slice Filter for Renderer 2

    //grad is a field F
    //Get NEW F

    //Hedgehog filter for Renderer 3
        //hedgehog filter -- needs to be pointed at the data
        //Set active attribute

    //Streamlines filter for Renderer 4


    //polyDataMapper
    //2 actors per isosurface


    //color map colors need to be pointed at data set by data set mapper
    //set scalar range -- for color map

    //print method?????




    //Make Shapes


    vtkSmartPointer<vtkRenderWindow> rendWindow = vtkSmartPointer<vtkRenderWindow>::New();
    vtkSmartPointer<vtkRenderWindowInteractor> iren = vtkSmartPointer<vtkRenderWindowInteractor>::New();

    //Set Up Renderers
    vtkSmartPointer<vtkRenderer> ren1 = vtkSmartPointer<vtkRenderer>::New();
    vtkSmartPointer<vtkRenderer> ren2 = vtkSmartPointer<vtkRenderer>::New();
    vtkSmartPointer<vtkRenderer> ren3 = vtkSmartPointer<vtkRenderer>::New();
    vtkSmartPointer<vtkRenderer> ren4 = vtkSmartPointer<vtkRenderer>::New();

    //Set Viewports
    ren1->SetViewport(0, 0, .5, .5);
    ren2->SetViewport(0, 0.5, 0, 0.5);
    ren3->SetViewport(0.5, 0.5, 1.0, 1.0);
    ren4->SetViewport(0.5, 1.0, 0.5, 1.0);

    //Add Renderers to Render Window
    rendWindow->AddRenderer(ren1);
    rendWindow->AddRenderer(ren2);
    rendWindow->AddRenderer(ren3);
    rendWindow->AddRenderer(ren4);

    //Add Interactor
    renderWindowInteractor->SetRenderWindow(rendWindow);

    //Add Actor to Renderer, Set Background and Size

    //Invokes Initial Render
    iren->Initiate();
    iren->Start();

    return EXIT_SUCCESS;

}
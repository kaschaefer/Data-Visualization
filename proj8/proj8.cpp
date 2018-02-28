#include <vtkDataSetReader.h>
#include <vtkContourFilter.h>
#include <vtkDataSetWriter.h>

int main(){
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
    

    //Add the stuff to its respective renderer

    //Invokes Initial Render
    iren->Initiate();
    iren->Start();

    return EXIT_SUCCESS;

}
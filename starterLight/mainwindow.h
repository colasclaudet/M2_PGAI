#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QFileDialog>
#include <QMainWindow>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>

namespace Ui {
class MainWindow;
}

using namespace OpenMesh;
using namespace OpenMesh::Attributes;

struct MyTraits : public OpenMesh::DefaultTraits
{
    // use vertex normals and vertex colors
    VertexAttributes( OpenMesh::Attributes::Normal | OpenMesh::Attributes::Color | OpenMesh::Attributes::Status);
    // store the previous halfedge
    HalfedgeAttributes( OpenMesh::Attributes::PrevHalfedge );
    // use face normals face colors
    FaceAttributes( OpenMesh::Attributes::Normal | OpenMesh::Attributes::Color | OpenMesh::Attributes::Status);
    EdgeAttributes( OpenMesh::Attributes::Color | OpenMesh::Attributes::Status );
    // vertex thickness
    VertexTraits{float thickness; float value;};
    // edge thickness
    EdgeTraits{float thickness;};
};
typedef OpenMesh::TriMesh_ArrayKernelT<MyTraits> MyMesh;

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:

    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

    void displayMesh(MyMesh *_mesh, bool isTemperatureMap = false, float mapRange = -1);
    void resetAllColorsAndThickness(MyMesh* _mesh);
    float compute_face_area(MyMesh* _mesh, int n_face);
    float compute_area(MyMesh* _mesh);
    float get_min_area(MyMesh* _mesh);
    float get_max_area(MyMesh* _mesh);
    float compute_face_normal(MyMesh* _mesh, int n_face);
    float compute_normal(MyMesh* _mesh, int n_face);


protected:
    void mousePressEvent(QMouseEvent * e);

private slots:
    void on_pushButton_chargement_clicked();

    void on_pushButton_generer_clicked();

    void on_pushButton_area_clicked();

    void on_pushButton_clicked();

private:

    MyMesh mesh;

    Ui::MainWindow *ui;
};

#endif // MAINWINDOW_H

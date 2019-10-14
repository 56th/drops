/// \author Alexander Zhiliakov alex@math.uh.edu

#ifndef VTK_WRITER_HPP
#define VTK_WRITER_HPP

#ifdef _VTK
    #include "vtkSmartPointer.h"
    #include "vtkPoints.h"
    #include "vtkLagrangeTetra.h"
    #include "vtkCellArray.h"
    #include "vtkUnstructuredGrid.h"
    #include "vtkXMLUnstructuredGridWriter.h"
    #include "vtkDoubleArray.h"
    #include "vtkPointData.h"
#endif
#include "geom/multigrid.h"

namespace DROPS {
    class VTKWriter {
    public:
        struct VTKVar {
            enum class Type { P2, vecP2, P1, vecP1 };
            std::string name;
            VectorCL* value;
            Type type;
        };
    private:
        bool binary;
        std::vector<VTKVar> vars;
        MultiGridCL const * mg;
        std::string path;
        std::unordered_map<VertexCL const *, size_t> vertexIndex;
        std::unordered_map<EdgeCL const *, size_t> edgeIndex;
        size_t frame = 0;
        #ifdef _VTK
            vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
        #endif
    public:
        VTKWriter(std::string const & path_, MultiGridCL const & mg_, bool binary_ = true) : path(path_), mg(&mg_), binary(binary_) {
            #ifdef _VTK
                size_t n = 0;
                auto points = vtkSmartPointer<vtkPoints>::New();
                for (auto it = mg->GetTriangVertexBegin(); it != mg->GetTriangVertexEnd(); ++it) {
                    points->InsertNextPoint(it->GetCoord()[0], it->GetCoord()[1], it->GetCoord()[2]);
                    vertexIndex[&*it] = n++;
                }
                for (auto it = mg->GetTriangEdgeBegin(); it != mg->GetTriangEdgeEnd(); ++it) {
                    auto p = GetBaryCenter(*it);
                    points->InsertNextPoint(p[0], p[1], p[2]);
                    edgeIndex[&*it] = n++;
                }
                auto cellArray = vtkSmartPointer<vtkCellArray>::New();
                for (auto it = mg->GetTetrasBegin(); it != mg->GetTetrasEnd(); ++it) {
                    auto tetra = vtkSmartPointer<vtkLagrangeTetra>::New();
                    tetra->GetPointIds()->SetNumberOfIds(10);
                    tetra->GetPoints()->SetNumberOfPoints(10);
                    tetra->Initialize();
                    for (size_t i : {0, 1, 2, 3})
                        tetra->GetPointIds()->SetId(i, vertexIndex[it->GetVertex(i)]);
                    // for nodes ordering, cf. https://vtk.org/doc/nightly/html/classvtkQuadraticTetra.html#details and https://blog.kitware.com/modeling-arbitrary-order-lagrange-finite-elements-in-the-visualization-toolkit/
                    tetra->GetPointIds()->SetId(4, edgeIndex[it->GetEdge(EdgeByVert(0, 1))]);
                    tetra->GetPointIds()->SetId(5, edgeIndex[it->GetEdge(EdgeByVert(1, 2))]);
                    tetra->GetPointIds()->SetId(6, edgeIndex[it->GetEdge(EdgeByVert(0, 2))]);
                    tetra->GetPointIds()->SetId(7, edgeIndex[it->GetEdge(EdgeByVert(0, 3))]);
                    tetra->GetPointIds()->SetId(8, edgeIndex[it->GetEdge(EdgeByVert(1, 3))]);
                    tetra->GetPointIds()->SetId(9, edgeIndex[it->GetEdge(EdgeByVert(2, 3))]);
                    cellArray->InsertNextCell(tetra);
                }
                unstructuredGrid->SetPoints(points);
                unstructuredGrid->SetCells(VTK_LAGRANGE_TETRAHEDRON, cellArray);
            #endif
        }
        VTKWriter& add(VTKVar const & var) {
            vars.push_back(var);
            return *this;
        }
        VTKWriter& write(double time) {
            std::string funcName = __func__;
            #ifdef _VTK
                // (1) update mesh, cf. https://lorensen.github.io/VTKExamples/site/Cxx/IO/WriteVTU/
                // vertexIndex.clear();
                // edgeIndex.clear();
                // ...
                // (2) update vars
                auto n = unstructuredGrid->GetNumberOfPoints();
                for (auto const & var : vars) {
                    if(var.type == VTKVar::Type::P2) {
                        if (var.value->size() != n)
                            throw std::invalid_argument(funcName + ": inconsistent number of d.o.f. for P2 interpolant");
                        auto array = vtkSmartPointer<vtkDoubleArray>::New();
                        array->SetArray(&var.value->operator[](0), n, 1);
                        array->SetName(var.name.c_str());
                        unstructuredGrid->GetPointData()->AddArray(array);
                    }
                    else if (var.type == VTKVar::Type::P1) {
                        auto  array = vtkSmartPointer<vtkDoubleArray>::New();
                        std::vector<double> value;
                        value.assign(&var.value->operator[](0), &var.value->operator[](0) + var.value->size());
                        for (auto it = mg->GetTriangEdgeBegin(); it != mg->GetTriangEdgeEnd(); ++it) {
                            auto i1 = vertexIndex[it->GetVertex(0)];
                            auto i2 = vertexIndex[it->GetVertex(1)];
                            value.push_back(.5 * (value[i1] + value[i2]));
                        }
                        if (value.size() != n) throw std::invalid_argument(funcName + ": inconsistent number of d.o.f. for P1 interpolant");
                        auto* ptr = new double[n];
                        std::copy(value.begin(), value.end(), ptr);
                        array->SetArray(ptr, n, 0 /* this will free the memory for ptr, cf. https://vtk.org/doc/release/5.6/html/a00505.html#d35ae5bee4aa873d543f1ab3eaf94454 */);
                        array->SetName(var.name.c_str());
                        unstructuredGrid->GetPointData()->AddArray(array);
                    }
                    else if (var.type == VTKVar::Type::vecP2) {
                        if (var.value->size() != 3 * n)
                            throw std::invalid_argument(funcName + ": inconsistent number of d.o.f. for vecP2 interpolant");
                        auto array = vtkSmartPointer<vtkDoubleArray>::New();
                        array->SetNumberOfComponents(3); // cf. https://vtk.org/Wiki/VTK/Examples/Cxx/PolyData/PolyDataCellNormals
                        array->SetArray(&var.value->operator[](0), 3 * n, 1);
                        array->SetName(var.name.c_str());
                        unstructuredGrid->GetPointData()->AddArray(array);
                    }
                    else if (var.type == VTKVar::Type::vecP1) {
                        std::vector<double> value;
                        value.assign(&var.value->operator[](0), &var.value->operator[](0) + var.value->size());
                        for (auto it = mg->GetTriangEdgeBegin(); it != mg->GetTriangEdgeEnd(); ++it) {
                            auto i1 = 3 * vertexIndex[it->GetVertex(0)];
                            auto i2 = 3 * vertexIndex[it->GetVertex(1)];
                            value.push_back(.5 * (value[i1] + value[i2]));
                            value.push_back(.5 * (value[i1 + 1] + value[i2 + 1]));
                            value.push_back(.5 * (value[i1 + 2] + value[i2 + 2]));
                        }
                        if (value.size() != 3 * n)
                            throw std::invalid_argument(funcName + ": inconsistent number of d.o.f. for vecP1 interpolant");
                        auto* ptr = new double[3 * n];
                        std::copy(value.begin(), value.end(), ptr);
                        auto  array = vtkSmartPointer<vtkDoubleArray>::New();
                        array->SetNumberOfComponents(3); // cf. https://vtk.org/Wiki/VTK/Examples/Cxx/PolyData/PolyDataCellNormals
                        array->SetArray(ptr, 3 * n, 0 /* this will free the memory for ptr, cf. https://vtk.org/doc/release/5.6/html/a00505.html#d35ae5bee4aa873d543f1ab3eaf94454 */);
                        array->SetName(var.name.c_str());
                        unstructuredGrid->GetPointData()->AddArray(array);
                    }
                    else
                        throw std::logic_error(funcName + ": export for " + typeid(var.type).name() + " not implemented");
                }
                // (3) write
                auto writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
                auto name = path + '_' + std::to_string(frame) + ".vtu";
                writer->SetFileName(name.c_str());
                writer->SetInputData(unstructuredGrid);
                if (!binary) writer->SetDataModeToAscii();
                writer->Write();
                // (4) update .pvd
                name = name.substr(name.find_last_of("/\\") + 1);
                if  (frame == 0) {
                    std::ofstream pvd(path + ".pvd");
                    pvd <<
                        "<?xml version=\"1.0\"?>\n"
                        "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
                        "<Collection>\n"
                            "\t<DataSet timestep=\""<< time <<"\" group=\"\" part=\"0\" file=\"" << name <<"\"/>\n"
                        "</Collection>\n"
                        "</VTKFile>";
                } else {
                    std::ofstream pvd;
                    pvd.open(path + ".pvd", std::ios_base::in);
                    if(!pvd.is_open()) std::logic_error(funcName + ": cannot reopen .pvd file");
                    pvd.seekp(-24, std::ios_base::end);
                    pvd <<
                            "\t<DataSet timestep=\""<< time <<"\" group=\"\" part=\"0\" file=\"" << name <<"\"/>\n"
                        "</Collection>\n"
                        "</VTKFile>";
                    pvd.close();
                }
                ++frame;
                return *this;
            #else
                throw std::logic_error(funcName + ": DROPS compiled w/o VTK lib");
            #endif
        }
    };
}

#endif // VTK_WRITER_HPP
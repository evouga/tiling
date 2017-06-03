#include "OpenMeshGlobals.h"
#include "IsotropicRemeshing.h"
#include <iostream>

using namespace std;


// Try it with
// Remeshing.exe bunny.obj bunny-coarse.obj 0.01
// Remeshing.exe bunny.obj bunny-fine.obj 0.005

int main(int argc, char *argv[])
{
    if(argc < 4 || argc > 5)
    {
        cerr << "Usage: Remeshing (input .obj) (output .obj) (desired edge lengths) [num iterations = 10]" << endl;
        return -1;
    }

    TriangleMesh inputMesh;
    if(!OpenMesh::IO::read_mesh<TriangleMesh>(inputMesh, argv[1]))
    {
        cerr << "Couldn't load input file " << argv[1] << endl;
        return -1;
    }

    inputMesh.update_face_normals();
    inputMesh.update_vertex_normals();

    inputMesh.request_vertex_status();
    inputMesh.request_edge_status();
    inputMesh.request_face_status();

    int numiters = 10;
    if(argc == 5)
    {
        numiters = strtol(argv[4], NULL, 10);
        if(numiters <= 0)
        {
            cerr << "Number of iterations must be positive" << endl;
            return -1;
        }
    }

    double elength = strtod(argv[3], NULL);
    if(elength <= 0)
    {
        cerr << "Desired edge length must be positive" << endl;
        return -1;
    }

    IsotropicRemeshing remesh(elength);
    remesh.remesh(&inputMesh, numiters);


    if(!OpenMesh::IO::write_mesh(inputMesh, argv[2]))
    {
        cerr << "Couldn't write output mesh to file " << argv[2] << endl;
        return -1;
    }
    return 0;
}

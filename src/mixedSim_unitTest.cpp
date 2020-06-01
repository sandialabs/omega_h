#include "MeshSim.h"
#include "SimDiscrete.h"
#include "SimMessages.h"
#include "SimError.h"
#include "SimErrorCodes.h"
#include "SimMeshingErrorCodes.h"
#include "SimDiscreteErrorCodes.h"
#include <iostream>
using namespace std;

#include "Omega_h_comm.hpp"
#include "Omega_h_mesh.hpp"
#include "Omega_h_file.hpp"
using namespace Omega_h;

void messageHandler(int type, const char *msg)
{
  switch (type) {
  case Sim_InfoMsg:
    cout<<"Info: "<<msg<<endl;
    break;
  case Sim_DebugMsg:
    cout<<"Debug: "<<msg<<endl;
    break;
  case Sim_WarningMsg:
    cout<<"Warning: "<<msg<<endl;
    break;
  case Sim_ErrorMsg:
    cout<<"Error: "<<msg<<endl;
    break;
  }
  return;
}

int main(int argc, char *argv[]) {
  auto lib = Library(&argc, &argv);
  auto comm = lib.world();
  auto mesh = Mesh(comm->library());
  try {
    Sim_logOn("importData1.log");  // start logging
    MS_init(); // Call before calling Sim_readLicenseFile
    Sim_readLicenseFile(0);

    pMesh meshtest; //mesh to load

/*
    // For tet on wedge, plus hex and pyramid adjacenct to wedge
    int numVerts = 12;
    double coords[12*3] = {0.000, 0.000, 0.000,
                           1.000, 0.000, 0.000,
                           0.500, 0.866, 0.000,
                           0.000, 0.000, 1.000, 
                           1.000, 0.000, 1.000,
                           0.500, 0.866, 1.000,
                           0.500, 0.289, 1.866,
                           1.000, 1.500, 0.500,
                           0.000, -1.000, 0.000,
                           1.000, -1.000, 0.000,
                           0.000, -1.000, 1.000,
                           1.000, -1.000, 1.000
                          };
    int numElems = 4;
    int elementType[4] = {12, 10, 11, 13};
    int elementData[6+4+5+8] = {0,1,2,3,4,5, 3,4,5,6, 5,4,1,2,7,
8,9,1,0,10,11,4,3};
    pVertex vReturn[12]; // array of created vertices
    pEntity eReturn[4]; // array of created entities
    //
*/
    //For pyramid on hex
    int numVerts = 9;
    double coords[9*3] =  {
                           0.000, 0.000, 0.000,
                           1.000, 0.000, 0.000,
                           1.000, 1.000, 0.000,
                           0.000, 1.000, 0.000,
                           0.000, 0.000, 1.000,
                           1.000, 0.000, 1.000,
                           1.000, 1.000, 1.000,
                           0.000, 1.000, 1.000,
                           0.500, 0.500, 2.000
                          };
    int numElems = 2;
    int elementType[2] = {13, 11};
    int elementData[8+5] = {0,1,2,3,4,5,6,7, 4,5,6,7,8};
    pVertex vReturn[9];
    pEntity eReturn[2];

/*
    //For tet on wedge
    int numVerts = 7;
    double coords[7*3] =  {0.000, 0.000, 0.000,
                           1.000, 0.000, 0.000,
                           0.500, 0.866, 0.000,
                           0.000, 0.000, 1.000,
                           1.000, 0.000, 1.000,
                           0.500, 0.866, 1.000,
                           0.500, 0.289, 1.866,
                          };
    int numElems = 2;
    int elementType[2] = {12, 10};
    int elementData[6+4] = {0,1,2,3,4,5,3,4,5,6};
    pVertex vReturn[7];
    pEntity eReturn[2];
*/

/*  For cube as per simmetrix example
    int numVerts = 8; // number of vertices  
    double coords[8*3] = {0.0,0.0,0.0,
                        1.0,0.0,0.0,
                        1.0,0.0,1.0,
                        0.0,0.0,1.0,
                        0.0,1.0,0.0,
                        1.0,1.0,0.0,
                        1.0,1.0,1.0,
                        0.0,1.0,1.0}; 
    int numElems = 12;
    int elementType[12] = {5,5,5,5,5,5,5,5,5,5,5,5}; 
    int elementData[12*3] =
{0,2,3,0,1,2,0,5,4,0,1,5,1,6,5,1,2,6,2,7,6,2,3,7,3,4,7,3,0,4,4,6,7,4,5,6}; 
    pVertex vReturn[8];
    pEntity eReturn[12];
*/
    SimDiscrete_start(0);
    Sim_setMessageHandler(messageHandler);
    pProgress progress = Progress_new();
    Progress_setDefaultCallback(progress);

    meshtest = M_new(0,0);
    if(M_importFromData(meshtest,numVerts,coords,numElems,
      elementType,elementData,vReturn,eReturn,progress)) {
      cerr<<"Error importing mesh data"<<endl;
      M_release(meshtest);
    }

    pDiscreteModel modeltest = DM_createFromMesh(meshtest, 0, progress);
    if(!modeltest) { //check for error
      cerr<<"Error creating Discrete model from mesh"<<endl;
      M_release(meshtest);
    }

    DM_findEdgesByFaceNormals(modeltest, 0, progress);
    DM_eliminateDanglingEdges(modeltest, progress);
    if(DM_completeTopology(modeltest, progress)) {
      cerr<<"Error completing Discrete model topology"<<endl;
      M_release(meshtest);
      GM_release(modeltest);
    }

    GM_write(modeltest,"/users/joshia5/simmodeler/Example_pyramid_hex.smd",0,progress);
    M_write(meshtest,"/users/joshia5/simmodeler/Example_pyramid_hex.sms", 0,progress);

    M_release(meshtest);
    GM_release(modeltest);
    Progress_delete(progress);
    SimDiscrete_stop(0);
    Sim_unregisterAllKeys();
    MS_exit();
    Sim_logOff();

  } catch (pSimError err) {
    cerr<<"SimModSuite error caught:"<<endl;
    cerr<<"  Error code: "<<SimError_code(err)<<endl;
    cerr<<"  Error string: "<<SimError_toString(err)<<endl;
    SimError_delete(err);
  } catch (...) {
    cerr<<"Unhandled exception caught"<<endl;
  }
  return 0;
}
